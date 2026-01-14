#include "SlamDump.h"
#include "SlamQuant.h"
#include "SlamCompat.h"
#include "SlamSolver.h"
#include "SlamQcOutput.h"
#include "SlamVarianceAnalysis.h"
#include "libem/slam_snp_em.h"
#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>

struct Args {
    std::string dumpPath;
    std::string bamPath;
    std::string gtfPath;
    std::string fastaPath;
    std::string outPrefix;
    std::string maskPath;
    bool snpMaskFromBam = false;
    std::string snpMaskOut;
    double snpMaskPval = 0.001;
    double snpMaskMinTcRatio = 0.3;
    uint32_t snpMaskMinCov = 6;
    uint32_t snpMaskMinAlt = 1;
    double snpMaskErr = -1.0;
    std::string autoTrimMode;
    std::string trimScope = "first";
    std::string strandness = "none";
    std::string qcReportPrefix;
    int trim5p = 0;
    int trim3p = 0;
    double errorRate = -1.0;
    double convRate = -1.0;
    uint32_t autoTrimMaxReads = 100000;
    uint32_t autoTrimMinReads = 1000;
    uint32_t autoTrimSmoothWindow = 5;
    uint32_t autoTrimSegMinLen = 3;
    uint32_t autoTrimMaxTrim = 15;
};

static bool parseArgs(int argc, char** argv, Args* args) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                exit(2);
            }
            return argv[++i];
        };
        if (a == "--dump") args->dumpPath = next("--dump");
        else if (a == "--bam") args->bamPath = next("--bam");
        else if (a == "--gtf") args->gtfPath = next("--gtf");
        else if (a == "--fasta") args->fastaPath = next("--fasta");
        else if (a == "--out") args->outPrefix = next("--out");
        else if (a == "--slamSnpMaskIn") args->maskPath = next("--slamSnpMaskIn");
        else if (a == "--snpMaskFromBam") args->snpMaskFromBam = true;
        else if (a == "--snpMaskOut") args->snpMaskOut = next("--snpMaskOut");
        else if (a == "--snpMaskPval") args->snpMaskPval = std::stod(next("--snpMaskPval"));
        else if (a == "--snpMaskMinTcRatio") args->snpMaskMinTcRatio = std::stod(next("--snpMaskMinTcRatio"));
        else if (a == "--snpMaskMinCov") args->snpMaskMinCov = static_cast<uint32_t>(std::stoul(next("--snpMaskMinCov")));
        else if (a == "--snpMaskMinAlt") args->snpMaskMinAlt = static_cast<uint32_t>(std::stoul(next("--snpMaskMinAlt")));
        else if (a == "--snpMaskErr") args->snpMaskErr = std::stod(next("--snpMaskErr"));
        else if (a == "--trim5p") args->trim5p = std::stoi(next("--trim5p"));
        else if (a == "--trim3p") args->trim3p = std::stoi(next("--trim3p"));
        else if (a == "--autoTrim") args->autoTrimMode = next("--autoTrim");
        else if (a == "--trimScope") args->trimScope = next("--trimScope");
        else if (a == "--strandness") args->strandness = next("--strandness");
        else if (a == "--slamQcReport") args->qcReportPrefix = next("--slamQcReport");
        else if (a == "--errorRate") args->errorRate = std::stod(next("--errorRate"));
        else if (a == "--convRate") args->convRate = std::stod(next("--convRate"));
        else if (a == "--autoTrimDetectionReads") args->autoTrimMaxReads = static_cast<uint32_t>(std::stoul(next("--autoTrimDetectionReads")));
        else if (a == "--autoTrimMinReads") args->autoTrimMinReads = static_cast<uint32_t>(std::stoul(next("--autoTrimMinReads")));
        else if (a == "--autoTrimSmoothWindow") args->autoTrimSmoothWindow = static_cast<uint32_t>(std::stoul(next("--autoTrimSmoothWindow")));
        else if (a == "--autoTrimSegMinLen") args->autoTrimSegMinLen = static_cast<uint32_t>(std::stoul(next("--autoTrimSegMinLen")));
        else if (a == "--autoTrimMaxTrim") args->autoTrimMaxTrim = static_cast<uint32_t>(std::stoul(next("--autoTrimMaxTrim")));
        else {
            std::cerr << "Unknown arg: " << a << "\n";
            return false;
        }
    }
    if (args->outPrefix.empty()) {
        std::cerr << "Usage: --dump <path> | --bam <path> --gtf <path> --fasta <path> "
                  << "--out <prefix> [--slamSnpMaskIn <bed.gz>] [--trim5p N --trim3p N]\n";
        return false;
    }
    if (args->dumpPath.empty() && args->bamPath.empty()) {
        std::cerr << "Either --dump or --bam is required.\n";
        return false;
    }
    if (!args->bamPath.empty() && (args->gtfPath.empty() || args->fastaPath.empty())) {
        std::cerr << "--bam requires --gtf and --fasta.\n";
        return false;
    }
    return true;
}

static int strandnessToInt(const std::string& s) {
    std::string x = s;
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x == "none" || x == "0") return 0;
    if (x == "sense" || x == "1") return 1;
    if (x == "antisense" || x == "2") return 2;
    return 0;
}

static void addMaskCount(std::unordered_map<uint64_t, uint32_t>& counts,
                         uint64_t pos,
                         bool isMismatch) {
    uint32_t& entry = counts[pos];
    uint32_t cov = entry >> 16;
    uint32_t mis = entry & 0xFFFF;
    if (cov < 0xFFFF) ++cov;
    if (isMismatch && mis < 0xFFFF) ++mis;
    entry = (cov << 16) | mis;
}

static std::unordered_set<uint64_t> buildMaskFromReads(const std::vector<SlamBufferedRead>& reads,
                                                       const Args& args,
                                                       double errorRate) {
    std::unordered_map<uint64_t, uint32_t> counts;
    counts.reserve(100000);
    for (const auto& r : reads) {
        for (const auto& p : r.positions) {
            bool isMis = (p.refBase != p.readBase);
            addMaskCount(counts, p.genomicPos, isMis);
        }
    }
    std::unordered_set<uint64_t> mask;
    const double pErr = (args.snpMaskErr > 0.0) ? args.snpMaskErr : errorRate;
    for (const auto& kv : counts) {
        uint32_t packed = kv.second;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        if (n < args.snpMaskMinCov) continue;
        if (k < args.snpMaskMinAlt) continue;
        double ratio = static_cast<double>(k) / static_cast<double>(n);
        if (ratio < args.snpMaskMinTcRatio) continue;
        double log_pval = log_binom_tail_cdf(n, k, pErr);
        double pval = std::exp(log_pval);
        if (pval < args.snpMaskPval) {
            mask.insert(kv.first);
        }
    }
    return mask;
}

static bool writeMaskBed(const std::string& path,
                         const std::unordered_set<uint64_t>& mask,
                         const SlamDumpMetadata& meta,
                         std::string* err) {
    if (path.empty()) return true;
    std::ofstream out(path.c_str());
    if (!out.good()) {
        if (err) *err = "Failed to write mask BED: " + path;
        return false;
    }
    // Build chromosome lookup by chrStart
    for (uint64_t pos : mask) {
        // find chr index
        auto it = std::upper_bound(meta.chrStart.begin(), meta.chrStart.end(), pos);
        size_t chrIdx = (it == meta.chrStart.begin()) ? 0 : static_cast<size_t>(std::distance(meta.chrStart.begin(), it) - 1);
        if (chrIdx >= meta.chrNames.size()) continue;
        uint64_t chrStart = meta.chrStart[chrIdx];
        uint64_t pos0 = pos - chrStart;
        out << meta.chrNames[chrIdx] << "\t" << pos0 << "\t" << (pos0 + 1) << "\n";
    }
    return true;
}

struct GeneInfo {
    std::string id;
    std::string name;
    char strand = '.';
    std::string chr;
    uint64_t start = 0;
    uint64_t end = 0;
};

struct ExonInterval {
    uint64_t start = 0;
    uint64_t end = 0;
    uint32_t geneIdx = 0;
};

static void parseGtfAttributes(const std::string& attrs, std::string* geneId, std::string* geneName) {
    // Very small parser for key "value"; pairs
    size_t pos = 0;
    while (pos < attrs.size()) {
        while (pos < attrs.size() && (attrs[pos] == ' ' || attrs[pos] == '\t' || attrs[pos] == ';')) pos++;
        size_t keyStart = pos;
        while (pos < attrs.size() && attrs[pos] != ' ' && attrs[pos] != '\t') pos++;
        if (pos == keyStart) break;
        std::string key = attrs.substr(keyStart, pos - keyStart);
        while (pos < attrs.size() && (attrs[pos] == ' ' || attrs[pos] == '\t')) pos++;
        if (pos >= attrs.size() || attrs[pos] != '"') {
            while (pos < attrs.size() && attrs[pos] != ';') pos++;
            continue;
        }
        ++pos; // skip quote
        size_t valStart = pos;
        while (pos < attrs.size() && attrs[pos] != '"') pos++;
        std::string val = attrs.substr(valStart, pos - valStart);
        while (pos < attrs.size() && attrs[pos] != ';') pos++;
        if (pos < attrs.size()) pos++;
        if (key == "gene_id") {
            *geneId = val;
        } else if (key == "gene_name") {
            *geneName = val;
        }
    }
}

static bool loadGtf(const std::string& path,
                    std::unordered_map<std::string, std::vector<ExonInterval>>& exonsByChr,
                    std::unordered_map<std::string, std::vector<ExonInterval>>& spansByChr,
                    std::vector<GeneInfo>& genes,
                    std::unordered_map<std::string, uint32_t>& geneIndex,
                    std::string* err) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        if (err) *err = "Failed to open GTF: " + path;
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cols;
        size_t start = 0;
        for (size_t i = 0; i < 8; ++i) {
            size_t tab = line.find('\t', start);
            if (tab == std::string::npos) break;
            cols.push_back(line.substr(start, tab - start));
            start = tab + 1;
        }
        cols.push_back(line.substr(start));
        if (cols.size() < 9) continue;
        if (cols[2] != "exon") continue;
        const std::string& chr = cols[0];
        uint64_t exonStart = std::stoull(cols[3]) - 1;
        uint64_t exonEnd = std::stoull(cols[4]);
        char strand = cols[6].empty() ? '.' : cols[6][0];
        std::string geneId, geneName;
        parseGtfAttributes(cols[8], &geneId, &geneName);
        if (geneId.empty()) continue;
        auto it = geneIndex.find(geneId);
        uint32_t gidx;
        if (it == geneIndex.end()) {
            gidx = static_cast<uint32_t>(genes.size());
            geneIndex[geneId] = gidx;
            GeneInfo gi;
            gi.id = geneId;
            gi.name = geneName.empty() ? geneId : geneName;
            gi.strand = strand;
            gi.chr = chr;
            gi.start = exonStart;
            gi.end = exonEnd;
            genes.push_back(std::move(gi));
        } else {
            gidx = it->second;
            if (genes[gidx].chr.empty()) {
                genes[gidx].chr = chr;
            }
            genes[gidx].start = std::min(genes[gidx].start, exonStart);
            genes[gidx].end = std::max(genes[gidx].end, exonEnd);
        }
        ExonInterval ex;
        ex.start = exonStart;
        ex.end = exonEnd;
        ex.geneIdx = gidx;
        exonsByChr[chr].push_back(ex);
    }
    // Build gene span intervals by chromosome
    for (uint32_t i = 0; i < genes.size(); ++i) {
        const GeneInfo& gi = genes[i];
        if (gi.chr.empty()) continue;
        ExonInterval span;
        span.start = gi.start;
        span.end = gi.end;
        span.geneIdx = i;
        spansByChr[gi.chr].push_back(span);
    }
    // Sort intervals
    for (auto& kv : exonsByChr) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end(), [](const ExonInterval& a, const ExonInterval& b) {
            return a.start < b.start;
        });
    }
    for (auto& kv : spansByChr) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end(), [](const ExonInterval& a, const ExonInterval& b) {
            return a.start < b.start;
        });
    }
    return true;
}

static void collectOverlaps(const std::vector<ExonInterval>& intervals,
                            uint64_t start, uint64_t end,
                            std::set<uint32_t>& outGenes) {
    // intervals sorted by start
    size_t lo = 0, hi = intervals.size();
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (intervals[mid].start < start) lo = mid + 1;
        else hi = mid;
    }
    size_t idx = (lo > 0) ? lo - 1 : 0;
    for (size_t i = idx; i < intervals.size(); ++i) {
        const auto& iv = intervals[i];
        if (iv.start > end) break;
        if (iv.end > start && iv.start < end) {
            outGenes.insert(iv.geneIdx);
        }
    }
}

static std::vector<SlamBufferedRead> buildReadsFromBam(const Args& args,
                                                       SlamDumpMetadata* meta,
                                                       std::string* err) {
    std::unordered_map<std::string, std::vector<ExonInterval>> exonsByChr;
    std::unordered_map<std::string, std::vector<ExonInterval>> spansByChr;
    std::vector<GeneInfo> genes;
    std::unordered_map<std::string, uint32_t> geneIndex;
    if (!loadGtf(args.gtfPath, exonsByChr, spansByChr, genes, geneIndex, err)) {
        return {};
    }
    meta->geneIds.clear();
    meta->geneNames.clear();
    for (const auto& gi : genes) {
        meta->geneIds.push_back(gi.id);
        meta->geneNames.push_back(gi.name);
    }

    faidx_t* fai = fai_load(args.fastaPath.c_str());
    if (!fai) {
        if (err) *err = "Failed to open FASTA index: " + args.fastaPath;
        return {};
    }

    samFile* fp = sam_open(args.bamPath.c_str(), "r");
    if (!fp) {
        if (err) *err = "Failed to open BAM: " + args.bamPath;
        fai_destroy(fai);
        return {};
    }
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        if (err) *err = "Failed to read BAM header";
        sam_close(fp);
        fai_destroy(fai);
        return {};
    }
    meta->chrNames.clear();
    meta->chrStart.clear();
    meta->chrStart.reserve(hdr->n_targets + 1);
    uint64_t running = 0;
    for (int i = 0; i < hdr->n_targets; ++i) {
        meta->chrNames.push_back(hdr->target_name[i]);
        meta->chrStart.push_back(running);
        running += hdr->target_len[i];
    }

    std::vector<SlamBufferedRead> reads;
    bam1_t* b = bam_init1();
    while (sam_read1(fp, hdr, b) >= 0) {
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if (b->core.tid < 0 || b->core.tid >= hdr->n_targets) continue;
        const char* chr = hdr->target_name[b->core.tid];
        auto exIt = exonsByChr.find(chr);
        auto spIt = spansByChr.find(chr);
        if (exIt == exonsByChr.end() || spIt == spansByChr.end()) continue;

        uint32_t* cigar = bam_get_cigar(b);
        int nCigar = b->core.n_cigar;
        int64_t refPos = b->core.pos; // 0-based
        int32_t readPos = 0;

        // compute ref span for fasta fetch
        int64_t refEnd = refPos;
        for (int i = 0; i < nCigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP) {
                refEnd += len;
            }
        }
        int seqLen = 0;
        char* refSeq = faidx_fetch_seq(fai, chr, static_cast<int>(refPos), static_cast<int>(refEnd - 1), &seqLen);
        if (!refSeq || seqLen <= 0) {
            if (refSeq) free(refSeq);
            continue;
        }

        std::set<uint32_t> exonicGenes;
        std::set<uint32_t> spanGenes;
        std::vector<SlamBufferedPosition> positions;
        positions.reserve(b->core.l_qseq);

        for (int i = 0; i < nCigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (int j = 0; j < len; ++j) {
                    int64_t gpos = refPos + j;
                    int32_t rpos = readPos + j;
                    if (rpos < 0 || rpos >= b->core.l_qseq) continue;
                    int refOff = static_cast<int>(gpos - refPos);
                    if (refOff < 0 || refOff >= seqLen) continue;
                    char refBase = refSeq[refOff];
                    uint8_t readBase = bam_seqi(bam_get_seq(b), rpos);
                    uint8_t qual = bam_get_qual(b)[rpos];
                    SlamBufferedPosition bp;
                    bp.readPos = static_cast<uint32_t>(rpos);
                    bp.genomicPos = meta->chrStart[b->core.tid] + static_cast<uint64_t>(gpos);
                    // Map bases to 0-3
                    auto mapBase = [](char c) -> uint8_t {
                        switch (std::toupper(c)) {
                            case 'A': return 0;
                            case 'C': return 1;
                            case 'G': return 2;
                            case 'T': return 3;
                            default: return 4;
                        }
                    };
                    bp.refBase = mapBase(refBase);
                    char rb = (readBase <= 15) ? seq_nt16_str[readBase] : 'N';
                    bp.readBase = mapBase(rb);
                    bp.qual = (qual == 255) ? 30 : qual;
                    bp.secondMate = (b->core.flag & BAM_FREAD2) != 0;
                    bp.overlap = false;
                    if (bp.refBase < 4 && bp.readBase < 4) {
                        positions.push_back(bp);
                    }
                }
                // overlap gene intervals for this block
                collectOverlaps(exIt->second, refPos, refPos + len, exonicGenes);
                collectOverlaps(spIt->second, refPos, refPos + len, spanGenes);
                refPos += len;
                readPos += len;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                collectOverlaps(spIt->second, refPos, refPos + len, spanGenes);
                refPos += len;
            } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
                readPos += len;
            } else if (op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                // no-op
            }
        }
        free(refSeq);

        std::set<uint32_t> geneIds = exonicGenes.empty() ? spanGenes : exonicGenes;
        if (geneIds.empty()) {
            continue;
        }
        bool isIntronic = exonicGenes.empty();

        bool readMinus = (b->core.flag & BAM_FREVERSE) != 0;
        bool anySense = false;
        bool anyOpp = false;
        for (uint32_t gid : geneIds) {
            if (gid >= genes.size()) continue;
            char gstrand = genes[gid].strand;
            if (gstrand == '.') continue;
            bool sense = (gstrand == (readMinus ? '-' : '+'));
            anySense = anySense || sense;
            anyOpp = anyOpp || !sense;
            if (anySense && anyOpp) break;
        }
        bool oppositeStrand = anySense ? false : anyOpp;

        SlamBufferedRead r;
        std::string qname = bam_get_qname(b);
        r.readName = qname;
        r.readLength0 = static_cast<uint32_t>(b->core.l_qseq);
        r.readLength1 = 0;
        r.isMinus = readMinus;
        r.oppositeStrand = oppositeStrand;
        r.isIntronic = isIntronic;
        r.fileIndex = 0;
        r.weight = geneIds.empty() ? 0.0 : (1.0 / static_cast<double>(geneIds.size()));
        r.geneIds.assign(geneIds.begin(), geneIds.end());
        r.positions = std::move(positions);
        reads.push_back(std::move(r));
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    fai_destroy(fai);
    return reads;
}

static void writeSlamOut(const std::string& outFile,
                         const std::vector<std::string>& geneIds,
                         const std::vector<std::string>& geneNames,
                         const SlamQuant& quant,
                         double errorRate,
                         double convRate) {
    std::ofstream out(outFile.c_str());
    if (!out.good()) return;
    out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tMAP\tSigma\tLogLikelihood\n";
    SlamSolver solver(errorRate, convRate);
    const auto& genes = quant.genes();
    for (size_t i = 0; i < genes.size(); ++i) {
        const SlamGeneStats& stats = genes[i];
        if (stats.readCount <= 0.0) continue;
        SlamResult res = solver.solve(stats.histogram);
        const std::string& gid = (i < geneIds.size()) ? geneIds[i] : std::string("GENE_") + std::to_string(i);
        const std::string& gname = (i < geneNames.size() && !geneNames[i].empty()) ? geneNames[i] : gid;
        out << gid << "\t" << gname << "\t"
            << stats.readCount << "\t"
            << stats.conversions << "\t"
            << stats.coverage << "\t"
            << res.ntr << "\t"
            << res.ntr << "\t"
            << res.sigma << "\t"
            << res.log_likelihood << "\n";
    }
}

static SlamVarianceTrimResult computeTrimForReads(const std::vector<SlamBufferedRead>& reads,
                                                  uint32_t maxReads,
                                                  uint32_t minReads,
                                                  uint32_t smoothWindow,
                                                  uint32_t minSegLen,
                                                  uint32_t maxTrim,
                                                  uint32_t readLength) {
    SlamVarianceAnalyzer analyzer(maxReads, minReads, smoothWindow, minSegLen, maxTrim);
    for (const auto& r : reads) {
        if (!analyzer.recordRead()) break;
        for (const auto& p : r.positions) {
            bool isT = false;
            bool isTc = false;
            if (!r.isMinus) {
                isT = (p.refBase == 3);
                isTc = (p.refBase == 3 && p.readBase == 1);
            } else {
                isT = (p.refBase == 0);
                isTc = (p.refBase == 0 && p.readBase == 2);
            }
            analyzer.recordPosition(p.readPos, p.qual, isT, isTc);
        }
    }
    return analyzer.computeTrim(readLength);
}

int main(int argc, char** argv) {
    Args args;
    if (!parseArgs(argc, argv, &args)) {
        return 2;
    }

    SlamDumpMetadata meta;
    std::vector<SlamBufferedRead> reads;
    std::string err;
    if (!args.dumpPath.empty()) {
        if (!readSlamDump(args.dumpPath, &meta, &reads, &err)) {
            std::cerr << "Failed to read dump: " << err << "\n";
            return 1;
        }
    } else {
        reads = buildReadsFromBam(args, &meta, &err);
        if (!err.empty()) {
            std::cerr << "Failed to build reads from BAM: " << err << "\n";
            return 1;
        }
    }
    if (meta.geneIds.empty()) {
        std::cerr << "Dump contains no genes\n";
        return 1;
    }

    double errorRate = (args.errorRate >= 0.0) ? args.errorRate : meta.errorRate;
    double convRate = (args.convRate >= 0.0) ? args.convRate : meta.convRate;
    if (errorRate <= 0.0) errorRate = 0.001;

    SlamSnpMask mask;
    SlamSnpMask* maskPtr = nullptr;
    if (!args.maskPath.empty()) {
        if (!mask.loadBedWithChrMap(args.maskPath, meta.chrNames, meta.chrStart, &err)) {
            std::cerr << "Failed to load mask: " << err << "\n";
            return 1;
        }
        maskPtr = &mask;
    } else if (args.snpMaskFromBam) {
        auto memMask = buildMaskFromReads(reads, args, errorRate);
        mask.loadFromPositions(memMask);
        maskPtr = &mask;
        writeMaskBed(args.snpMaskOut.empty() ? (args.outPrefix + "snp_mask.bed") : args.snpMaskOut,
                     memMask, meta, &err);
    }

    // Determine read length (for trim computation)
    uint32_t readLength = 0;
    for (const auto& r : reads) {
        readLength = std::max(readLength, static_cast<uint32_t>(r.readLength0 + r.readLength1));
    }
    if (readLength == 0) readLength = 100;

    // Compute trims if requested
    std::unordered_map<uint32_t, SlamVarianceTrimResult> perFileTrim;
    SlamVarianceTrimResult globalTrim;
    if (args.autoTrimMode == "variance") {
        if (args.trimScope == "per-file") {
            std::unordered_map<uint32_t, std::vector<SlamBufferedRead>> byFile;
            for (const auto& r : reads) {
                byFile[r.fileIndex].push_back(r);
            }
            for (const auto& kv : byFile) {
                perFileTrim[kv.first] = computeTrimForReads(kv.second,
                                                           args.autoTrimMaxReads,
                                                           args.autoTrimMinReads,
                                                           args.autoTrimSmoothWindow,
                                                           args.autoTrimSegMinLen,
                                                           args.autoTrimMaxTrim,
                                                           readLength);
            }
        } else {
            globalTrim = computeTrimForReads(reads,
                                             args.autoTrimMaxReads,
                                             args.autoTrimMinReads,
                                             args.autoTrimSmoothWindow,
                                             args.autoTrimSegMinLen,
                                             args.autoTrimMaxTrim,
                                             readLength);
            if (globalTrim.success) {
                args.trim5p = globalTrim.trim5p;
                args.trim3p = globalTrim.trim3p;
            }
        }
    }

    auto runPartition = [&](const std::vector<SlamBufferedRead>& partReads,
                            int trim5p, int trim3p) {
        SlamQuant q(meta.geneIds.size(), false);
        q.enableReadBuffer(static_cast<uint64_t>(partReads.size()));
        for (const auto& r : partReads) {
            q.bufferRead(SlamBufferedRead(r));
        }
        SlamCompatConfig cfg;
        cfg.trim5p = trim5p;
        cfg.trim3p = trim3p;
        SlamCompat compat(cfg, {}, {});
        q.replayBufferedReads(&compat, maskPtr, strandnessToInt(args.strandness));
        return q;
    };

    SlamQuant merged(meta.geneIds.size(), false);
    if (args.trimScope == "per-file") {
        std::unordered_map<uint32_t, std::vector<SlamBufferedRead>> byFile;
        for (const auto& r : reads) {
            byFile[r.fileIndex].push_back(r);
        }
        for (const auto& kv : byFile) {
            int t5 = args.trim5p;
            int t3 = args.trim3p;
            auto it = perFileTrim.find(kv.first);
            if (it != perFileTrim.end() && it->second.success) {
                t5 = it->second.trim5p;
                t3 = it->second.trim3p;
            }
            SlamQuant q = runPartition(kv.second, t5, t3);
            merged.merge(q);
        }
    } else {
        SlamQuant q = runPartition(reads, args.trim5p, args.trim3p);
        merged.merge(q);
    }

    std::string outBase = args.outPrefix;
    writeSlamOut(outBase + "SlamQuant.out", meta.geneIds, meta.geneNames, merged, errorRate, convRate);
    merged.writeDiagnostics(outBase + "SlamQuant.out.diagnostics");
    merged.writeTransitions(outBase + "SlamQuant.out.transitions.tsv");
    merged.writeMismatches(outBase + "SlamQuant.out.mismatches.tsv", outBase);
    merged.writeMismatchDetails(outBase + "SlamQuant.out.mismatchdetails.tsv");

    if (!args.qcReportPrefix.empty()) {
        std::string jsonPath = args.qcReportPrefix + ".slam_qc.json";
        std::string htmlPath = args.qcReportPrefix + ".slam_qc.html";
        SlamVarianceTrimResult* trimPtr = nullptr;
        if (args.autoTrimMode == "variance" && globalTrim.success) {
            trimPtr = &globalTrim;
        }
        writeSlamQcComprehensiveJson(merged, jsonPath, args.trim5p, args.trim3p, trimPtr);
        writeSlamQcComprehensiveHtml(jsonPath, htmlPath);
    }

    std::cerr << "Requant complete: " << outBase << "SlamQuant.out\n";
    return 0;
}
