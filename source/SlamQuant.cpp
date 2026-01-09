#include "SlamQuant.h"

#include "Genome.h"
#include "Transcriptome.h"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {
constexpr uint32_t kSnpMinCoverage = 10;
constexpr double kSnpMinMismatchFrac = 0.5;
const char* debugDropReasonName(SlamDebugDropReason reason) {
    switch (reason) {
        case SlamDebugDropReason::None:
            return "PASS";
        case SlamDebugDropReason::NoGenes:
            return "DROP_NO_GENES";
        case SlamDebugDropReason::SnpMask:
            return "DROP_SNP_MASK";
        case SlamDebugDropReason::Strandness:
            return "DROP_STRANDNESS";
        case SlamDebugDropReason::ZeroWeight:
            return "DROP_ZERO_WEIGHT";
        default:
            return "UNKNOWN";
    }
}
}

SlamQuant::SlamQuant(uint32_t nGenes, bool snpDetect)
    : geneStats_(nGenes), snpDetectEnabled_(snpDetect) {}

SlamQuant::SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes, bool snpDetect)
    : geneStats_(nGenes), snpDetectEnabled_(snpDetect),
      allowedGenes_(std::move(allowedGenes)) {}

void SlamQuant::initDebug(const Transcriptome& tr,
                          const std::unordered_set<std::string>& debugGenes,
                          const std::unordered_set<std::string>& debugReads,
                          size_t maxReads,
                          const std::string& outPrefix) {
    debugEnabled_ = !debugGenes.empty() || !debugReads.empty();
    debugMaxReads_ = maxReads;
    debugOutPrefix_ = outPrefix;
    debugReadRecords_.clear();
    debugReadSet_.clear();
    debugGeneMask_.clear();
    debugGeneStats_.clear();
    if (!debugEnabled_) {
        return;
    }
    debugReadSet_ = debugReads;
    if (!debugGenes.empty()) {
        debugGeneMask_.assign(tr.nGe, 0);
        debugGeneStats_.assign(tr.nGe, SlamDebugGeneStats());
        std::unordered_set<uint32_t> numericGenes;
        for (const auto& name : debugGenes) {
            if (name.empty()) {
                continue;
            }
            bool numeric = true;
            for (char c : name) {
                if (c < '0' || c > '9') {
                    numeric = false;
                    break;
                }
            }
            if (numeric) {
                uint32_t gid = static_cast<uint32_t>(std::stoul(name));
                if (gid < tr.nGe) {
                    numericGenes.insert(gid);
                }
            }
        }
        for (uint32_t i = 0; i < tr.nGe; ++i) {
            if (numericGenes.count(i) > 0) {
                debugGeneMask_[i] = 1;
                continue;
            }
            const std::string& gid = tr.geID[i];
            const std::string& gname = tr.geName[i];
            const std::string& gcanon = (i < tr.geIDCanonical.size()) ? tr.geIDCanonical[i] : std::string();
            if ((!gid.empty() && debugGenes.count(gid) > 0) ||
                (!gname.empty() && debugGenes.count(gname) > 0) ||
                (!gcanon.empty() && debugGenes.count(gcanon) > 0)) {
                debugGeneMask_[i] = 1;
            }
        }
    }
}

bool SlamQuant::debugGeneEnabled(uint32_t geneId) const {
    return debugEnabled_ && geneId < debugGeneMask_.size() && debugGeneMask_[geneId];
}

bool SlamQuant::debugReadMatch(const char* readName) const {
    if (!debugEnabled_ || debugReadSet_.empty() || readName == nullptr) {
        return false;
    }
    std::string name(readName);
    size_t end = name.find_first_of(" \t");
    if (end != std::string::npos) {
        name = name.substr(0, end);
    }
    if (!name.empty() && name[0] == '@') {
        name.erase(0, 1);
    }
    return debugReadSet_.count(name) > 0;
}

void SlamQuant::debugCountDrop(uint32_t geneId, SlamDebugDropReason reason) {
    if (!debugGeneEnabled(geneId)) {
        return;
    }
    SlamDebugGeneStats& stats = debugGeneStats_[geneId];
    if (reason == SlamDebugDropReason::SnpMask) {
        stats.dropsSnpMask++;
    } else if (reason == SlamDebugDropReason::Strandness) {
        stats.dropsStrandness++;
    }
}

void SlamQuant::debugAddAssignment(uint32_t geneId, double weight, bool intronic,
                                   bool oppositeStrand, uint16_t nT, uint8_t k) {
    if (!debugGeneEnabled(geneId)) {
        return;
    }
    SlamDebugGeneStats& stats = debugGeneStats_[geneId];
    stats.readsAssigned++;
    if (intronic) {
        stats.readsAssignedIntronic++;
        stats.intronicWeight += weight;
    } else {
        stats.exonicWeight += weight;
        stats.coverage += static_cast<double>(nT) * weight;
        stats.conversions += static_cast<double>(k) * weight;
    }
    if (oppositeStrand) {
        stats.antisenseWeight += weight;
    } else {
        stats.senseWeight += weight;
    }
}

void SlamQuant::debugLogRead(const SlamDebugReadRecord& record) {
    if (!debugEnabled_ || debugMaxReads_ == 0) {
        return;
    }
    if (debugReadRecords_.size() >= debugMaxReads_) {
        return;
    }
    debugReadRecords_.push_back(record);
}

const char* slamMismatchCategoryName(SlamMismatchCategory cat) {
    switch (cat) {
        case SlamMismatchCategory::Exonic:
            return "Exonic";
        case SlamMismatchCategory::ExonicSense:
            return "ExonicSense";
        case SlamMismatchCategory::Intronic:
            return "Intronic";
        case SlamMismatchCategory::IntronicSense:
            return "IntronicSense";
        default:
            return "Unknown";
    }
}

void SlamQuant::addRead(uint32_t geneId, uint16_t nT, uint8_t k, double weight) {
    if (geneId >= geneStats_.size() || weight <= 0.0) {
        return;
    }
    if (!allowedGenes_.empty() && geneId < allowedGenes_.size() && allowedGenes_[geneId] == 0) {
        return;
    }
    if (nT > 511) nT = 511;
    if (k > 255) k = 255;
    uint16_t key = static_cast<uint16_t>((nT << 8) | k);
    SlamGeneStats& stats = geneStats_[geneId];
    stats.histogram[key] += weight;
    stats.readCount += weight;
    stats.conversions += weight * static_cast<double>(k);
    stats.coverage += weight * static_cast<double>(nT);
}

void SlamQuant::recordSnpObservation(uint64_t pos, bool isMismatch) {
    if (!snpDetectEnabled_) {
        return;
    }
    uint32_t &entry = snpMask_[pos];
    uint32_t cov = entry >> 16;
    uint32_t mis = entry & 0xFFFF;
    if (cov < 0xFFFF) {
        ++cov;
    }
    if (isMismatch && mis < 0xFFFF) {
        ++mis;
    }
    entry = (cov << 16) | mis;
}

void SlamQuant::bufferSnpRead(uint32_t geneId, uint16_t nT,
                              const std::vector<uint32_t>& mismatchPositions, double weight) {
    if (!snpDetectEnabled_ || mismatchPositions.empty() || weight <= 0.0) {
        return;
    }
    if (!allowedGenes_.empty() && geneId < allowedGenes_.size() && allowedGenes_[geneId] == 0) {
        return;
    }
    uint16_t k = static_cast<uint16_t>(mismatchPositions.size() > 0xFFFF
                                           ? 0xFFFF
                                           : mismatchPositions.size());
    snpReadBuffer_.push_back(geneId);
    snpReadBuffer_.push_back((static_cast<uint32_t>(k) << 16) | nT);
    for (uint16_t i = 0; i < k; ++i) {
        snpReadBuffer_.push_back(mismatchPositions[i]);
    }
    snpReadWeights_.push_back(weight);
}

void SlamQuant::finalizeSnpMask(SlamSnpBufferStats* outStats) {
    if (outStats) {
        *outStats = SlamSnpBufferStats();
        outStats->bufferedReads = snpReadWeights_.size();
        outStats->bufferEntries = snpReadBuffer_.size();
        outStats->maskEntries = snpMask_.size();
        outStats->bufferBytes =
            static_cast<uint64_t>(snpReadBuffer_.size()) * sizeof(uint32_t) +
            static_cast<uint64_t>(snpReadWeights_.size()) * sizeof(double);
    }
    if (!snpDetectEnabled_ || snpFinalized_) {
        return;
    }
    snpFinalized_ = true;
    if (snpReadBuffer_.empty()) {
        return;
    }
    std::unordered_set<uint64_t> blacklist;
    blacklist.reserve(snpMask_.size());
    for (const auto &kv : snpMask_) {
        uint32_t cov = kv.second >> 16;
        uint32_t mis = kv.second & 0xFFFF;
        if (cov > kSnpMinCoverage && cov > 0) {
            double rate = static_cast<double>(mis) / static_cast<double>(cov);
            if (rate > kSnpMinMismatchFrac) {
                blacklist.insert(kv.first);
            }
        }
    }
    if (outStats) {
        outStats->blacklistEntries = blacklist.size();
    }
    uint64_t totalMismatches = 0;
    uint64_t totalMismatchesKept = 0;
    size_t idx = 0;
    size_t readIdx = 0;
    while (idx + 1 < snpReadBuffer_.size() && readIdx < snpReadWeights_.size()) {
        uint32_t geneId = snpReadBuffer_[idx++];
        uint32_t header = snpReadBuffer_[idx++];
        uint16_t nT = static_cast<uint16_t>(header & 0xFFFF);
        uint16_t k = static_cast<uint16_t>(header >> 16);
        uint16_t validK = 0;
        totalMismatches += k;
        for (uint16_t i = 0; i < k && idx < snpReadBuffer_.size(); ++i, ++idx) {
            uint32_t pos = snpReadBuffer_[idx];
            if (blacklist.find(pos) == blacklist.end()) {
                ++validK;
            }
        }
        totalMismatchesKept += validK;
        uint8_t k8 = static_cast<uint8_t>(validK > 255 ? 255 : validK);
        addRead(geneId, nT, k8, snpReadWeights_[readIdx]);
        ++readIdx;
    }
    if (outStats) {
        outStats->bufferedMismatches = totalMismatches;
        outStats->bufferedMismatchesKept = totalMismatchesKept;
        uint64_t reads = snpReadWeights_.size();
        if (reads > 0) {
            outStats->avgMismatches = static_cast<double>(totalMismatches) /
                static_cast<double>(reads);
            outStats->avgMismatchesKept = static_cast<double>(totalMismatchesKept) /
                static_cast<double>(reads);
        }
    }
    snpReadBuffer_.clear();
    snpReadWeights_.clear();
    snpMask_.clear();
}

void SlamQuant::addTransitionBase(SlamMismatchCategory category, uint32_t readPos, bool secondMate,
                                  bool overlap, bool opposite, int genomicBase, int readBase, double weight) {
    if (weight <= 0.0 || genomicBase < 0 || genomicBase > 3 || readBase < 0 || readBase > 3) {
        return;
    }
    const size_t cat = static_cast<size_t>(category);
    transitions_[cat].coverage[genomicBase] += weight;
    if (genomicBase != readBase) {
        transitions_[cat].mismatches[genomicBase][readBase] += weight;
    }
    if (!overlap) {
        SlamTransitionStats& orient = secondMate ? transitionsSecond_[cat] : transitionsFirst_[cat];
        orient.coverage[genomicBase] += weight;
        if (genomicBase != readBase) {
            orient.mismatches[genomicBase][readBase] += weight;
        }
    }
    if (readPos >= positionTransitions_[cat].size()) {
        positionTransitions_[cat].resize(readPos + 1);
    }
    SlamPositionStats& pos = positionTransitions_[cat][readPos];
    const size_t ov = overlap ? 1 : 0;
    const size_t opp = opposite ? 1 : 0;
    pos.coverage[ov][opp][genomicBase] += weight;
    if (genomicBase != readBase) {
        pos.mismatches[ov][opp][genomicBase][readBase] += weight;
    }
}

void SlamQuant::merge(const SlamQuant& other) {
    if (other.geneStats_.size() != geneStats_.size()) {
        return;
    }
    if (other.snpDetectEnabled_) {
        snpDetectEnabled_ = true;
    }
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        SlamGeneStats& dst = geneStats_[i];
        const SlamGeneStats& src = other.geneStats_[i];
        dst.readCount += src.readCount;
        dst.conversions += src.conversions;
        dst.coverage += src.coverage;
        for (const auto& kv : src.histogram) {
            dst.histogram[kv.first] += kv.second;
        }
    }
    // Merge diagnostics
    diag_.readsDroppedSnpMask += other.diag_.readsDroppedSnpMask;
    diag_.readsDroppedStrandness += other.diag_.readsDroppedStrandness;
    diag_.readsZeroGenes += other.diag_.readsZeroGenes;
    diag_.readsProcessed += other.diag_.readsProcessed;
    diag_.readsNAlignWithGeneZero += other.diag_.readsNAlignWithGeneZero;
    diag_.compatAlignsReclassifiedIntronic += other.diag_.compatAlignsReclassifiedIntronic;
    diag_.compatAlignsLenientAccepted += other.diag_.compatAlignsLenientAccepted;
    diag_.compatAlignsOverlapWeightApplied += other.diag_.compatAlignsOverlapWeightApplied;
    diag_.compatPositionsSkippedOverlap += other.diag_.compatPositionsSkippedOverlap;
    diag_.compatPositionsSkippedTrim += other.diag_.compatPositionsSkippedTrim;
    diag_.readsSumWeightLessThanOne += other.diag_.readsSumWeightLessThanOne;
    for (const auto& kv : other.diag_.nTrDistribution) {
        diag_.nTrDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.geneSetSizeDistribution) {
        diag_.geneSetSizeDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.nAlignWithGeneDistribution) {
        diag_.nAlignWithGeneDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.sumWeightDistribution) {
        diag_.sumWeightDistribution[kv.first] += kv.second;
    }
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        for (int g = 0; g < 4; ++g) {
            transitions_[c].coverage[g] += other.transitions_[c].coverage[g];
            transitionsFirst_[c].coverage[g] += other.transitionsFirst_[c].coverage[g];
            transitionsSecond_[c].coverage[g] += other.transitionsSecond_[c].coverage[g];
            for (int r = 0; r < 4; ++r) {
                transitions_[c].mismatches[g][r] += other.transitions_[c].mismatches[g][r];
                transitionsFirst_[c].mismatches[g][r] += other.transitionsFirst_[c].mismatches[g][r];
                transitionsSecond_[c].mismatches[g][r] += other.transitionsSecond_[c].mismatches[g][r];
            }
        }
        if (other.positionTransitions_[c].size() > positionTransitions_[c].size()) {
            positionTransitions_[c].resize(other.positionTransitions_[c].size());
        }
        for (size_t i = 0; i < other.positionTransitions_[c].size(); ++i) {
            for (int ov = 0; ov < 2; ++ov) {
                for (int opp = 0; opp < 2; ++opp) {
                    for (int g = 0; g < 4; ++g) {
                        positionTransitions_[c][i].coverage[ov][opp][g] +=
                            other.positionTransitions_[c][i].coverage[ov][opp][g];
                        for (int r = 0; r < 4; ++r) {
                            positionTransitions_[c][i].mismatches[ov][opp][g][r] +=
                                other.positionTransitions_[c][i].mismatches[ov][opp][g][r];
                        }
                    }
                }
            }
        }
    }
    if (!other.snpMask_.empty()) {
        for (const auto &kv : other.snpMask_) {
            uint32_t &entry = snpMask_[kv.first];
            uint32_t cov = entry >> 16;
            uint32_t mis = entry & 0xFFFF;
            uint32_t addCov = kv.second >> 16;
            uint32_t addMis = kv.second & 0xFFFF;
            cov = std::min<uint32_t>(0xFFFF, cov + addCov);
            mis = std::min<uint32_t>(0xFFFF, mis + addMis);
            entry = (cov << 16) | mis;
        }
    }
    if (!other.snpReadBuffer_.empty()) {
        snpReadBuffer_.insert(snpReadBuffer_.end(), other.snpReadBuffer_.begin(), other.snpReadBuffer_.end());
        snpReadWeights_.insert(snpReadWeights_.end(), other.snpReadWeights_.begin(), other.snpReadWeights_.end());
    }
    if (debugEnabled_ && other.debugEnabled_) {
        if (debugGeneStats_.size() == other.debugGeneStats_.size()) {
            for (size_t i = 0; i < debugGeneStats_.size(); ++i) {
                SlamDebugGeneStats& dst = debugGeneStats_[i];
                const SlamDebugGeneStats& src = other.debugGeneStats_[i];
                dst.exonicWeight += src.exonicWeight;
                dst.intronicWeight += src.intronicWeight;
                dst.senseWeight += src.senseWeight;
                dst.antisenseWeight += src.antisenseWeight;
                dst.conversions += src.conversions;
                dst.coverage += src.coverage;
                dst.readsAssigned += src.readsAssigned;
                dst.readsAssignedIntronic += src.readsAssignedIntronic;
                dst.dropsSnpMask += src.dropsSnpMask;
                dst.dropsStrandness += src.dropsStrandness;
            }
        }
        if (debugMaxReads_ != 0 && !other.debugReadRecords_.empty()) {
            for (const auto& rec : other.debugReadRecords_) {
                if (debugReadRecords_.size() >= debugMaxReads_) {
                    break;
                }
                debugReadRecords_.push_back(rec);
            }
        }
    }
}

void SlamQuant::write(const Transcriptome& tr, const std::string& outFile,
                      double errorRate, double convRate) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tMAP\tSigma\tLogLikelihood\n";

    SlamSolver solver(errorRate, convRate);
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        const SlamGeneStats& stats = geneStats_[i];
        if (stats.readCount <= 0.0) {
            continue;
        }
        SlamResult res = solver.solve(stats.histogram);
        const std::string& geneId = tr.geID[i];
        const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
        out << geneId << "\t"
            << geneName << "\t"
            << stats.readCount << "\t"
            << stats.conversions << "\t"
            << stats.coverage << "\t"
            << res.ntr << "\t"
            << res.ntr << "\t"
            << res.sigma << "\t"
            << res.log_likelihood << "\n";
    }
}

void SlamQuant::writeDiagnostics(const std::string& diagFile) const {
    std::ofstream out(diagFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Metric\tValue\n";
    out << "readsProcessed\t" << diag_.readsProcessed << "\n";
    out << "readsDroppedSnpMask\t" << diag_.readsDroppedSnpMask << "\n";
    out << "readsDroppedStrandness\t" << diag_.readsDroppedStrandness << "\n";
    out << "readsZeroGenes\t" << diag_.readsZeroGenes << "\n";
    out << "readsNAlignWithGeneZero\t" << diag_.readsNAlignWithGeneZero << "\n";
    out << "readsSumWeightLessThanOne\t" << diag_.readsSumWeightLessThanOne << "\n";
    out << "\nnTrDistribution:\n";
    for (const auto& kv : diag_.nTrDistribution) {
        out << "nTr_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\ngeneSetSizeDistribution:\n";
    for (const auto& kv : diag_.geneSetSizeDistribution) {
        out << "geneSetSize_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\nnAlignWithGeneDistribution:\n";
    for (const auto& kv : diag_.nAlignWithGeneDistribution) {
        out << "nAlignWithGene_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\nsumWeightDistribution (bucketed):\n";
    for (const auto& kv : diag_.sumWeightDistribution) {
        double bucketStart = kv.first * 0.1;
        double bucketEnd = (kv.first == 10) ? 999.0 : (kv.first + 1) * 0.1;
        out << "sumWeight_" << bucketStart << "_to_" << bucketEnd << "\t" << kv.second << "\n";
    }
    
    // Compat mode counters
    if (diag_.compatAlignsReclassifiedIntronic > 0 || diag_.compatAlignsLenientAccepted > 0 ||
        diag_.compatAlignsOverlapWeightApplied > 0 || diag_.compatPositionsSkippedOverlap > 0 ||
        diag_.compatPositionsSkippedTrim > 0) {
        out << "\nCompatModeCounters:\n";
        out << "compatAlignsReclassifiedIntronic\t" << diag_.compatAlignsReclassifiedIntronic << "\n";
        out << "compatAlignsLenientAccepted\t" << diag_.compatAlignsLenientAccepted << "\n";
        out << "compatAlignsOverlapWeightApplied\t" << diag_.compatAlignsOverlapWeightApplied << "\n";
        out << "compatPositionsSkippedOverlap\t" << diag_.compatPositionsSkippedOverlap << "\n";
        out << "compatPositionsSkippedTrim\t" << diag_.compatPositionsSkippedTrim << "\n";
    }
}

void SlamQuant::writeDebug(const Transcriptome& tr, double errorRate, double convRate) const {
    if (!debugEnabled_) {
        return;
    }
    std::string base = debugOutPrefix_.empty() ? "SlamQuant.debug" : debugOutPrefix_;
    if (!debugGeneStats_.empty()) {
        std::ofstream out((base + ".gene.tsv").c_str());
        if (out.good()) {
            out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tk_over_nT\t"
                << "DebugExonicWeight\tDebugIntronicWeight\tDebugSenseWeight\tDebugAntisenseWeight\t"
                << "DebugReadsAssigned\tDebugReadsIntronic\tDropsSnpMask\tDropsStrandness\n";
            SlamSolver solver(errorRate, convRate);
            for (size_t i = 0; i < debugGeneStats_.size(); ++i) {
                if (debugGeneMask_.empty() || debugGeneMask_[i] == 0) {
                    continue;
                }
                const SlamGeneStats& stats = geneStats_[i];
                SlamResult res = solver.solve(stats.histogram);
                const std::string& geneId = tr.geID[i];
                const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
                const SlamDebugGeneStats& dbg = debugGeneStats_[i];
                double kOverNT = (stats.coverage > 0.0) ? (stats.conversions / stats.coverage) : 0.0;
                out << geneId << "\t"
                    << geneName << "\t"
                    << stats.readCount << "\t"
                    << stats.conversions << "\t"
                    << stats.coverage << "\t"
                    << res.ntr << "\t"
                    << kOverNT << "\t"
                    << dbg.exonicWeight << "\t"
                    << dbg.intronicWeight << "\t"
                    << dbg.senseWeight << "\t"
                    << dbg.antisenseWeight << "\t"
                    << dbg.readsAssigned << "\t"
                    << dbg.readsAssignedIntronic << "\t"
                    << dbg.dropsSnpMask << "\t"
                    << dbg.dropsStrandness << "\n";
            }
        }
    }
    if (!debugReadRecords_.empty() && debugMaxReads_ != 0) {
        std::ofstream out((base + ".reads.tsv").c_str());
        if (out.good()) {
            out << "ReadName\tReadLoc\tGene\tSymbol\tStatus\tIntronic\tOppositeStrand\tWeight\t"
                << "nT\tk\tReadLength\tSnpBuffered\tConvReadPos\tConvGenPos\n";
            for (const auto& rec : debugReadRecords_) {
                std::string geneId = "NA";
                std::string geneName = "NA";
                if (rec.geneId < tr.geID.size()) {
                    geneId = tr.geID[rec.geneId];
                    const std::string& gname = tr.geName[rec.geneId];
                    geneName = gname.empty() ? geneId : gname;
                }
                out << rec.readName << "\t"
                    << rec.readLoc << "\t"
                    << geneId << "\t"
                    << geneName << "\t"
                    << debugDropReasonName(rec.status) << "\t"
                    << (rec.intronic ? 1 : 0) << "\t"
                    << (rec.oppositeStrand ? 1 : 0) << "\t"
                    << rec.weight << "\t"
                    << rec.nT << "\t"
                    << rec.k << "\t"
                    << rec.readLength << "\t"
                    << (rec.snpBuffered ? 1 : 0) << "\t"
                    << rec.convReadPos << "\t"
                    << rec.convGenPos << "\n";
            }
        }
    }
}

void SlamQuant::writeTransitions(const std::string& outFile) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Genomic\tRead\tCoverage\tMismatches\n";
    const SlamTransitionStats& stats = transitions_[static_cast<size_t>(SlamMismatchCategory::Exonic)];
    for (int g = 0; g < 4; ++g) {
        for (int r = 0; r < 4; ++r) {
            if (g == r) {
                continue;
            }
            out << bases[g] << "\t" << bases[r] << "\t"
                << stats.coverage[g] << "\t"
                << stats.mismatches[g][r] << "\n";
        }
    }
}

void SlamQuant::writeMismatches(const std::string& outFile, const std::string& condition) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        // Try to create parent directory if it doesn't exist
        size_t lastSlash = outFile.find_last_of("/\\");
        if (lastSlash != std::string::npos) {
            std::string dir = outFile.substr(0, lastSlash);
            system(("mkdir -p " + dir).c_str());
            out.open(outFile.c_str());
        }
        if (!out.good()) {
            return;
        }
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Category\tCondition\tOrientation\tGenomic\tRead\tCoverage\tMismatches\n";
    const char* labels[2] = {"First", "Second"};
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        const char* category = slamMismatchCategoryName(static_cast<SlamMismatchCategory>(c));
        const SlamTransitionStats* orientations[2] = {&transitionsFirst_[c], &transitionsSecond_[c]};
        for (int o = 0; o < 2; ++o) {
            const SlamTransitionStats& stats = *orientations[o];
            for (int g = 0; g < 4; ++g) {
                for (int r = 0; r < 4; ++r) {
                    if (g == r) {
                        continue;
                    }
                    out << category << "\t" << condition << "\t" << labels[o] << "\t"
                        << bases[g] << "\t" << bases[r] << "\t"
                        << stats.coverage[g] << "\t" << stats.mismatches[g][r] << "\n";
                }
            }
        }
    }
}

void SlamQuant::writeMismatchDetails(const std::string& outFile) const {
    std::ofstream out(outFile.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        return;
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Category\tGenomic\tRead\tPosition\tOverlap\tOpposite\tCoverage\tMismatches\n";
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        const char* category = slamMismatchCategoryName(static_cast<SlamMismatchCategory>(c));
        const auto& positions = positionTransitions_[c];
        for (size_t pos = 0; pos < positions.size(); ++pos) {
            const SlamPositionStats& stats = positions[pos];
            for (int ov = 0; ov < 2; ++ov) {
                for (int opp = 0; opp < 2; ++opp) {
                    for (int g = 0; g < 4; ++g) {
                        double cov = stats.coverage[ov][opp][g];
                        if (cov <= 0.0) {
                            continue;
                        }
                        for (int r = 0; r < 4; ++r) {
                            if (g == r) {
                                continue;
                            }
                            out << category << "\t" << bases[g] << "\t" << bases[r] << "\t"
                                << pos << "\t" << ov << "\t" << opp << "\t" << cov << "\t"
                                << stats.mismatches[ov][opp][g][r] << "\n";
                        }
                    }
                }
            }
        }
    }
    out.close();
}

void SlamQuant::writeTopMismatches(const Transcriptome& tr, const std::string& refFile,
                                   const std::string& mismatchFile, size_t topN) const {
    // Parse reference file (handle gzipped files)
    std::unordered_map<std::string, std::map<std::string, double>> refData;
    std::ifstream ref(refFile.c_str(), std::ios::binary);
    if (!ref.good()) {
        return;  // Skip if reference file not found
    }
    
    // Check if gzipped by reading magic bytes
    unsigned char magic[2];
    ref.read(reinterpret_cast<char*>(magic), 2);
    ref.seekg(0);
    bool isGzip = (magic[0] == 0x1f && magic[1] == 0x8b);
    
    if (isGzip) {
        // For gzipped files, use zcat via popen
        std::string cmd = "zcat " + refFile;
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            return;
        }
        char buffer[4096];
        std::string line;
        bool headerSkipped = false;
        while (fgets(buffer, sizeof(buffer), pipe)) {
            line = buffer;
            if (line.back() == '\n') line.pop_back();
            if (!headerSkipped) {
                headerSkipped = true;
                continue;
            }
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string geneId, symbol, readCount, conversions, coverage, ntr;
            if (iss >> geneId >> symbol >> readCount >> conversions >> coverage >> ntr) {
                try {
                    refData[geneId]["ReadCount"] = std::stod(readCount);
                    refData[geneId]["Conversions"] = std::stod(conversions);
                    refData[geneId]["Coverage"] = std::stod(coverage);
                    refData[geneId]["NTR"] = std::stod(ntr);
                } catch (...) {
                    continue;  // Skip malformed lines
                }
            }
        }
        pclose(pipe);
    } else {
        // Plain text file
        std::string line;
        std::getline(ref, line);  // Skip header
        while (std::getline(ref, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string geneId, symbol, readCount, conversions, coverage, ntr;
            if (iss >> geneId >> symbol >> readCount >> conversions >> coverage >> ntr) {
                try {
                    refData[geneId]["ReadCount"] = std::stod(readCount);
                    refData[geneId]["Conversions"] = std::stod(conversions);
                    refData[geneId]["Coverage"] = std::stod(coverage);
                    refData[geneId]["NTR"] = std::stod(ntr);
                } catch (...) {
                    continue;  // Skip malformed lines
                }
            }
        }
    }

    // Collect mismatches
    std::vector<std::pair<double, size_t>> mismatches;
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        const SlamGeneStats& stats = geneStats_[i];
        const std::string& geneId = tr.geID[i];
        auto it = refData.find(geneId);
        if (it != refData.end()) {
            double delta = std::abs(stats.readCount - it->second.at("ReadCount"));
            if (delta > 0.1) {
                mismatches.push_back({delta, i});
            }
        }
    }
    std::sort(mismatches.rbegin(), mismatches.rend());
    if (mismatches.size() > topN) {
        mismatches.resize(topN);
    }

    // Write mismatch file
    std::ofstream out(mismatchFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Gene\tSymbol\tReadCount_ref\tReadCount_test\tReadCount_delta\t"
        << "Conversions_ref\tConversions_test\tConversions_delta\t"
        << "Coverage_ref\tCoverage_test\tCoverage_delta\t"
        << "NTR_ref\tNTR_test\tNTR_delta\n";
    
    SlamSolver solver(0.001, 0.05);  // Default rates
    for (const auto& mismatch : mismatches) {
        size_t i = mismatch.second;
        const SlamGeneStats& stats = geneStats_[i];
        const std::string& geneId = tr.geID[i];
        const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
        auto it = refData.find(geneId);
        if (it == refData.end()) continue;
        
        SlamResult res = solver.solve(stats.histogram);
        const auto& ref = it->second;
        
        out << geneId << "\t" << geneName << "\t"
            << ref.at("ReadCount") << "\t" << stats.readCount << "\t"
            << (stats.readCount - ref.at("ReadCount")) << "\t"
            << ref.at("Conversions") << "\t" << stats.conversions << "\t"
            << (stats.conversions - ref.at("Conversions")) << "\t"
            << ref.at("Coverage") << "\t" << stats.coverage << "\t"
            << (stats.coverage - ref.at("Coverage")) << "\t"
            << ref.at("NTR") << "\t" << res.ntr << "\t"
            << (res.ntr - ref.at("NTR")) << "\n";
    }
}

bool SlamSnpMask::loadBed(const std::string& path, const Genome& genome, std::string* err) {
    positions_.clear();
    std::ifstream in(path.c_str());
    if (!in.good()) {
        if (err) *err = "Failed to open SNP BED: " + path;
        return false;
    }

    std::unordered_map<std::string, uint32_t> chrMap;
    chrMap.reserve(genome.chrName.size());
    for (uint32_t i = 0; i < genome.chrName.size(); ++i) {
        chrMap[genome.chrName[i]] = i;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string chr;
        uint64_t start = 0;
        uint64_t end = 0;
        if (!(iss >> chr >> start >> end)) {
            continue;
        }
        auto it = chrMap.find(chr);
        if (it == chrMap.end()) {
            continue;
        }
        uint64_t chrStart = genome.chrStart[it->second];
        if (end <= start) {
            continue;
        }
        uint64_t maxLen = end - start;
        if (maxLen > 1000) {
            // Avoid pathological regions; SNP BEDs should be single-base.
            continue;
        }
        for (uint64_t pos = start; pos < end; ++pos) {
            positions_.insert(chrStart + pos);
        }
    }
    return true;
}
