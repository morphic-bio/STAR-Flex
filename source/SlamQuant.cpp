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

SlamQuant::SlamQuant(uint32_t nGenes) : geneStats_(nGenes) {}

SlamQuant::SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes)
    : geneStats_(nGenes), allowedGenes_(std::move(allowedGenes)) {}

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

void SlamQuant::addTransitions(const double coverage[4], const double mismatches[4][4], double weight) {
    if (weight <= 0.0) {
        return;
    }
    for (int g = 0; g < 4; ++g) {
        transitions_.coverage[g] += coverage[g] * weight;
        for (int r = 0; r < 4; ++r) {
            if (g == r) {
                continue;
            }
            transitions_.mismatches[g][r] += mismatches[g][r] * weight;
        }
    }
}

void SlamQuant::merge(const SlamQuant& other) {
    if (other.geneStats_.size() != geneStats_.size()) {
        return;
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
    diag_.readsZeroGenes += other.diag_.readsZeroGenes;
    diag_.readsProcessed += other.diag_.readsProcessed;
    diag_.readsNAlignWithGeneZero += other.diag_.readsNAlignWithGeneZero;
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
    for (int g = 0; g < 4; ++g) {
        transitions_.coverage[g] += other.transitions_.coverage[g];
        for (int r = 0; r < 4; ++r) {
            transitions_.mismatches[g][r] += other.transitions_.mismatches[g][r];
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
}

void SlamQuant::writeTransitions(const std::string& outFile) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Genomic\tRead\tCoverage\tMismatches\n";
    for (int g = 0; g < 4; ++g) {
        for (int r = 0; r < 4; ++r) {
            if (g == r) {
                continue;
            }
            out << bases[g] << "\t" << bases[r] << "\t"
                << transitions_.coverage[g] << "\t"
                << transitions_.mismatches[g][r] << "\n";
        }
    }
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
