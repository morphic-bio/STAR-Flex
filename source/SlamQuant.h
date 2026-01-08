#ifndef SLAM_QUANT_H
#define SLAM_QUANT_H

#include "SlamSolver.h"

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

class Genome;
class Transcriptome;

class SlamSnpMask {
public:
    bool loadBed(const std::string& path, const Genome& genome, std::string* err);
    bool contains(uint64_t pos) const { return positions_.count(pos) > 0; }
    size_t size() const { return positions_.size(); }

private:
    std::unordered_set<uint64_t> positions_;
};

struct SlamGeneStats {
    MismatchHistogram histogram;
    double readCount = 0.0;
    double conversions = 0.0;
    double coverage = 0.0;
};

struct SlamDiagnostics {
    uint64_t readsDroppedSnpMask = 0;
    uint64_t readsZeroGenes = 0;
    uint64_t readsProcessed = 0;
    uint64_t readsNAlignWithGeneZero = 0;  // reads where nAlignWithGene == 0
    uint64_t readsSumWeightLessThanOne = 0;  // reads where sumWeight < 1.0
    std::map<size_t, uint64_t> nTrDistribution;  // nTr -> count
    std::map<size_t, uint64_t> geneSetSizeDistribution;  // gene set size -> count
    std::map<size_t, uint64_t> nAlignWithGeneDistribution;  // nAlignWithGene -> count
    std::map<double, uint64_t> sumWeightDistribution;  // sumWeight bucket -> count (bucketed)
};

struct SlamTransitionStats {
    double coverage[4] = {0.0, 0.0, 0.0, 0.0};
    double mismatches[4][4] = {{0.0}};
};

class SlamQuant {
public:
    explicit SlamQuant(uint32_t nGenes);
    SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes);

    void addRead(uint32_t geneId, uint16_t nT, uint8_t k, double weight);
    void addTransitions(const double coverage[4], const double mismatches[4][4], double weight);
    void merge(const SlamQuant& other);
    void write(const Transcriptome& tr, const std::string& outFile,
               double errorRate, double convRate) const;
    void writeDiagnostics(const std::string& diagFile) const;
    void writeTransitions(const std::string& outFile) const;
    void writeTopMismatches(const Transcriptome& tr, const std::string& refFile,
                           const std::string& mismatchFile, size_t topN) const;

    const std::vector<SlamGeneStats>& genes() const { return geneStats_; }
    SlamDiagnostics& diagnostics() { return diag_; }
    const SlamDiagnostics& diagnostics() const { return diag_; }

private:
    std::vector<SlamGeneStats> geneStats_;
    SlamDiagnostics diag_;
    SlamTransitionStats transitions_;
    std::vector<uint8_t> allowedGenes_;
};

#endif
