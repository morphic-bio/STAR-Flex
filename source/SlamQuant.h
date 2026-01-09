#ifndef SLAM_QUANT_H
#define SLAM_QUANT_H

#include "SlamSolver.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <array>
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

struct SlamSnpBufferStats {
    uint64_t bufferedReads = 0;
    uint64_t bufferedMismatches = 0;
    uint64_t bufferedMismatchesKept = 0;
    uint64_t bufferEntries = 0;
    uint64_t bufferBytes = 0;
    uint64_t maskEntries = 0;
    uint64_t blacklistEntries = 0;
    double avgMismatches = 0.0;
    double avgMismatchesKept = 0.0;
};

struct SlamTransitionStats {
    double coverage[4] = {0.0, 0.0, 0.0, 0.0};
    double mismatches[4][4] = {{0.0}};
};

struct SlamPositionStats {
    double coverage[2][2][4] = {};
    double mismatches[2][2][4][4] = {};
};

enum class SlamMismatchCategory : uint8_t {
    Exonic = 0,
    ExonicSense = 1,
    Intronic = 2,
    IntronicSense = 3,
    Count
};

constexpr size_t kSlamMismatchCategoryCount = static_cast<size_t>(SlamMismatchCategory::Count);
const char* slamMismatchCategoryName(SlamMismatchCategory cat);

class SlamQuant {
public:
    explicit SlamQuant(uint32_t nGenes, bool snpDetect = false);
    SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes, bool snpDetect = false);

    void addRead(uint32_t geneId, uint16_t nT, uint8_t k, double weight);
    void addTransitionBase(SlamMismatchCategory category, uint32_t readPos, bool secondMate,
                           bool overlap, bool opposite, int genomicBase, int readBase, double weight);
    bool snpDetectEnabled() const { return snpDetectEnabled_; }
    void recordSnpObservation(uint64_t pos, bool isMismatch);
    void bufferSnpRead(uint32_t geneId, uint16_t nT,
                       const std::vector<uint32_t>& mismatchPositions, double weight);
    void finalizeSnpMask(SlamSnpBufferStats* outStats = nullptr);
    void merge(const SlamQuant& other);
    void write(const Transcriptome& tr, const std::string& outFile,
               double errorRate, double convRate) const;
    void writeDiagnostics(const std::string& diagFile) const;
    void writeTransitions(const std::string& outFile) const;
    void writeMismatches(const std::string& outFile, const std::string& condition) const;
    void writeMismatchDetails(const std::string& outFile) const;
    void writeTopMismatches(const Transcriptome& tr, const std::string& refFile,
                           const std::string& mismatchFile, size_t topN) const;

    const std::vector<SlamGeneStats>& genes() const { return geneStats_; }
    SlamDiagnostics& diagnostics() { return diag_; }
    const SlamDiagnostics& diagnostics() const { return diag_; }

private:
    std::vector<SlamGeneStats> geneStats_;
    SlamDiagnostics diag_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitions_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitionsFirst_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitionsSecond_;
    std::array<std::vector<SlamPositionStats>, kSlamMismatchCategoryCount> positionTransitions_;
    bool snpDetectEnabled_ = false;
    bool snpFinalized_ = false;
    std::unordered_map<uint64_t, uint32_t> snpMask_;
    std::vector<uint32_t> snpReadBuffer_;
    std::vector<double> snpReadWeights_;
    std::vector<uint8_t> allowedGenes_;
};

#endif
