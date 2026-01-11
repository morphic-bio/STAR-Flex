#ifndef SLAM_VARIANCE_ANALYSIS_H
#define SLAM_VARIANCE_ANALYSIS_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

// Per-position statistics for variance analysis
struct SlamPositionVarianceStats {
    uint64_t readCount = 0;           // Number of reads covering this position
    uint64_t tCount = 0;              // Number of T bases at this position
    uint64_t tcCount = 0;             // Number of T→C conversions at this position
    double qualSum = 0.0;              // Sum of quality scores
    uint64_t qualCount = 0;            // Count of quality scores
    double qualSumSq = 0.0;            // Sum of squares for variance calculation
    double tcRateSum = 0.0;           // Sum of T→C rates (for per-read averaging)
    uint64_t tcRateCount = 0;         // Count of reads with T bases
    
    // Computed statistics
    double meanQual() const {
        return qualCount > 0 ? qualSum / qualCount : 0.0;
    }
    
    double varianceQual() const {
        if (qualCount < 2) return 0.0;
        double mean = meanQual();
        return (qualSumSq / qualCount) - (mean * mean);
    }
    
    double stddevQual() const {
        return std::sqrt(varianceQual());
    }
    
    double meanTcRate() const {
        return tcRateCount > 0 ? tcRateSum / tcRateCount : 0.0;
    }
};

// Variance analysis results
struct SlamVarianceTrimResult {
    int trim5p = 0;
    int trim3p = 0;
    bool success = false;
    std::string mode;                  // "auto", "auto_fallback", "manual"
    uint64_t readsAnalyzed = 0;
    uint32_t kneeBin5p = 0;
    uint32_t kneeBin3p = 0;
};

// Variance analyzer for auto-trim
class SlamVarianceAnalyzer {
public:
    SlamVarianceAnalyzer(uint32_t maxReads = 100000, uint32_t minReads = 1000);
    
    // Record a read (call once per read before recording positions)
    bool recordRead();
    
    // Record a position observation (only if readsAnalyzed_ < maxReads_)
    void recordPosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc);
    
    // Compute trim values from variance curves
    SlamVarianceTrimResult computeTrim(uint32_t readLength);
    
    // Get per-position statistics
    const std::unordered_map<uint32_t, SlamPositionVarianceStats>& getStats() const {
        return positionStats_;
    }
    
    // Get number of reads analyzed
    uint64_t readsAnalyzed() const { return readsAnalyzed_; }
    
    // Reset for new file
    void reset();
    
    // Merge stats from another analyzer (for thread aggregation)
    void merge(const SlamVarianceAnalyzer& other);
    
private:
    std::unordered_map<uint32_t, SlamPositionVarianceStats> positionStats_;
    uint64_t readsAnalyzed_;
    uint32_t maxReads_;
    uint32_t minReads_;
    
    // Knee detection helper (similar to SNP threshold estimation)
    static std::pair<uint32_t, uint32_t> detectKnee(
        const std::vector<double>& varianceCurve,
        uint32_t maxPos);
};

#endif // SLAM_VARIANCE_ANALYSIS_H
