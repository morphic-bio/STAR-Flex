#include "SlamVarianceAnalysis.h"
#include <algorithm>
#include <cmath>

namespace {
constexpr uint32_t kTrimHistogramBins = 100;
constexpr double kTrimKneeEpsilon = 0.02;
constexpr uint32_t kTrimMinClamp = 0;
constexpr uint32_t kTrimMaxClamp = 50;  // Max trim of 50 bases
constexpr double kTrimMinVariance = 0.1;  // Minimum variance to consider
}

SlamVarianceAnalyzer::SlamVarianceAnalyzer(uint32_t maxReads, uint32_t minReads)
    : readsAnalyzed_(0), maxReads_(maxReads), minReads_(minReads) {
}

bool SlamVarianceAnalyzer::recordRead() {
    // Check if we've reached max reads
    if (maxReads_ > 0 && readsAnalyzed_ >= maxReads_) {
        return false; // Stop collecting
    }
    readsAnalyzed_++;
    return true; // Continue collecting
}

void SlamVarianceAnalyzer::recordPosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc) {
    // Gate on readsAnalyzed_ <= maxReads_ (allow Nth read's positions)
    if (maxReads_ > 0 && readsAnalyzed_ > maxReads_) {
        return; // Stop collecting positions
    }
    
    auto& stats = positionStats_[readPos];
    stats.readCount++;
    
    // Record quality statistics
    double q = static_cast<double>(qual);
    stats.qualSum += q;
    stats.qualCount++;
    stats.qualSumSq += q * q;
    
    // Record T→C statistics
    if (isT) {
        stats.tCount++;
        if (isTc) {
            stats.tcCount++;
        }
        // Track per-T-base T→C indicator (0 or 1) for variance calculation
        double tcRate = isTc ? 1.0 : 0.0;
        stats.tcRateSum += tcRate;
        stats.tcRateSumSq += tcRate * tcRate;  // For variance: sum of squares
        stats.tcRateCount++;
    }
}

void SlamVarianceAnalyzer::reset() {
    positionStats_.clear();
    readsAnalyzed_ = 0;
}

std::pair<uint32_t, uint32_t> SlamVarianceAnalyzer::detectKnee(
    const std::vector<double>& varianceCurve,
    uint32_t maxPos) {
    
    if (varianceCurve.empty() || maxPos == 0) {
        return {0, 0};
    }
    
    // Build cumulative distribution from high variance to low
    // Similar to SNP threshold estimation, but for variance
    std::vector<double> cumulative(varianceCurve.size(), 0.0);
    double maxVar = 0.0;
    for (size_t i = 0; i < varianceCurve.size(); ++i) {
        if (varianceCurve[i] > maxVar) {
            maxVar = varianceCurve[i];
        }
    }
    
    if (maxVar < kTrimMinVariance) {
        return {0, 0};  // Too low variance, fallback
    }
    
    // Normalize and compute cumulative
    for (size_t i = 0; i < varianceCurve.size(); ++i) {
        double normVar = varianceCurve[i] / maxVar;
        cumulative[i] = normVar;
    }
    
    // Find knee using distance from diagonal (Kneedle algorithm)
    // For 5' trim: look for high variance at start (positions 0-N)
    // For 3' trim: look for high variance at end (positions L-N to L)
    
    double maxDist5p = -1.0;
    uint32_t knee5p = 0;
    
    // Scan from start for 5' trim
    for (uint32_t i = 0; i < std::min(static_cast<uint32_t>(varianceCurve.size()), maxPos / 2); ++i) {
        double xNorm = static_cast<double>(i) / (maxPos / 2);
        double yNorm = cumulative[i];
        double dist = yNorm + xNorm - 1.0;  // Distance from diagonal
        if (dist > maxDist5p) {
            maxDist5p = dist;
            knee5p = i;
        }
    }
    
    double maxDist3p = -1.0;
    uint32_t knee3p = 0;
    
    // Scan from end for 3' trim
    uint32_t start3p = varianceCurve.size() > maxPos / 2 
        ? static_cast<uint32_t>(varianceCurve.size() - maxPos / 2) 
        : 0;
    for (uint32_t i = start3p; i < varianceCurve.size(); ++i) {
        uint32_t posFromEnd = static_cast<uint32_t>(varianceCurve.size() - 1 - i);
        double xNorm = static_cast<double>(posFromEnd) / (maxPos / 2);
        double yNorm = cumulative[i];
        double dist = yNorm + xNorm - 1.0;
        if (dist > maxDist3p) {
            maxDist3p = dist;
            knee3p = i;
        }
    }
    
    // Check knee strength
    if (maxDist5p < kTrimKneeEpsilon) {
        knee5p = 0;
    }
    if (maxDist3p < kTrimKneeEpsilon) {
        knee3p = 0;
    }
    
    // Clamp to valid range
    knee5p = std::min(knee5p, kTrimMaxClamp);
    knee3p = std::min(static_cast<uint32_t>(varianceCurve.size() - knee3p), kTrimMaxClamp);
    
    return {knee5p, knee3p};
}

void SlamVarianceAnalyzer::merge(const SlamVarianceAnalyzer& other) {
    // Merge read count (sum across threads, but cap at maxReads_)
    readsAnalyzed_ = std::min(readsAnalyzed_ + other.readsAnalyzed_, 
                              static_cast<uint64_t>(maxReads_ > 0 ? maxReads_ : UINT64_MAX));
    
    // Merge position stats
    for (const auto& kv : other.positionStats_) {
        auto& stats = positionStats_[kv.first];
        const auto& otherStats = kv.second;
        
        stats.readCount += otherStats.readCount;
        stats.tCount += otherStats.tCount;
        stats.tcCount += otherStats.tcCount;
        stats.qualSum += otherStats.qualSum;
        stats.qualCount += otherStats.qualCount;
        stats.qualSumSq += otherStats.qualSumSq;
        stats.tcRateSum += otherStats.tcRateSum;
        stats.tcRateSumSq += otherStats.tcRateSumSq;
        stats.tcRateCount += otherStats.tcRateCount;
    }
}

SlamVarianceTrimResult SlamVarianceAnalyzer::computeTrim(uint32_t readLength) {
    SlamVarianceTrimResult result;
    result.readsAnalyzed = readsAnalyzed_;
    
    if (readsAnalyzed_ < minReads_) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    if (positionStats_.empty()) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Build variance curve from position statistics
    // We use the VARIANCE (or stdev) of T→C rate across reads at each position
    // High variance indicates inconsistent conversion - likely artifacts
    // Low variance indicates consistent signal - either true conversions or no conversions
    std::vector<double> varianceCurve(readLength, 0.0);
    uint32_t maxPos = 0;
    
    for (const auto& kv : positionStats_) {
        uint32_t pos = kv.first;
        if (pos >= readLength) {
            continue;
        }
        const auto& stats = kv.second;
        
        // Primary signal: standard deviation of T→C rate at this position
        // High stdev = inconsistent conversion = artifact
        // Low stdev = consistent (either all convert or none convert) = signal
        double tcStddev = stats.stddevTcRate();
        
        // Secondary signal: quality score variance (optional, lower weight)
        // High quality variance can also indicate problematic positions
        double qualStddev = stats.stddevQual();
        
        // Combined signal: primarily T→C stdev, with small quality contribution
        // T→C stdev is on 0-0.5 scale (max variance for Bernoulli is 0.25, stdev = 0.5)
        // Quality stdev is on 0-~10 scale typically
        double combinedVariance = tcStddev + (qualStddev / 100.0);  // Normalize quality contribution
        
        varianceCurve[pos] = combinedVariance;
        if (pos > maxPos) {
            maxPos = pos;
        }
    }
    
    if (maxPos == 0) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Detect knees
    auto [trim5p, trim3pPos] = detectKnee(varianceCurve, maxPos + 1);
    
    // Convert 3' position to trim amount
    uint32_t trim3p = (trim3pPos < readLength) ? (readLength - trim3pPos) : 0;
    
    result.trim5p = static_cast<int>(trim5p);
    result.trim3p = static_cast<int>(trim3p);
    result.kneeBin5p = trim5p;
    result.kneeBin3p = trim3pPos;
    result.mode = "auto";
    result.success = true;
    
    return result;
}
