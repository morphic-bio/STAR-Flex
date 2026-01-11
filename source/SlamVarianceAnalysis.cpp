#include "SlamVarianceAnalysis.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

SlamVarianceAnalyzer::SlamVarianceAnalyzer(uint32_t maxReads, uint32_t minReads,
                                           uint32_t smoothWindow, uint32_t minSegLen, uint32_t maxTrim)
    : readsAnalyzed_(0), maxReads_(maxReads), minReads_(minReads),
      smoothWindow_(smoothWindow), minSegLen_(minSegLen), maxTrim_(maxTrim) {
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

std::vector<double> SlamVarianceAnalyzer::smoothMedian(const std::vector<double>& values, uint32_t window) {
    if (window <= 1 || values.empty()) {
        return values;
    }
    
    std::vector<double> result(values.size());
    int half = static_cast<int>(window / 2);
    
    for (size_t i = 0; i < values.size(); ++i) {
        int start = std::max(0, static_cast<int>(i) - half);
        int end = std::min(static_cast<int>(values.size()), static_cast<int>(i) + half + 1);
        
        // Collect non-NaN values in window
        std::vector<double> chunk;
        for (int j = start; j < end; ++j) {
            if (!std::isnan(values[j])) {
                chunk.push_back(values[j]);
            }
        }
        
        if (chunk.empty()) {
            result[i] = std::nan("");
        } else {
            // Compute median
            std::sort(chunk.begin(), chunk.end());
            size_t mid = chunk.size() / 2;
            if (chunk.size() % 2 == 0) {
                result[i] = (chunk[mid - 1] + chunk[mid]) / 2.0;
            } else {
                result[i] = chunk[mid];
            }
        }
    }
    
    return result;
}

void SlamVarianceAnalyzer::interpolateMissing(std::vector<double>& values) {
    if (values.empty()) return;
    
    // Find valid indices
    std::vector<size_t> validIdx;
    for (size_t i = 0; i < values.size(); ++i) {
        if (!std::isnan(values[i])) {
            validIdx.push_back(i);
        }
    }
    
    if (validIdx.size() < 2) {
        // Not enough valid values for interpolation, fill with 0
        for (auto& v : values) {
            if (std::isnan(v)) v = 0.0;
        }
        return;
    }
    
    // Linear interpolation for NaN values
    for (size_t i = 0; i < values.size(); ++i) {
        if (std::isnan(values[i])) {
            // Find surrounding valid indices
            size_t left = 0, right = validIdx.size() - 1;
            for (size_t j = 0; j < validIdx.size(); ++j) {
                if (validIdx[j] < i) left = j;
                if (validIdx[j] > i) { right = j; break; }
            }
            
            size_t li = validIdx[left];
            size_t ri = validIdx[right];
            
            if (li == ri) {
                values[i] = values[li];
            } else {
                // Linear interpolation
                double t = static_cast<double>(i - li) / static_cast<double>(ri - li);
                values[i] = values[li] + t * (values[ri] - values[li]);
            }
        }
    }
}

SegmentFit SlamVarianceAnalyzer::fitSegment(
    const std::vector<double>& prefixX,
    const std::vector<double>& prefixXX,
    const std::vector<double>& prefixY,
    const std::vector<double>& prefixXY,
    const std::vector<double>& prefixYY,
    const std::vector<double>& y,
    uint32_t start, uint32_t end) {
    
    SegmentFit fit;
    
    if (end < start) {
        return fit;
    }
    
    uint32_t n = end - start + 1;
    if (n < 2) {
        fit.intercept = y[start];
        fit.slope = 0.0;
        fit.sse = 0.0;
        return fit;
    }
    
    // Use prefix sums: sum from [start, end] = prefix[end+1] - prefix[start]
    double Sx = prefixX[end + 1] - prefixX[start];
    double Sxx = prefixXX[end + 1] - prefixXX[start];
    double Sy = prefixY[end + 1] - prefixY[start];
    double Sxy = prefixXY[end + 1] - prefixXY[start];
    double Syy = prefixYY[end + 1] - prefixYY[start];
    
    double den = n * Sxx - Sx * Sx;
    if (std::abs(den) < 1e-9) {
        fit.slope = 0.0;
    } else {
        fit.slope = (n * Sxy - Sx * Sy) / den;
    }
    fit.intercept = (Sy - fit.slope * Sx) / n;
    
    // SSE = sum(y^2) + m^2*sum(x^2) + n*b^2 + 2*m*b*sum(x) - 2*m*sum(xy) - 2*b*sum(y)
    fit.sse = Syy + fit.slope * fit.slope * Sxx + n * fit.intercept * fit.intercept
              + 2.0 * fit.slope * fit.intercept * Sx
              - 2.0 * fit.slope * Sxy
              - 2.0 * fit.intercept * Sy;
    
    // Clamp negative SSE (numerical precision)
    if (fit.sse < 0.0) fit.sse = 0.0;
    
    return fit;
}

std::tuple<uint32_t, uint32_t, double, SegmentFit, SegmentFit, SegmentFit>
SlamVarianceAnalyzer::segmentedRegression(const std::vector<double>& y, uint32_t minSegLen, uint32_t maxTrim) {
    uint32_t n = static_cast<uint32_t>(y.size());
    
    // Default return: no breakpoints found
    if (n < minSegLen * 3 + 1) {
        return {0, n - 1, std::numeric_limits<double>::max(), {}, {}, {}};
    }
    
    // Build prefix sums for efficient segment fitting
    std::vector<double> prefixX(n + 1, 0.0);
    std::vector<double> prefixXX(n + 1, 0.0);
    std::vector<double> prefixY(n + 1, 0.0);
    std::vector<double> prefixXY(n + 1, 0.0);
    std::vector<double> prefixYY(n + 1, 0.0);
    
    for (uint32_t i = 0; i < n; ++i) {
        double x = static_cast<double>(i);
        prefixX[i + 1] = prefixX[i] + x;
        prefixXX[i + 1] = prefixXX[i] + x * x;
        prefixY[i + 1] = prefixY[i] + y[i];
        prefixXY[i + 1] = prefixXY[i] + x * y[i];
        prefixYY[i + 1] = prefixYY[i] + y[i] * y[i];
    }
    
    // Candidate breakpoint ranges
    // b1: first breakpoint (end of segment 1, segment 2 starts at b1)
    // b2: second breakpoint (end of segment 2, segment 3 starts at b2+1)
    // Segments: [0, b1-1], [b1, b2], [b2+1, n-1]
    
    uint32_t b1Min = minSegLen;
    uint32_t b1Max = std::min(maxTrim, n - 2 * minSegLen - 1);
    uint32_t b2MinFloor = std::max(minSegLen - 1, n - 1 - maxTrim);
    uint32_t b2Max = n - minSegLen - 1;
    
    if (b1Max < b1Min || b2MinFloor > b2Max) {
        return {0, n - 1, std::numeric_limits<double>::max(), {}, {}, {}};
    }
    
    double bestSSE = std::numeric_limits<double>::max();
    uint32_t bestB1 = 0, bestB2 = n - 1;
    SegmentFit bestSeg1, bestSeg2, bestSeg3;
    
    for (uint32_t b1 = b1Min; b1 <= b1Max; ++b1) {
        uint32_t b2Min = std::max(b1 + minSegLen - 1, b2MinFloor);
        
        for (uint32_t b2 = b2Min; b2 <= b2Max; ++b2) {
            // Segment 1: [0, b1-1]
            SegmentFit seg1 = fitSegment(prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, 0, b1 - 1);
            // Segment 2: [b1, b2]
            SegmentFit seg2 = fitSegment(prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b1, b2);
            // Segment 3: [b2+1, n-1]
            SegmentFit seg3 = fitSegment(prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b2 + 1, n - 1);
            
            double totalSSE = seg1.sse + seg2.sse + seg3.sse;
            
            if (totalSSE < bestSSE) {
                bestSSE = totalSSE;
                bestB1 = b1;
                bestB2 = b2;
                bestSeg1 = seg1;
                bestSeg2 = seg2;
                bestSeg3 = seg3;
            }
        }
    }
    
    return {bestB1, bestB2, bestSSE, bestSeg1, bestSeg2, bestSeg3};
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
    
    // Build stdev curve from position statistics
    std::vector<double> stdevCurve(readLength, std::nan(""));
    uint32_t maxPos = 0;
    
    for (const auto& kv : positionStats_) {
        uint32_t pos = kv.first;
        if (pos >= readLength) {
            continue;
        }
        const auto& stats = kv.second;
        
        // Use T→C rate stdev as the signal
        double tcStddev = stats.stddevTcRate();
        stdevCurve[pos] = tcStddev;
        
        if (pos > maxPos) {
            maxPos = pos;
        }
    }
    
    if (maxPos == 0) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Truncate to actual read length observed
    stdevCurve.resize(maxPos + 1);
    
    // Interpolate missing values
    interpolateMissing(stdevCurve);
    
    // Smooth with median window
    std::vector<double> smoothed = smoothMedian(stdevCurve, smoothWindow_);
    result.smoothedCurve = smoothed;
    
    // Segmented regression with 2 breakpoints
    auto [b1, b2, sse, seg1, seg2, seg3] = segmentedRegression(smoothed, minSegLen_, maxTrim_);
    
    if (b1 == 0 && b2 == smoothed.size() - 1) {
        // No valid segmentation found
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Set trim values
    // b1 = first position of middle segment (trim5p = b1)
    // b2 = last position of middle segment (trim3p = readLen - 1 - b2)
    uint32_t n = static_cast<uint32_t>(smoothed.size());
    result.trim5p = std::min(b1, maxTrim_);
    result.trim3p = std::min(n - 1 - b2, maxTrim_);
    
    // Store breakpoints for backward compatibility and QC
    result.kneeBin5p = b1;
    result.kneeBin3p = b2;
    result.totalSSE = sse;
    result.seg1 = seg1;
    result.seg2 = seg2;
    result.seg3 = seg3;
    
    result.mode = "auto_segmented";
    result.success = true;
    
    return result;
}
