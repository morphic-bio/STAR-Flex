#include "SnpMaskBuild.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Transcriptome.h"
#include "ErrorWarning.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/tbx.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <memory>
#include <limits>

SnpMaskBuild::SnpMaskBuild(Parameters& P, const Genome& genome)
    : P_(P), genome_(genome) {
    emModel_.reset(new SlamSnpEM(P.quant.slamSnpMask.maxIter, P.quant.slamSnpMask.convergeRelLL));
}

bool SnpMaskBuild::parseFofn(const std::string& fofnPath,
                              std::vector<std::pair<std::string, std::string>>& fastqPairs,
                              std::string* err) {
    fastqPairs.clear();
    std::ifstream in(fofnPath.c_str());
    if (!in.good()) {
        if (err) *err = "Cannot open FOFN file: " + fofnPath;
        return false;
    }
    
    std::string line;
    size_t lineNum = 0;
    while (std::getline(in, line)) {
        ++lineNum;
        
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) {
            continue;  // Blank line
        }
        size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start, end - start + 1);
        
        // Skip comment lines
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Parse tab-separated or space-separated
        size_t tabPos = line.find('\t');
        size_t spacePos = line.find(' ');
        
        if (tabPos != std::string::npos) {
            // Paired-end: R1<TAB>R2
            std::string r1 = line.substr(0, tabPos);
            std::string r2 = line.substr(tabPos + 1);
            
            // Trim whitespace from each
            size_t r1_end = r1.find_last_not_of(" \t");
            if (r1_end != std::string::npos) {
                r1 = r1.substr(0, r1_end + 1);
            }
            size_t r2_start = r2.find_first_not_of(" \t");
            if (r2_start != std::string::npos) {
                r2 = r2.substr(r2_start);
            }
            size_t r2_end = r2.find_last_not_of(" \t");
            if (r2_end != std::string::npos) {
                r2 = r2.substr(0, r2_end + 1);
            }
            
            if (!r1.empty() && !r2.empty()) {
                fastqPairs.push_back({r1, r2});
            } else {
                if (err) {
                    std::ostringstream oss;
                    oss << "FOFN line " << lineNum << ": invalid PE format (expected R1<TAB>R2)";
                    *err = oss.str();
                }
                return false;
            }
        } else if (spacePos != std::string::npos) {
            // Try space-separated PE
            std::istringstream iss(line);
            std::string r1, r2;
            if (iss >> r1 >> r2) {
                fastqPairs.push_back({r1, r2});
            } else {
                // Single-end: just one path
                fastqPairs.push_back({line, ""});
            }
        } else {
            // Single-end: one path
            fastqPairs.push_back({line, ""});
        }
    }
    
    if (fastqPairs.empty()) {
        if (err) *err = "FOFN file contains no valid FASTQ paths";
        return false;
    }
    
    return true;
}

void SnpMaskBuild::recordObservation(uint64_t pos, bool isMismatch) {
    uint32_t& entry = positionCounts_[pos];
    uint32_t cov = entry >> 16;
    uint32_t mis = entry & 0xFFFF;
    
    // Cap at 65535
    if (cov < 0xFFFF) {
        ++cov;
    } else {
        // Overflow: increment overflow counter (will be tracked separately)
        // Still increment mismatch if applicable, but don't overflow coverage
    }
    
    if (isMismatch && mis < 0xFFFF) {
        ++mis;
    }
    
    entry = (cov << 16) | mis;
    touchedPositions_.insert(pos);
}

SnpHistogram SnpMaskBuild::buildHistogram() const {
    SnpHistogram hist;
    
    for (const auto& kv : positionCounts_) {
        uint32_t packed = kv.second;
        uint32_t n = packed >> 16;
        // Note: k = packed & 0xFFFF (mismatches) is encoded in packed
        
        // Apply minCov filter
        if (n >= P_.quant.slamSnpMask.minCov) {
            hist[packed] += 1;  // One site per position
        }
    }
    
    return hist;
}

void SnpMaskBuild::filterAndComputePosteriors() {
    maskedPositions_.clear();
    
    if (!emModel_ || emResult_.iterations == 0) {
        return;  // EM not fitted
    }
    
    for (const auto& kv : positionCounts_) {
        uint64_t pos = kv.first;
        uint32_t packed = kv.second;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        
        // Apply filters
        if (n < P_.quant.slamSnpMask.minCov) {
            continue;
        }
        if (k < P_.quant.slamSnpMask.minAlt) {
            continue;
        }
        
        // Compute posterior
        double post_snp = emModel_->posterior(n, k);
        
        if (post_snp >= P_.quant.slamSnpMask.posterior) {
            maskedPositions_.insert(pos);
        }
    }
}

char SnpMaskBuild::getRefBase(uint64_t pos) const {
    if (pos >= genome_.nGenome) {
        return 'N';
    }
    uint8_t base = genome_.G[pos];
    const char bases[] = {'A', 'C', 'G', 'T'};
    if (base < 4) {
        return bases[base];
    }
    return 'N';
}

char SnpMaskBuild::getAltBase(uint64_t /* pos */) const {
    // For now, return 'C' as placeholder (most common T->C mismatch)
    // In a full implementation, we'd track the actual alt base from reads
    return 'C';
}

std::string SnpMaskBuild::getComponent(uint32_t n, uint32_t k) const {
    if (!emModel_ || emResult_.iterations == 0) {
        return "ERR";
    }
    
    // Compute responsibilities to determine component
    // Note: log_binom_pmf is in slam_snp_em.cpp, but we'll compute directly here
    auto log_binom = [](uint32_t n, uint32_t k, double p) -> double {
        if (n == 0 || k > n) return -std::numeric_limits<double>::infinity();
        if (p <= 0.0) return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
        if (p >= 1.0) return (k == n) ? 0.0 : -std::numeric_limits<double>::infinity();
        double log_n_choose_k = (k > 0 && k < n) ?
            std::lgamma(static_cast<double>(n + 1)) -
            std::lgamma(static_cast<double>(k + 1)) -
            std::lgamma(static_cast<double>(n - k + 1)) : 0.0;
        return log_n_choose_k + static_cast<double>(k) * std::log(p) +
               static_cast<double>(n - k) * std::log1p(-p);
    };
    
    double log_prob_ERR = log_binom(n, k, emModel_->get_p_ERR());
    double log_prob_HET = log_binom(n, k, 0.5);
    double log_prob_HOM = log_binom(n, k, 0.95);
    
    double log_weight_ERR = std::log(emModel_->get_pi_ERR()) + log_prob_ERR;
    double log_weight_HET = std::log(emModel_->get_pi_HET()) + log_prob_HET;
    double log_weight_HOM = std::log(emModel_->get_pi_HOM()) + log_prob_HOM;
    
    double max_weight = std::max({log_weight_ERR, log_weight_HET, log_weight_HOM});
    
    if (max_weight == log_weight_HET) {
        return "HET";
    } else if (max_weight == log_weight_HOM) {
        return "HOM";
    } else {
        return "ERR";
    }
}

bool SnpMaskBuild::buildMask(const std::vector<std::pair<std::string, std::string>>& /* fastqPairs */,
                             Transcriptome* /* transcriptome */,
                             SnpMaskBuildStats* stats,
                             std::string* err) {
    // Data should already be extracted via extractFromSlamQuant() before this is called
    // (done in STAR.cpp after alignment completes)
    
    if (positionCounts_.empty()) {
        if (err) *err = "No SNP observations collected - extractFromSlamQuant() must be called first";
        return false;
    }
    
    // Build histogram
    SnpHistogram histogram = buildHistogram();
    
    if (histogram.empty()) {
        if (err) *err = "No candidate sites after filtering (minCov=" + 
                        std::to_string(P_.quant.slamSnpMask.minCov) + ")";
        return false;
    }
    
    // Fit EM model
    emResult_ = emModel_->fit(histogram);
    
    if (!emResult_.converged && emResult_.iterations >= P_.quant.slamSnpMask.maxIter) {
        // Warning but not fatal - use the best fit we have
        if (err) {
            *err = "EM did not converge after " + std::to_string(emResult_.iterations) + " iterations";
        }
    }
    
    // Filter and compute posteriors
    filterAndComputePosteriors();
    
    // Compute statistics
    if (stats) {
        stats->totalSites = touchedPositions_.size();
        
        // Count actual candidate sites (not just unique (n,k) pairs)
        uint64_t candidateCount = 0;
        for (const auto& kv : histogram) {
            candidateCount += kv.second;  // Each histogram entry has count of sites
        }
        stats->candidateSites = candidateCount;
        stats->maskedSites = maskedPositions_.size();
        
        // Compute global baseline before/after masking
        uint64_t total_k_before = 0;
        uint64_t total_n_before = 0;
        uint64_t total_k_after = 0;
        uint64_t total_n_after = 0;
        
        for (const auto& kv : positionCounts_) {
            uint32_t packed = kv.second;
            uint32_t n = packed >> 16;
            uint32_t k = packed & 0xFFFF;
            
            total_k_before += k;
            total_n_before += n;
            
            if (maskedPositions_.count(kv.first) == 0) {
                total_k_after += k;
                total_n_after += n;
            }
        }
        
        stats->globalBaselineBefore = (total_n_before > 0) ? 
            (static_cast<double>(total_k_before) / static_cast<double>(total_n_before)) : 0.0;
        stats->globalBaselineAfter = (total_n_after > 0) ? 
            (static_cast<double>(total_k_after) / static_cast<double>(total_n_after)) : 0.0;
        
        // Count coverage overflows (sites with coverage >= 65535)
        stats->coverageOverflowCount = 0;
        for (const auto& kv : positionCounts_) {
            uint32_t packed = kv.second;
            uint32_t n = packed >> 16;
            if (n >= 0xFFFF) {
                stats->coverageOverflowCount++;
            }
        }
    }
    
    return true;
}

void SnpMaskBuild::extractFromSlamQuant(SlamQuant* slamQuant) {
    if (!slamQuant) {
        return;
    }
    
    // Extract SNP mask data from SlamQuant
    positionCounts_ = slamQuant->getSnpMaskData();
    
    // Build touched positions set
    touchedPositions_.clear();
    for (const auto& kv : positionCounts_) {
        touchedPositions_.insert(kv.first);
    }
}

bool SnpMaskBuild::writeBed(const std::string& bedPath, const Genome& genome, std::string* err) {
    // Build list of masked positions with metadata
    struct BedEntry {
        uint32_t chrIdx;
        uint64_t pos;
        uint32_t n;
        uint32_t k;
        double f;
        double post_snp;
        std::string comp;
    };
    
    std::vector<BedEntry> entries;
    entries.reserve(maskedPositions_.size());
    
    // Helper lambda for binary search to find chromosome index
    auto findChrIdx = [&genome](uint64_t pos) -> uint32_t {
        // Binary search in chrStart
        const auto& chrStart = genome.chrStart;
        if (chrStart.empty()) return 0;
        
        // Find the largest chrStart[i] <= pos
        auto it = std::upper_bound(chrStart.begin(), chrStart.end(), pos);
        if (it == chrStart.begin()) {
            return 0;
        }
        return static_cast<uint32_t>(std::distance(chrStart.begin(), it) - 1);
    };
    
    for (uint64_t pos : maskedPositions_) {
        uint32_t packed = positionCounts_.at(pos);
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        double f = (n > 0) ? (static_cast<double>(k) / static_cast<double>(n)) : 0.0;
        double post_snp = emModel_->posterior(n, k);
        
        // Find chromosome index using binary search
        uint32_t chrIdx = findChrIdx(pos);
        
        BedEntry entry;
        entry.chrIdx = chrIdx;
        entry.pos = pos;
        entry.n = n;
        entry.k = k;
        entry.f = f;
        entry.post_snp = post_snp;
        entry.comp = getComponent(n, k);
        entries.push_back(entry);
    }
    
    // Sort by chromosome and position
    std::sort(entries.begin(), entries.end(), [&](const BedEntry& a, const BedEntry& b) {
        if (a.chrIdx != b.chrIdx) {
            return a.chrIdx < b.chrIdx;
        }
        return a.pos < b.pos;
    });
    
    // Write BED file (bgzip-compressed)
    BGZF* bgzf = bgzf_open(bedPath.c_str(), "w");
    if (!bgzf) {
        if (err) *err = "Cannot open BED output file: " + bedPath;
        return false;
    }
    
    // Write header comment
    std::string header = "#chrom\tstart\tend\tref\talt\tn\tk\tf\tpost_snp\tcomp\n";
    bgzf_write(bgzf, header.c_str(), header.size());
    
    for (const BedEntry& entry : entries) {
        if (entry.chrIdx >= genome.chrName.size()) {
            continue;
        }
        
        uint64_t chrStart = genome.chrStart[entry.chrIdx];
        uint64_t pos0 = entry.pos - chrStart;
        char ref = getRefBase(entry.pos);
        char alt = getAltBase(entry.pos);
        
        std::ostringstream oss;
        oss << genome.chrName[entry.chrIdx] << "\t"
            << pos0 << "\t"
            << (pos0 + 1) << "\t"
            << ref << "\t"
            << alt << "\t"
            << entry.n << "\t"
            << entry.k << "\t"
            << std::fixed << std::setprecision(6) << entry.f << "\t"
            << std::fixed << std::setprecision(6) << entry.post_snp << "\t"
            << entry.comp << "\n";
        
        std::string line = oss.str();
        bgzf_write(bgzf, line.c_str(), line.size());
    }
    
    bgzf_close(bgzf);
    
    // Build tabix index
    int ret = tbx_index_build(bedPath.c_str(), 0, &tbx_conf_bed);
    if (ret != 0) {
        if (err) *err = "Failed to build tabix index for: " + bedPath;
        return false;
    }
    
    return true;
}

bool SnpMaskBuild::writeSummary(const std::string& summaryPath, const SnpMaskBuildStats& stats, std::string* err) {
    std::ofstream out(summaryPath.c_str());
    if (!out.good()) {
        if (err) *err = "Cannot open summary output file: " + summaryPath;
        return false;
    }
    
    out << "metric\tvalue\n";
    out << "p_ERR\t" << std::fixed << std::setprecision(8) << emResult_.p_ERR << "\n";
    out << "pi_ERR\t" << std::fixed << std::setprecision(8) << emResult_.pi_ERR << "\n";
    out << "pi_HET\t" << std::fixed << std::setprecision(8) << emResult_.pi_HET << "\n";
    out << "pi_HOM\t" << std::fixed << std::setprecision(8) << emResult_.pi_HOM << "\n";
    out << "candidates_total\t" << stats.totalSites << "\n";
    out << "candidates_passing\t" << stats.candidateSites << "\n";
    out << "masked_sites\t" << stats.maskedSites << "\n";
    out << "global_baseline_before\t" << std::fixed << std::setprecision(8) << stats.globalBaselineBefore << "\n";
    out << "global_baseline_after\t" << std::fixed << std::setprecision(8) << stats.globalBaselineAfter << "\n";
    out << "iterations\t" << emResult_.iterations << "\n";
    out << "final_log_likelihood\t" << std::fixed << std::setprecision(8) << emResult_.final_log_likelihood << "\n";
    out << "coverage_overflow_count\t" << stats.coverageOverflowCount << "\n";
    
    return true;
}
