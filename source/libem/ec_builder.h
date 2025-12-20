#ifndef EC_BUILDER_H
#define EC_BUILDER_H

#include "em_types.h"
#include "alignment_filter.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>

// Log-space constants
constexpr double LOG_0 = -std::numeric_limits<double>::infinity();
constexpr double LOG_1 = 0.0;
// LOG_EPSILON: Very small log probability used for orphan penalties (Salmon default)
// log(0.375e-10) ≈ -24.006680182952184
constexpr double LOG_EPSILON = -24.006680182952184;

// Log-sum-exp: log(exp(a) + exp(b)) = a + log(1 + exp(b - a))
// Numerically stable when a >= b
inline double logAdd(double a, double b) {
    if (a == LOG_0) return b;
    if (b == LOG_0) return a;
    if (a > b) {
        return a + std::log(1.0 + std::exp(b - a));
    } else {
        return b + std::log(1.0 + std::exp(a - b));
    }
}

// Library format components (matching Salmon)
enum class ReadType : uint8_t { SINGLE_END = 0, PAIRED_END = 1 };
enum class ReadOrientation : uint8_t {
    SAME = 0,      // M (matching orientation)
    AWAY = 1,      // O (outward)
    TOWARD = 2,    // I (inward)
    NONE = 3       // For single-end
};
enum class ReadStrandedness : uint8_t {
    SA = 0,        // Stranded, antisense
    AS = 1,         // Antisense
    S = 2,          // Stranded, sense
    A = 3,          // Sense
    U = 4           // Unstranded
};

// Library format for compatibility checking (matching Salmon's LibraryFormat)
struct LibraryFormat {
    ReadType type;
    ReadOrientation orientation;
    ReadStrandedness strandedness;
    
    LibraryFormat() : type(ReadType::PAIRED_END), orientation(ReadOrientation::TOWARD), strandedness(ReadStrandedness::U) {}
    LibraryFormat(ReadType t, ReadOrientation o, ReadStrandedness s) : type(t), orientation(o), strandedness(s) {}
    
    // Common presets
    static LibraryFormat IU() { return {ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U}; }  // Inward Unstranded
    static LibraryFormat ISF() { return {ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA}; } // Inward Stranded Forward
    static LibraryFormat ISR() { return {ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS}; } // Inward Stranded Reverse
    static LibraryFormat U() { return {ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U}; }      // Single-end Unstranded
    
    bool operator==(const LibraryFormat& other) const {
        return type == other.type && orientation == other.orientation && strandedness == other.strandedness;
    }
    
    // Generate a unique ID for hashing/comparison
    uint8_t typeId() const {
        return static_cast<uint8_t>(type) | 
               (static_cast<uint8_t>(orientation) << 1) |
               (static_cast<uint8_t>(strandedness) << 3);
    }
    
    // Create LibraryFormat from ID
    static LibraryFormat formatFromID(uint8_t id) {
        ReadType rt = static_cast<ReadType>(id & 0x01);
        ReadOrientation ro = static_cast<ReadOrientation>((id >> 1) & 0x3);
        ReadStrandedness rs = static_cast<ReadStrandedness>((id >> 3) & 0x7);
        return LibraryFormat(rt, ro, rs);
    }
};

// Configuration for EC building
struct ECBuilderParams {
    bool use_range_factorization = false;
    uint32_t range_factorization_bins = 4;   // Salmon default
    bool use_rank_eq_classes = false;        // Sort by conditional prob
    bool use_frag_len_dist = false;          // Fragment length distribution (default false for parity)
    bool use_error_model = false;            // Error model (default false for parity)
    bool use_as_without_cigar = false;       // AS-based likelihood (RapMap/Pufferfish)
    
    // Alignment-mode auxProb parameters (matching Salmon)
    double score_exp = 1.0;                  // Score exponent for AS-based likelihood (Salmon default)
    double incompat_prior = 0.0;             // Incompatibility prior (linear space; 0.0 → ignore_incompat=true)
    bool ignore_incompat = true;             // If true, skip incompatible alignments (set when incompat_prior == 0.0)
    LibraryFormat lib_format = LibraryFormat::IU();  // Library format for compat check (default: IU = inward unstranded)
    bool discard_orphans = false;            // If true, discard orphan reads (default false)
};

// Per-alignment trace information
struct AlignmentTrace {
    uint32_t transcript_id;
    int32_t as_tag;
    int32_t best_as;
    double log_frag_prob;
    double err_like;
    double log_compat_prob;
    bool dropped_incompat;
    bool is_orphan;
    double aux_prob;           // Before normalization
    double weight;             // After normalization
};

// Intermediate representation before EC aggregation
struct ReadMapping {
    std::vector<uint32_t> transcript_ids;
    std::vector<double> aux_probs;     // Log-space initially, then normalized
    std::vector<size_t> alignment_indices;   // Indices into the original alignment list
    double aux_denom;                   // log-sum for normalization
    
    // Trace information (populated if tracing enabled)
    std::vector<AlignmentTrace> trace_info;
    int32_t best_as;                    // Best AS across all alignments
    size_t num_dropped_incompat;        // Count of dropped incompatible alignments
    
    ReadMapping() : aux_denom(LOG_0), best_as(0), num_dropped_incompat(0) {}
};

// Compatibility checking functions (ported from Salmon)
bool compatibleHit(const LibraryFormat& expected, bool isForward, MateStatus ms);
bool compatibleHit(const LibraryFormat& expected, const LibraryFormat& observed);
bool isCompatible(const RawAlignment& aln, const LibraryFormat& expected_format, bool isForward);

// Compute auxProb for each transcript in a read's alignments
// Matches Salmon's SalmonQuantify.cpp lines 506-530
ReadMapping computeAuxProbs(
    const std::vector<RawAlignment>& alignments,
    const ECBuilderParams& params,
    bool enable_trace = false  // If true, populate trace_info
);

// Apply range factorization by appending bin IDs to transcript list
// Matches Salmon's SalmonQuantify.cpp lines 578-586
void applyRangeFactorization(
    std::vector<uint32_t>& txp_ids,
    const std::vector<double>& aux_probs,
    uint32_t range_factorization_bins
);

// Sort transcripts by conditional probability for rank-based ECs
// Matches Salmon's SalmonQuantify.cpp lines 557-575
void applyRankEqClasses(
    std::vector<uint32_t>& txp_ids,
    std::vector<double>& aux_probs
);

// Build equivalence classes from filtered alignments
// This is the main entry point that combines all the steps
ECTable buildEquivalenceClasses(
    const std::vector<std::vector<RawAlignment>>& read_alignments,
    const ECBuilderParams& params,
    size_t num_transcripts
);

#endif // EC_BUILDER_H
