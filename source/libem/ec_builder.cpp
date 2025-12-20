#include "ec_builder.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

// Check if alignment is compatible with library format (for single-end or orphaned reads)
// Ported from Salmon's compatibleHit() function
bool compatibleHit(const LibraryFormat& expected, bool isForward, MateStatus ms) {
    auto expectedStrand = expected.strandedness;
    auto expectedType = expected.type;
    
    switch (ms) {
    case MateStatus::SINGLE_END:
        if (isForward) { // U, SF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::S);
        } else { // U, SR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::A);
        }
        break;
        
    case MateStatus::PAIRED_END_LEFT:
        if (expected.orientation == ReadOrientation::SAME) {
            return (expectedStrand == ReadStrandedness::U ||
                    (expectedStrand == ReadStrandedness::S && isForward) ||
                    (expectedStrand == ReadStrandedness::A && !isForward));
        } else if (isForward) { // IU, ISF, OU, OSF, MU, MSF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::SA);
        } else { // IU, ISR, OU, OSR, MU, MSR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::AS);
        }
        break;
        
    case MateStatus::PAIRED_END_RIGHT:
        if (expected.orientation == ReadOrientation::SAME) {
            return (expectedStrand == ReadStrandedness::U ||
                    (expectedStrand == ReadStrandedness::S && isForward) ||
                    (expectedStrand == ReadStrandedness::A && !isForward));
        } else if (isForward) { // IU, ISR, OU, OSR, MU, MSR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::AS);
        } else { // IU, ISF, OU, OSF, MU, MSF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::SA);
        }
        break;
        
    default:
        return false;
    }
}

// Check if paired-end alignment is compatible (check observed vs expected format)
bool compatibleHit(const LibraryFormat& expected, const LibraryFormat& observed) {
    if (observed.type != ReadType::PAIRED_END) {
        return false;
    }
    
    auto es = expected.strandedness;
    auto eo = expected.orientation;
    auto os = observed.strandedness;
    auto oo = observed.orientation;
    
    // If orientations are different, incompatible
    if (eo != oo) {
        return false;
    } else {
        // Orientations match - check strandedness
        return (es == ReadStrandedness::U || es == os);
    }
}

// Main compatibility check function (matching Salmon's isCompatible)
bool isCompatible(const RawAlignment& aln, const LibraryFormat& expected_format, bool isForward) {
    // For unstranded libraries, all alignments are compatible
    if (expected_format.strandedness == ReadStrandedness::U) {
        return true;
    }
    
    // For paired-end paired reads, we'd need both mates' orientations to determine observed format
    // For now, assume compatible (proper pairs are usually compatible)
    // TODO: Enhance to check observed format from both mates when available
    if (aln.mate_status == MateStatus::PAIRED_END_PAIRED) {
        // Proper pairs are typically compatible - could enhance with observed format check
        return true;
    } else {
        // Single-end or orphan - use the compatibility check
        return compatibleHit(expected_format, isForward, aln.mate_status);
    }
}

// Compute auxProb for each transcript in a read's alignments
// Matches Salmon's SalmonQuantifyAlignments.cpp alignment-mode calculation
ReadMapping computeAuxProbs(
    const std::vector<RawAlignment>& alignments,
    const ECBuilderParams& params,
    bool enable_trace
) {
    ReadMapping mapping;
    mapping.aux_denom = LOG_0;  // -infinity in log space
    
    if (alignments.empty()) {
        return mapping;
    }
    
    // Step 1: Find bestAS across all alignments for this read
    // For paired-end: AS is sum of both mates; for single-end: AS tag value
    int32_t best_as = std::numeric_limits<int32_t>::min();
    for (const auto& aln : alignments) {
        int32_t aln_as = aln.score;  // Already contains AS (sum for paired, single AS for single)
        if (aln_as > best_as) {
            best_as = aln_as;
        }
    }
    mapping.best_as = best_as;
    
    // Step 2: Compute auxProb for each alignment
    for (size_t aln_idx = 0; aln_idx < alignments.size(); ++aln_idx) {
        const auto& aln = alignments[aln_idx];
        AlignmentTrace trace;
        trace.transcript_id = aln.transcript_id;
        trace.as_tag = aln.score;
        trace.best_as = best_as;
        
        // Skip decoys (if needed - for now, process all)
        if (aln.is_decoy) {
            // Could skip decoys here if desired
        }
        
        // Check compatibility
        bool is_compat = isCompatible(aln, params.lib_format, aln.is_forward);
        trace.dropped_incompat = false;
        
        // 3. logAlignCompatProb: Compatibility probability (compute before skipping)
        double log_compat_prob = 0.0;
        if (!is_compat) {
            if (params.incompat_prior > 0.0) {
                // Use log(incompatPrior) as penalty
                log_compat_prob = std::log(params.incompat_prior);
            } else {
                // If incompatPrior == 0.0, use LOG_0 (will cause skip)
                log_compat_prob = LOG_0;
            }
        }
        trace.log_compat_prob = log_compat_prob;
        
        // If ignore_incompat=true and not compatible, skip this alignment
        if (params.ignore_incompat && !is_compat) {
            trace.dropped_incompat = true;
            mapping.num_dropped_incompat++;
            if (enable_trace) {
                mapping.trace_info.push_back(trace);
            }
            continue;
        }
        
        // 1. logFragProb: Fragment length probability
        double log_frag_prob = 0.0;
        bool is_orphan = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                         aln.mate_status == MateStatus::PAIRED_END_RIGHT);
        trace.is_orphan = is_orphan;
        
        if (is_orphan) {
            // For orphans: use LOG_EPSILON (orphanProb = LOG_EPSILON when discardOrphansAln=false)
            log_frag_prob = LOG_EPSILON;
        } else if (params.use_frag_len_dist) {
            // If FLD enabled and paired/inward, use PMF (fragment length distribution)
            // For parity with --noFragLengthDist, this branch is skipped
            log_frag_prob = aln.log_frag_prob;  // Would be computed from PMF
        } else {
            // For parity with --noFragLengthDist: set to 0.0
            log_frag_prob = 0.0;
        }
        trace.log_frag_prob = log_frag_prob;
        
        // 2. errLike: Error model likelihood
        // If pre-computed from CIGAR (has_err_like), use that
        // Otherwise fall back to AS-based: errLike = -scoreExp * (bestAS - alnAS)
        int32_t aln_as = aln.score;
        double err_like = 0.0;
        if (params.use_error_model && aln.has_err_like) {
            // Use pre-computed CIGAR-based error likelihood
            err_like = aln.err_like;
        } else if (params.use_as_without_cigar) {
            // AS-based likelihood (RapMap/Pufferfish style)
            double score_diff = static_cast<double>(best_as - aln_as);
            err_like = -params.score_exp * score_diff;
        }
        trace.err_like = err_like;
        
        // Combine: auxProb = logFragProb + errLike + logAlignCompatProb
        double aux_prob = log_frag_prob + err_like + log_compat_prob;
        trace.aux_prob = aux_prob;
        
        // Skip if aux_prob is LOG_0 (incompatible and incompatPrior == 0.0)
        if (aux_prob == LOG_0) {
            if (enable_trace) {
                mapping.trace_info.push_back(trace);
            }
            continue;
        }
        
        // Match Salmon's filtering checks:
        // 1. transcriptLogCount != LOG_0 (transcript mass check)
        //    In initial round, all transcripts should have prior mass, so this always passes.
        //    We don't track transcript mass, so we skip this check.
        //    TODO: If we add transcript mass tracking, add: if (transcriptLogCount == LOG_0) continue;
        
        // 2. startPosProb != LOG_0 (start position probability)
        //    With --noLengthCorrection: startPosProb = -logRefLength = -1.0 (not LOG_0)
        //    So this check always passes. We skip it since we don't use length correction.
        //    TODO: If we add length correction, add: if (startPosProb == LOG_0) continue;
        
        mapping.transcript_ids.push_back(aln.transcript_id);
        mapping.aux_probs.push_back(aux_prob);
        mapping.alignment_indices.push_back(aln_idx);
        mapping.aux_denom = logAdd(mapping.aux_denom, aux_prob);  // log-sum-exp
        
        if (enable_trace) {
            mapping.trace_info.push_back(trace);
        }
    }
    
    // Step 3: Normalize weights: weight = exp(auxProb - auxDenom)
    // Note: Salmon keeps duplicate transcript IDs for multi-mapping reads
    // (each alignment location is a separate entry in the EC)
    if (mapping.aux_denom != LOG_0) {
        for (size_t i = 0; i < mapping.aux_probs.size(); ++i) {
            double normalized_weight = std::exp(mapping.aux_probs[i] - mapping.aux_denom);
            mapping.aux_probs[i] = normalized_weight;
        }
    }
    
    return mapping;
}

// Apply range factorization by appending bin IDs to transcript list
// Matches Salmon's SalmonQuantify.cpp lines 578-586
void applyRangeFactorization(
    std::vector<uint32_t>& txp_ids,
    const std::vector<double>& aux_probs,
    uint32_t range_factorization_bins
) {
    if (range_factorization_bins == 0) return;
    
    int32_t txps_size = static_cast<int32_t>(txp_ids.size());
    // rangeCount = sqrt(n) + bins (Salmon's formula)
    int32_t range_count = static_cast<int32_t>(std::sqrt(static_cast<double>(txps_size))) + static_cast<int32_t>(range_factorization_bins);
    
    // Append range bin numbers as pseudo-transcript IDs
    for (int32_t i = 0; i < txps_size; i++) {
        int32_t range_number = static_cast<int32_t>(aux_probs[i] * range_count);
        txp_ids.push_back(static_cast<uint32_t>(range_number));
    }
}

// Sort transcripts by conditional probability for rank-based ECs
// Matches Salmon's SalmonQuantify.cpp lines 557-575
void applyRankEqClasses(
    std::vector<uint32_t>& txp_ids,
    std::vector<double>& aux_probs
) {
    if (txp_ids.size() <= 1) return;
    
    std::vector<size_t> inds(txp_ids.size());
    std::iota(inds.begin(), inds.end(), 0);
    
    // Sort indices by aux_probs (ascending order)
    std::sort(inds.begin(), inds.end(), [&aux_probs](size_t i, size_t j) {
        return aux_probs[i] < aux_probs[j];
    });
    
    std::vector<uint32_t> txp_ids_new(txp_ids.size());
    std::vector<double> aux_probs_new(aux_probs.size());
    for (size_t r = 0; r < inds.size(); ++r) {
        txp_ids_new[r] = txp_ids[inds[r]];
        aux_probs_new[r] = aux_probs[inds[r]];
    }
    
    std::swap(txp_ids, txp_ids_new);
    std::swap(aux_probs, aux_probs_new);
}

// Build equivalence classes from filtered alignments
// This is the main entry point that combines all the steps
ECTable buildEquivalenceClasses(
    const std::vector<std::vector<RawAlignment>>& read_alignments,
    const ECBuilderParams& params,
    size_t num_transcripts
) {
    ECTable ec_table;
    ec_table.n_transcripts = num_transcripts;
    
    // Map from EC label (order-sensitive) to EC index
    // Salmon treats TranscriptGroup ordering as significant.
    std::unordered_map<std::string, size_t> ec_map;
    
    for (const auto& alignments : read_alignments) {
        if (alignments.empty()) continue;
        
        // Compute auxProbs for this read (tracing disabled in EC building)
        ReadMapping mapping = computeAuxProbs(alignments, params, false);
        
        if (mapping.transcript_ids.empty()) continue;
        
        // Apply rank-based ECs if enabled
        if (params.use_rank_eq_classes && mapping.transcript_ids.size() > 1) {
            applyRankEqClasses(mapping.transcript_ids, mapping.aux_probs);
        }
        
        // Apply range factorization if enabled
        if (params.use_range_factorization) {
            applyRangeFactorization(mapping.transcript_ids, mapping.aux_probs, params.range_factorization_bins);
        }
        
        std::string ec_key;
        for (size_t i = 0; i < mapping.transcript_ids.size(); ++i) {
            if (i > 0) ec_key += ",";
            ec_key += std::to_string(mapping.transcript_ids[i]);
        }
        
        // Find or create EC
        auto it = ec_map.find(ec_key);
        if (it == ec_map.end()) {
            // Create new EC
            EC ec;
            ec.transcript_ids = mapping.transcript_ids;
            ec.weights = mapping.aux_probs;
            ec.count = 1.0;
            ec_table.ecs.push_back(ec);
            ec_map[ec_key] = ec_table.ecs.size() - 1;
        } else {
            // Increment count and accumulate weights (Salmon sums weights per EC)
            EC& ec = ec_table.ecs[it->second];
            ec.count += 1.0;
            if (ec.weights.size() == mapping.aux_probs.size()) {
                for (size_t i = 0; i < ec.weights.size(); ++i) {
                    ec.weights[i] += mapping.aux_probs[i];
                }
            }
        }
    }
    
    ec_table.n_ecs = ec_table.ecs.size();
    return ec_table;
}
