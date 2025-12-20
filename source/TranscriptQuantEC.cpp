#include "TranscriptQuantEC.h"
#include "Transcriptome.h"
#include "ec_builder.h"
#include <algorithm>
#include <cmath>

TranscriptQuantEC::TranscriptQuantEC(uint num_transcripts) 
    : num_transcripts_(num_transcripts) {
    ecTable_.n_transcripts = static_cast<size_t>(num_transcripts);
    ecTable_.n_ecs = 0;
}

RawAlignment TranscriptQuantEC::transcriptToRawAlignment(Transcript* tr, uint32_t transcript_id) {
    RawAlignment aln;
    aln.transcript_id = transcript_id;
    aln.pos = tr->exons[0][EX_G];  // Transcriptomic start position
    aln.score = tr->maxScore;
    aln.est_aln_prob = 1.0;  // Will be normalized later
    aln.log_frag_prob = 0.0;  // Not using fragment length dist for now
    aln.log_compat_prob = 0.0;  // Compatibility already checked
    aln.err_like = 0.0;
    aln.has_err_like = false;
    aln.is_decoy = false;  // TODO: check if transcript is decoy
    aln.is_forward = (tr->Str == 0);
    
    // Determine mate status
    if (tr->readNmates == 1) {
        aln.mate_status = MateStatus::SINGLE_END;
    } else {
        // Check if this is a proper pair
        bool proper_pair = (tr->exons[0][EX_iFrag] != tr->exons[tr->nExons-1][EX_iFrag]);
        if (proper_pair) {
            aln.mate_status = MateStatus::PAIRED_END_PAIRED;
        } else {
            // Determine which mate
            if (tr->exons[0][EX_iFrag] == 0) {
                aln.mate_status = MateStatus::PAIRED_END_LEFT;
            } else {
                aln.mate_status = MateStatus::PAIRED_END_RIGHT;
            }
        }
    }
    
    return aln;
}

void TranscriptQuantEC::addReadAlignments(Transcript** trMult, uint nTr, const std::vector<double>& auxProbs) {
    if (nTr == 0) return;
    
    current_read_alignments_.clear();
    current_read_alignments_.reserve(nTr);
    
    // Convert Transcript objects to RawAlignment
    for (uint i = 0; i < nTr; ++i) {
        Transcript* tr = trMult[i];
        uint transcript_id = tr->Chr;  // Chr field contains transcript index
        
        if (transcript_id >= num_transcripts_) {
            continue;  // Skip invalid transcript IDs
        }
        
        RawAlignment aln = transcriptToRawAlignment(tr, static_cast<uint32_t>(transcript_id));
        
        // Set alignment probability if provided
        if (i < auxProbs.size()) {
            aln.est_aln_prob = auxProbs[i];
        }
        
        current_read_alignments_.push_back(aln);
    }
    
    if (current_read_alignments_.empty()) return;
    
    // Build EC from this read's alignments
    // Use ECBuilderParams with defaults matching Salmon alignment mode
    ECBuilderParams params;
    params.use_range_factorization = false;
    params.use_rank_eq_classes = false;
    params.use_frag_len_dist = false;
    params.use_error_model = false;
    params.ignore_incompat = true;
    params.incompat_prior = 0.0;
    params.lib_format = LibraryFormat::IU();  // Default: Inward Unstranded
    
    // Compute auxProbs for this read
    ReadMapping mapping = computeAuxProbs(current_read_alignments_, params, false);
    
    if (mapping.transcript_ids.empty()) return;
    
    // Create EC signature (sorted transcript IDs + weights)
    ECSignature sig;
    sig.transcript_ids = mapping.transcript_ids;
    sig.weights.reserve(mapping.aux_probs.size());
    
    // Convert log-probs to linear weights (for EC signature)
    for (double log_prob : mapping.aux_probs) {
        sig.weights.push_back(std::exp(log_prob));
    }
    
    // Sort transcript IDs for canonical EC representation
    std::vector<size_t> sort_indices(sig.transcript_ids.size());
    for (size_t i = 0; i < sort_indices.size(); ++i) {
        sort_indices[i] = i;
    }
    std::sort(sort_indices.begin(), sort_indices.end(), 
        [&](size_t a, size_t b) {
            return sig.transcript_ids[a] < sig.transcript_ids[b];
        });
    
    // Reorder transcript IDs and weights
    std::vector<uint32_t> sorted_ids;
    std::vector<double> sorted_weights;
    sorted_ids.reserve(sig.transcript_ids.size());
    sorted_weights.reserve(sig.weights.size());
    for (size_t idx : sort_indices) {
        sorted_ids.push_back(sig.transcript_ids[idx]);
        sorted_weights.push_back(sig.weights[idx]);
    }
    sig.transcript_ids = sorted_ids;
    sig.weights = sorted_weights;
    
    // Find or create EC
    auto it = ec_signature_map_.find(sig);
    if (it != ec_signature_map_.end()) {
        // Existing EC: increment count
        size_t ec_idx = it->second;
        ecTable_.ecs[ec_idx].count += 1.0;
    } else {
        // New EC: create entry
        EC ec;
        ec.transcript_ids = sig.transcript_ids;
        ec.weights = sig.weights;  // Store weights for weighted format
        ec.count = 1.0;
        
        size_t ec_idx = ecTable_.ecs.size();
        ecTable_.ecs.push_back(ec);
        ecTable_.n_ecs++;
        ec_signature_map_[sig] = ec_idx;
    }
}

void TranscriptQuantEC::addReadAlignmentsSimple(const std::vector<uint32_t>& transcriptIds, const std::vector<double>& weights) {
    if (transcriptIds.empty()) return;
    
    // Create EC signature (sorted transcript IDs + weights)
    ECSignature sig;
    
    // Create sorted indices
    std::vector<size_t> sort_indices(transcriptIds.size());
    for (size_t i = 0; i < sort_indices.size(); ++i) {
        sort_indices[i] = i;
    }
    std::sort(sort_indices.begin(), sort_indices.end(), 
        [&](size_t a, size_t b) {
            return transcriptIds[a] < transcriptIds[b];
        });
    
    // Build sorted signature
    sig.transcript_ids.reserve(transcriptIds.size());
    sig.weights.reserve(weights.size());
    for (size_t idx : sort_indices) {
        sig.transcript_ids.push_back(transcriptIds[idx]);
        if (idx < weights.size()) {
            sig.weights.push_back(weights[idx]);
        } else {
            sig.weights.push_back(1.0 / transcriptIds.size());  // Default uniform
        }
    }
    
    // Find or create EC
    auto it = ec_signature_map_.find(sig);
    if (it != ec_signature_map_.end()) {
        // Existing EC: increment count
        size_t ec_idx = it->second;
        ecTable_.ecs[ec_idx].count += 1.0;
    } else {
        // New EC: create entry
        EC ec;
        ec.transcript_ids = sig.transcript_ids;
        ec.weights = sig.weights;
        ec.count = 1.0;
        
        size_t ec_idx = ecTable_.ecs.size();
        ecTable_.ecs.push_back(ec);
        ecTable_.n_ecs++;
        ec_signature_map_[sig] = ec_idx;
    }
}

void TranscriptQuantEC::addGCObservation(int32_t frag_start, int32_t frag_end, int32_t gc_pct, double weight) {
    observedGC_.inc(gc_pct, weight);
}

void TranscriptQuantEC::finalize() {
    // Normalize observed GC if needed
    // (Will be normalized later when computing bias ratios)
}

void TranscriptQuantEC::merge(const TranscriptQuantEC& other) {
    // Merge EC tables
    for (const EC& other_ec : other.ecTable_.ecs) {
        // Create signature for this EC
        ECSignature sig;
        sig.transcript_ids = other_ec.transcript_ids;
        sig.weights = other_ec.weights;
        
        // Sort for canonical representation
        std::vector<size_t> sort_indices(sig.transcript_ids.size());
        for (size_t i = 0; i < sort_indices.size(); ++i) {
            sort_indices[i] = i;
        }
        std::sort(sort_indices.begin(), sort_indices.end(),
            [&](size_t a, size_t b) {
                return sig.transcript_ids[a] < sig.transcript_ids[b];
            });
        
        std::vector<uint32_t> sorted_ids;
        std::vector<double> sorted_weights;
        sorted_ids.reserve(sig.transcript_ids.size());
        sorted_weights.reserve(sig.weights.size());
        for (size_t idx : sort_indices) {
            sorted_ids.push_back(sig.transcript_ids[idx]);
            sorted_weights.push_back(sig.weights[idx]);
        }
        sig.transcript_ids = sorted_ids;
        sig.weights = sorted_weights;
        
        // Find or create EC
        auto it = ec_signature_map_.find(sig);
        if (it != ec_signature_map_.end()) {
            // Existing EC: add counts
            size_t ec_idx = it->second;
            ecTable_.ecs[ec_idx].count += other_ec.count;
        } else {
            // New EC: add entry
            EC ec = other_ec;
            size_t ec_idx = ecTable_.ecs.size();
            ecTable_.ecs.push_back(ec);
            ecTable_.n_ecs++;
            ec_signature_map_[sig] = ec_idx;
        }
    }
    
    // Merge GC observations
    observedGC_.combineCounts(other.observedGC_);
}

