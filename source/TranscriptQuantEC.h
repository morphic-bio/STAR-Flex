#ifndef CODE_TranscriptQuantEC
#define CODE_TranscriptQuantEC

#include "IncludeDefine.h"
#include "Transcript.h"
#include "em_types.h"
#include "ec_builder.h"
#include "gc_bias.h"
#include <vector>
#include <unordered_map>
#include <string>

// Wrapper class for building equivalence classes during STAR alignment
// Maintains thread-local EC table and GC observations
class TranscriptQuantEC {
public:
    TranscriptQuantEC(uint num_transcripts);
    
    // Add alignments for a read (called after all alignments for a read are computed)
    // trMult: array of transcript pointers (alignments)
    // nTr: number of alignments
    // auxProbs: alignment probabilities (log-space, normalized)
    void addReadAlignments(Transcript** trMult, uint nTr, const std::vector<double>& auxProbs);
    
    // Simplified method: add alignments using just transcript IDs and weights
    // transcriptIds: vector of transcript IDs this read maps to
    // weights: vector of alignment weights (uniform or from alignment scores)
    void addReadAlignmentsSimple(const std::vector<uint32_t>& transcriptIds, const std::vector<double>& weights);
    
    // Add GC observation from a properly-paired fragment
    // frag_start, frag_end: fragment boundaries in transcript coordinates
    // gc_pct: GC percentage (0-100)
    // weight: alignment probability weight (log-space)
    void addGCObservation(int32_t frag_start, int32_t frag_end, int32_t gc_pct, double weight);
    
    // Finalize EC table (called after all reads processed)
    void finalize();
    
    // Merge EC table from another thread
    void merge(const TranscriptQuantEC& other);
    
    // Get EC table (for quantification)
    const ECTable& getECTable() const { return ecTable_; }
    
    // Get observed GC model (for GC bias correction)
    const GCFragModel& getObservedGC() const { return observedGC_; }
    GCFragModel& getObservedGC() { return observedGC_; }
    
private:
    ECTable ecTable_;
    GCFragModel observedGC_;
    uint num_transcripts_;
    
    // Temporary storage for current read's alignments
    std::vector<RawAlignment> current_read_alignments_;
    
    // Convert STAR Transcript to RawAlignment for EC building
    RawAlignment transcriptToRawAlignment(Transcript* tr, uint32_t transcript_id);
    
    // Hash function for EC signature (sorted transcript IDs)
    struct ECSignature {
        std::vector<uint32_t> transcript_ids;
        std::vector<double> weights;
        
        bool operator==(const ECSignature& other) const {
            if (transcript_ids.size() != other.transcript_ids.size()) return false;
            if (weights.size() != other.weights.size()) return false;
            for (size_t i = 0; i < transcript_ids.size(); ++i) {
                if (transcript_ids[i] != other.transcript_ids[i]) return false;
                if (std::abs(weights[i] - other.weights[i]) > 1e-10) return false;
            }
            return true;
        }
    };
    
    struct ECSignatureHash {
        size_t operator()(const ECSignature& sig) const {
            size_t h = 0;
            for (uint32_t id : sig.transcript_ids) {
                h ^= std::hash<uint32_t>{}(id) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
    
    // Map from EC signature to EC index
    std::unordered_map<ECSignature, size_t, ECSignatureHash> ec_signature_map_;
};

#endif // CODE_TranscriptQuantEC

