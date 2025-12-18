#include "trim.h"
#include <cstring>
#include <algorithm>
#include <climits>
#include <vector>
#include <cmath>

/*
 * Cutadapt Semiglobal Alignment Algorithm
 * ========================================
 * 
 * This file contains a verbatim port of cutadapt's semiglobal alignment algorithm
 * from cutadapt v5.1 (_align.pyx, Aligner.locate() method).
 * 
 * Source Reference:
 *   - cutadapt version: 5.1
 *   - Source file: src/cutadapt/_align.pyx
 *   - Method: Aligner.locate() (lines 298-587)
 *   - Downloaded via: pip download cutadapt==5.1 --no-binary :all:
 * 
 * Algorithm Details:
 *   - Uses dynamic programming with edit distance (Levenshtein)
 *   - Tracks cost (errors), score (matches - mismatches - indels), and origin (alignment start)
 *   - Supports semiglobal alignment via flags: QUERY_START, QUERY_STOP, REFERENCE_END
 *   - For 3' adapters: uses QUERY_START | QUERY_STOP | REFERENCE_END (flags = 14)
 * 
 * Scoring Constants (from cutadapt _align.pyx):
 *   - MATCH_SCORE = +1
 *   - MISMATCH_SCORE = -1
 *   - INSERTION_SCORE = -2
 *   - DELETION_SCORE = -2
 * 
 * Error Threshold:
 *   - errors <= floor(aligned_adapter_length * max_error_rate)
 *   - aligned_adapter_length is the adapter portion actually matched (ref_stop - ref_start)
 * 
 * IMPORTANT: If cutadapt is upgraded, this code must be re-validated against the new version.
 * Run: make test_trim_parity to verify parity with Trim Galore/cutadapt output.
 */

// Flags (matching cutadapt's EndSkip enum from align.py):
#define FLAG_START_IN_REFERENCE 1
#define FLAG_START_IN_QUERY     2
#define FLAG_STOP_IN_REFERENCE  4
#define FLAG_STOP_IN_QUERY      8

// Scoring constants (from cutadapt _align.pyx)
#define MATCH_SCORE      1
#define MISMATCH_SCORE  -1
#define INSERTION_SCORE -2
#define DELETION_SCORE  -2
#define INSERTION_COST   1
#define DELETION_COST    1

// DP matrix entry (from cutadapt _align.pyx)
struct Entry {
    int cost;    // edit distance
    int score;   // alignment score
    int origin;  // where alignment started: negative = ref position, positive = query position
};

// Match result
struct Match {
    int origin;
    int cost;
    int score;
    int ref_stop;
    int query_stop;
};

// Semiglobal alignment - exact port of cutadapt's Aligner.locate()
// Returns: (ref_start, ref_stop, query_start, query_stop, score, errors) or all -1 if no match
static void semiglobal_locate(
    const char* reference, int m,  // adapter
    const char* query, int n,       // read
    int flags,
    double max_error_rate,
    int min_overlap,
    // Output parameters
    int* out_ref_start, int* out_ref_stop,
    int* out_query_start, int* out_query_stop,
    int* out_score, int* out_errors
) {
    // Initialize output to "no match"
    *out_ref_start = *out_ref_stop = *out_query_start = *out_query_stop = -1;
    *out_score = 0;
    *out_errors = -1;
    
    if (m == 0 || n == 0) return;
    
    // Parse flags
    bool start_in_reference = (flags & FLAG_START_IN_REFERENCE) != 0;
    bool start_in_query = (flags & FLAG_START_IN_QUERY) != 0;
    bool stop_in_reference = (flags & FLAG_STOP_IN_REFERENCE) != 0;
    bool stop_in_query = (flags & FLAG_STOP_IN_QUERY) != 0;
    
    // Maximum number of errors
    int k = (int)(max_error_rate * m);
    
    // Determine largest and smallest column we need to compute
    int max_n = n;
    int min_n = 0;
    if (!start_in_query) {
        max_n = std::min(n, m + k);
    }
    if (!stop_in_query) {
        min_n = std::max(0, n - m - k);
    }
    
    // Allocate column
    std::vector<Entry> column(m + 1);
    
    // Fill column min_n (first column initialization)
    // Four cases based on start_in_reference and start_in_query
    if (!start_in_reference && !start_in_query) {
        for (int i = 0; i <= m; i++) {
            column[i].score = i * DELETION_SCORE;
            column[i].cost = std::max(i, min_n) * DELETION_COST;
            column[i].origin = 0;
        }
    } else if (start_in_reference && !start_in_query) {
        for (int i = 0; i <= m; i++) {
            column[i].score = 0;
            column[i].cost = min_n * DELETION_COST;
            column[i].origin = std::min(0, min_n - i);
        }
    } else if (!start_in_reference && start_in_query) {
        for (int i = 0; i <= m; i++) {
            column[i].score = i * DELETION_SCORE;
            column[i].cost = i * DELETION_COST;
            column[i].origin = std::max(0, min_n - i);
        }
    } else {  // start_in_reference && start_in_query
        for (int i = 0; i <= m; i++) {
            column[i].score = 0;
            column[i].cost = std::min(i, min_n) * DELETION_COST;
            column[i].origin = min_n - i;
        }
    }
    
    // Initialize best match
    Match best;
    best.ref_stop = m;
    best.query_stop = n;
    best.cost = m + n + 1;  // sentinel value indicating "no match yet"
    best.origin = 0;
    best.score = 0;
    
    // Ukkonen's trick: index of the last cell that is at most k
    int last = std::min(m, k + 1);
    if (start_in_reference) {
        last = m;
    }
    
    int origin_increment = start_in_query ? 1 : 0;
    int insertion_cost_increment = start_in_query ? 0 : INSERTION_COST;
    int insertion_score_increment = start_in_query ? 0 : INSERTION_SCORE;
    
    int last_filled_i = 0;
    
    // Iterate over columns
    for (int j = min_n + 1; j <= max_n; j++) {
        // Remember first entry before overwriting
        Entry diag_entry = column[0];
        
        // Fill in first entry in this column
        column[0].origin += origin_increment;
        column[0].cost += insertion_cost_increment;
        column[0].score += insertion_score_increment;
        
        for (int i = 1; i <= last; i++) {
            bool characters_equal = (reference[i-1] == query[j-1]);
            
            int cost, origin, score;
            
            if (characters_equal) {
                // Match - can skip computing costs for insertion and deletion
                cost = diag_entry.cost;
                origin = diag_entry.origin;
                score = diag_entry.score + MATCH_SCORE;
            } else {
                // Characters do not match
                Entry current_entry = column[i];
                Entry previous_entry = column[i-1];
                
                int cost_diag = diag_entry.cost + 1;
                int cost_insertion = current_entry.cost + INSERTION_COST;
                int cost_deletion = previous_entry.cost + DELETION_COST;
                
                if (cost_diag <= cost_deletion && cost_diag <= cost_insertion) {
                    // MISMATCH
                    cost = cost_diag;
                    origin = diag_entry.origin;
                    score = diag_entry.score + MISMATCH_SCORE;
                } else if (cost_deletion <= cost_insertion) {
                    // DELETION
                    cost = cost_deletion;
                    origin = previous_entry.origin;
                    score = previous_entry.score + DELETION_SCORE;
                } else {
                    // INSERTION
                    cost = cost_insertion;
                    origin = current_entry.origin;
                    score = current_entry.score + INSERTION_SCORE;
                }
            }
            
            // Remember current cell for next iteration
            diag_entry = column[i];
            
            column[i].cost = cost;
            column[i].origin = origin;
            column[i].score = score;
        }
        
        last_filled_i = last;
        
        // Ukkonen's trick: shrink last
        while (last >= 0 && column[last].cost > k) {
            last--;
        }
        if (last < m) {
            last++;
        } else if (stop_in_query) {
            // Found a match. If requested, find best match in last row.
            int cost = column[m].cost;
            int score = column[m].score;
            int origin = column[m].origin;
            int length = m + std::min(origin, 0);
            int cur_effective_length = length;
            
            bool is_acceptable = (
                length >= min_overlap &&
                cost <= (int)(cur_effective_length * max_error_rate)
            );
            
            int best_length = m + std::min(best.origin, 0);
            
            // Update best if:
            // - this is the first occurrence
            // - or this occurrence overlaps the previous best one and has a higher score
            // - or if this occurrence has greater length and higher score
            if (is_acceptable && (
                (best.cost == m + n + 1) ||  // No best match recorded so far
                (origin <= best.origin + m / 2 && score > best.score) ||  // Overlaps and better score
                (length > best_length && score > best.score)  // Longer and better score
            )) {
                best.score = score;
                best.cost = cost;
                best.origin = origin;
                best.ref_stop = m;
                best.query_stop = j;
                
                // Exact match, stop early
                if (cost == 0 && origin >= 0) {
                    break;
                }
            }
        }
    }
    
    // Search in last column (if stop_in_reference is set)
    if (max_n == n) {
        int first_i = stop_in_reference ? 0 : m;
        
        for (int i = last_filled_i; i >= first_i; i--) {
            int length = i + std::min(column[i].origin, 0);
            int cost = column[i].cost;
            int score = column[i].score;
            int cur_effective_length = length;
            
            bool is_acceptable = (
                length >= min_overlap &&
                cost <= (int)(cur_effective_length * max_error_rate)
            );
            
            int best_length = best.ref_stop + std::min(best.origin, 0);
            
            if (is_acceptable && (
                (best.cost == m + n + 1) ||
                (column[i].origin <= best.origin + m / 2 && score > best.score) ||
                (length > best_length && score > best.score)
            )) {
                best.score = score;
                best.cost = cost;
                best.origin = column[i].origin;
                best.ref_stop = i;
                best.query_stop = n;
            }
        }
    }
    
    // Check if we found a valid match
    if (best.cost == m + n + 1) {
        return;  // No valid match found
    }
    
    // Compute final positions
    int ref_start, query_start;
    if (best.origin >= 0) {
        ref_start = 0;
        query_start = best.origin;
    } else {
        ref_start = -best.origin;
        query_start = 0;
    }
    
    *out_ref_start = ref_start;
    *out_ref_stop = best.ref_stop;
    *out_query_start = query_start;
    *out_query_stop = best.query_stop;
    *out_score = best.score;
    *out_errors = best.cost;
}

// Find best 3' adapter match position using cutadapt's semiglobal alignment
// Returns position where read should be trimmed (seq_len if no match)
uint32_t find_adapter_3p(const char* seq, uint32_t seq_len,
                         const char* adapter, uint32_t adapter_len,
                         uint32_t min_overlap, double max_error_rate) {
    if (seq_len == 0 || adapter_len == 0) return seq_len;
    
    // Flags for 3' adapter (BackAdapter): QUERY_START | QUERY_STOP | REFERENCE_END
    int flags = FLAG_START_IN_QUERY | FLAG_STOP_IN_QUERY | FLAG_STOP_IN_REFERENCE;
    
    int ref_start, ref_stop, query_start, query_stop, score, errors;
    
    semiglobal_locate(
        adapter, (int)adapter_len,
        seq, (int)seq_len,
        flags,
        max_error_rate,
        (int)min_overlap,
        &ref_start, &ref_stop, &query_start, &query_stop, &score, &errors
    );
    
    if (errors < 0) {
        return seq_len;  // No match
    }
    
    // For 3' adapter, trim at query_start (where adapter match begins in read)
    return (uint32_t)query_start;
}
