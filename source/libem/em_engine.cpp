#include "em_engine.h"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cstring>

void initialize_abundances(TranscriptState& state, const EMParams& params) {
    if (params.init_by_length) {
        // Weight by transcript length
        double total_length = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_length += state.lengths[i];
        }
        
        if (total_length > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = state.lengths[i] / total_length;
            }
        } else {
            // Fallback to uniform if all lengths are zero
            double uniform = 1.0 / state.n;
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = uniform;
            }
        }
    } else {
        // Uniform initialization
        double uniform = 1.0 / state.n;
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = uniform;
        }
    }
}

// Compute log-likelihood using effective-length-weighted probabilities
// P(read | EC) = sum_i (abundance[i] / eff_length[i]) / sum_j (abundance[j] / eff_length[j])
// This matches Salmon's approach where theta_i / effLen_i gives the probability
double compute_log_likelihood(const ECTable& ecs, const double* abundances, const double* eff_lengths) {
    double ll = 0.0;
    for (const EC& ec : ecs.ecs) {
        double prob = 0.0;
        for (uint32_t tid : ec.transcript_ids) {
            // Weight by abundance / effective_length (Salmon's approach)
            if (eff_lengths[tid] > 0) {
                prob += abundances[tid] / eff_lengths[tid];
            }
        }
        if (prob > 0) {
            ll += ec.count * std::log(prob);
        }
    }
    return ll;
}

void zero_low_abundance(std::vector<double>& counts, TranscriptState& state, double threshold) {
    // Zero out transcripts with counts below threshold (Salmon-like behavior)
    size_t n_zeros = 0;
    double total_counts = 0.0;
    
    for (size_t i = 0; i < state.n; ++i) {
        if (counts[i] < threshold) {
            counts[i] = 0.0;
            state.abundances[i] = 0.0;
            n_zeros++;
        } else {
            total_counts += counts[i];
        }
    }
    
    // Renormalize abundances for non-zero transcripts
    if (total_counts > 0 && n_zeros > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (counts[i] > 0) {
                state.abundances[i] = counts[i] / total_counts;
            }
        }
    }
}

EMResult run_em(const ECTable& ecs, TranscriptState& state, const EMParams& params) {
    EMResult result;
    result.counts.resize(state.n);
    result.tpm.resize(state.n);
    
    // Initialize abundances
    initialize_abundances(state, params);
    
    // Set number of threads
    int num_threads = params.threads;
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
    
    // NOTE: EM runs after alignment is complete, so all threads are available.
    // This avoids nested OpenMP parallelism (alignment uses threads, EM uses threads separately).
    
    // Allocate expected counts array
    std::vector<double> expected_counts(state.n, 0.0);
    
    // Thread-local storage for expected_counts: flat layout [num_threads * n_transcripts] for cache friendliness
    // Each thread writes to its own buffer, then we reduce deterministically
    std::vector<double> expected_counts_tls(num_threads * state.n, 0.0);
    
    // Compute initial log-likelihood (using effective-length weighting)
    double prev_ll = compute_log_likelihood(ecs, state.abundances.data(), state.eff_lengths.data());
    result.final_ll = prev_ll;
    
    // EM iterations
    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        // E-step: compute expected counts using alpha * aux (Salmon's exact approach)
        // Clear thread-local buffers (reuse per iteration to avoid reallocation)
        std::fill(expected_counts_tls.begin(), expected_counts_tls.end(), 0.0);
        std::memset(expected_counts.data(), 0, state.n * sizeof(double));
        
        // Parallel EC loop: each thread writes to its own TLS buffer (no atomics)
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
            int thread_id = omp_get_thread_num();
            double* thread_counts = expected_counts_tls.data() + thread_id * state.n;
            const EC& ec = ecs.ecs[ec_idx];
            size_t groupSize = ec.transcript_ids.size();
            
            // Single-transcript fast path (Salmon behavior: full count)
            if (groupSize == 1) {
                uint32_t tid = ec.transcript_ids[0];
                thread_counts[tid] += ec.count;  // No atomic - each thread has its own buffer
                continue;
            }
            
            // Multi-transcript: compute denominator using alpha * aux
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                // Use alpha * aux (Salmon's EMUpdate approach)
                // If weights available, use them; otherwise fallback to 1/effLen
                double aux = ec.has_weights()
                    ? ec.weights[i]
                    : (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                denom += state.abundances[tid] * aux;
            }
            
            if (denom > 0) {
                double invDenom = ec.count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                    uint32_t tid = ec.transcript_ids[i];
                    double aux = ec.has_weights()
                        ? ec.weights[i]
                        : (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    double contribution = state.abundances[tid] * aux * invDenom;
                    thread_counts[tid] += contribution;  // No atomic - each thread has its own buffer
                }
            }
        }
        
        // Deterministic reduction: sum across thread-local buffers into expected_counts
        // Use fixed thread order (0..num_threads-1) for determinism
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (size_t i = 0; i < state.n; ++i) {
            double sum = 0.0;
            for (int t = 0; t < num_threads; ++t) {
                sum += expected_counts_tls[t * state.n + i];
            }
            expected_counts[i] = sum;
        }
        
        // M-step: update abundances (abundances are proportional to counts)
        double total_counts = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_counts += expected_counts[i];
        }
        
        if (total_counts > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = expected_counts[i] / total_counts;
            }
        }
        
        // Compute log-likelihood (using effective-length weighting)
        double curr_ll = compute_log_likelihood(ecs, state.abundances.data(), state.eff_lengths.data());
        
        // Check convergence
        double rel_change = std::abs(curr_ll - prev_ll) / (std::abs(prev_ll) + 1e-10);
        
        prev_ll = curr_ll;
        result.final_ll = curr_ll;
        result.iterations = iter + 1;
        
        if (rel_change < params.tolerance) {
            result.converged = true;
            break;
        }
    }
    
    // Copy final expected counts (NumReads in Salmon terminology)
    for (size_t i = 0; i < state.n; ++i) {
        result.counts[i] = expected_counts[i];
    }
    
    // Zero out low-abundance transcripts (Salmon-like behavior)
    if (params.zero_threshold > 0) {
        zero_low_abundance(result.counts, state, params.zero_threshold);
    }
    
    // Compute TPM the Salmon way: TPM_i = (count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6
    double total_normalized = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        if (state.eff_lengths[i] > 0) {
            total_normalized += result.counts[i] / state.eff_lengths[i];
        }
    }
    
    if (total_normalized > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (state.eff_lengths[i] > 0) {
                result.tpm[i] = (result.counts[i] / state.eff_lengths[i]) / total_normalized * 1e6;
            } else {
                result.tpm[i] = 0.0;
            }
        }
    }
    
    return result;
}
