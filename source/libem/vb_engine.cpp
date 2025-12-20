#include "vb_engine.h"
#include "em_engine.h"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <limits>

// A standalone, high-precision implementation of Digamma (Psi)
// Logic: Uses the recurrence relation psi(x) = psi(x+1) - 1/x to shift 
// small arguments to a stable region (x >= 6), then uses asymptotic expansion.
double robust_digamma(double x) {
    double result = 0.0;

    // 1. Handle the "Dirichlet Cliff" (x near 0) and negative numbers
    // The recurrence relation shifts us to the stable region.
    while (x < 6.0) {
        // If x is extremely close to zero, 1/x will be massive.
        // We rely on double precision to handle this naturally.
        if (x == 0.0) return -std::numeric_limits<double>::infinity(); // Singularity
        result -= 1.0 / x;
        x += 1.0;
    }

    // 2. Stable Region (x >= 6)
    // Use the asymptotic expansion:
    // psi(x) ~ ln(x) - 1/2x - 1/12x^2 + 1/120x^4 - 1/252x^6 ...
    
    // We compute 1/x^2 once to reuse it
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;

    result += std::log(x);
    result -= 0.5 * inv_x;

    // Bernoulli number terms
    double term = inv_x2;       // 1/x^2
    result -= (1.0/12.0) * term; 
    
    term *= inv_x2;             // 1/x^4
    result += (1.0/120.0) * term; 
    
    term *= inv_x2;             // 1/x^6
    result -= (1.0/252.0) * term; 
    
    // Usually sufficient precision for double (term is now ~1e-15 for x=6)
    // Add one more term if strict parity with Boost is required:
    // term *= inv_x2;             // 1/x^8
    // result += (1.0/240.0) * term; 

    return result;
}

// Constants matching Salmon's implementation
constexpr double digammaMin = 1e-10;
constexpr double minEQClassWeight = std::numeric_limits<double>::min();
constexpr double minAlpha = 1e-8;
constexpr double alphaCheckCutoff = 1e-2;

// Compute unique EC counts per transcript (for VB initialization)
std::vector<double> compute_unique_counts(const ECTable& ecs) {
    std::vector<double> unique_counts(ecs.n_transcripts, 0.0);
    
    for (const EC& ec : ecs.ecs) {
        if (ec.transcript_ids.size() == 1 && ec.count > 0) {
            // Single-transcript EC = unique evidence
            unique_counts[ec.transcript_ids[0]] += ec.count;
        }
    }
    
    return unique_counts;
}

// Check if transcript has unique evidence (appears in single-transcript ECs)
std::vector<bool> compute_unique_evidence(const ECTable& ecs) {
    std::vector<bool> has_unique(ecs.n_transcripts, false);
    
    for (const EC& ec : ecs.ecs) {
        if (ec.transcript_ids.size() == 1 && ec.count > 0) {
            has_unique[ec.transcript_ids[0]] = true;
        }
    }
    
    return has_unique;
}

double compute_elbo(const ECTable& ecs, const double* abundances, const double* eff_lengths, double vb_prior) {
    double ll = compute_log_likelihood(ecs, abundances, eff_lengths);
    
    // Add Dirichlet prior term: sum((alpha-1) * log(theta))
    // For uniform Dirichlet: alpha_i = vb_prior for all i
    double prior_term = 0.0;
    for (size_t i = 0; i < ecs.n_transcripts; ++i) {
        if (abundances[i] > 0) {
            prior_term += (vb_prior - 1.0) * std::log(abundances[i]);
        }
    }
    
    return ll + prior_term;
}

EMResult run_vb(const ECTable& ecs, TranscriptState& state, const EMParams& params) {
    EMResult result;
    result.counts.resize(state.n);
    result.tpm.resize(state.n);
    
    // Set number of threads
    int num_threads = params.threads;
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
    
    // Initialize priorAlphas based on mode (per-transcript vs per-nucleotide)
    std::vector<double> priorAlphas(state.n, params.vb_prior);
    if (!params.per_transcript_prior) {
        for (size_t i = 0; i < state.n; ++i) {
            priorAlphas[i] = params.vb_prior * state.eff_lengths[i];
        }
    }
    
    // Initialize alpha_i
    // Option 1 (default, backward compatible): prior + unique_counts_i
    // Option 2 (Salmon's approach): uniform initialization matching Salmon's uniformPrior
    // Salmon uses: uniformPrior = totalWeight / numActive, then alpha = uniformPrior (when fracObserved ≈ 0)
    // When fracObserved ≈ 0: alphas[i] ≈ uniformPrior = totalWeight / numActive
    // This matches Salmon's behavior when projectedCounts are small relative to numRequiredFragments
    std::vector<double> alpha(state.n);
    std::vector<double> unique_counts(state.n, 0.0); // Initialize to zero, compute if needed
    if (params.use_uniform_init) {
        // Uniform initialization matching Salmon's uniformPrior approach
        // Salmon computes: uniformPrior = totalWeight / numActive
        // where totalWeight = sum(txp.projectedCounts) and numActive = number of transcripts
        // When fracObserved ≈ 0 (totalWeight << numRequiredFragments), alphas[i] ≈ uniformPrior
        // 
        // For em_quant: compute totalWeight from EC counts (equivalent to projectedCounts)
        // Then set alpha[i] = uniformPrior = totalWeight / numActive
        double totalWeight = 0.0;
        for (const EC& ec : ecs.ecs) {
            totalWeight += ec.count;
        }
        
        // Compute uniformPrior = totalWeight / numActive (Salmon's approach)
        // If totalWeight is very small, use a small uniform value to match Salmon's behavior
        // when projectedCounts are near zero
        double uniformPrior = totalWeight / static_cast<double>(state.n);
        
        // For per-nucleotide prior: scale uniformPrior by effective length
        // For per-transcript prior: use uniformPrior directly
        if (!params.per_transcript_prior) {
            // Per-nucleotide prior: scale by effective length
            for (size_t i = 0; i < state.n; ++i) {
                alpha[i] = uniformPrior * state.eff_lengths[i];
            }
        } else {
            // Per-transcript prior: uniformPrior is constant per transcript
            for (size_t i = 0; i < state.n; ++i) {
                alpha[i] = uniformPrior;
            }
        }
    } else {
        // Legacy initialization: prior + unique_counts_i
        // IMPORTANT: No normalization of alpha - alpha values are absolute
        unique_counts = compute_unique_counts(ecs);
        for (size_t i = 0; i < state.n; ++i) {
            // Alpha seed: exactly prior + unique_counts (no normalization)
            alpha[i] = priorAlphas[i] + unique_counts[i];
        }
    }
    
    // Initialize abundances from alpha (normalized for E-step, but alpha itself is not normalized)
    double total_alpha = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        total_alpha += alpha[i];
    }
    if (total_alpha > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = alpha[i] / total_alpha;
        }
    } else {
        // Fallback to uniform if all zero
        double uniform = 1.0 / state.n;
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = uniform;
            alpha[i] = params.vb_prior; // Reset alpha to prior
        }
    }
    
    // Allocate expected counts array
    std::vector<double> expected_counts(state.n, 0.0);
    
    // Debug instrumentation: track selected transcripts
    std::vector<size_t> debug_indices;
    std::ofstream debug_stream;
    if (params.debug_trace && !params.debug_file.empty()) {
        debug_stream.open(params.debug_file);
        if (debug_stream.is_open()) {
            debug_stream << "iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count\n";
            debug_stream << "# EC-level trace: EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution\n";
            // Find indices for debug transcripts
            for (const std::string& txp_name : params.debug_transcripts) {
                for (size_t i = 0; i < state.n; ++i) {
                    if (state.names[i] == txp_name) {
                        debug_indices.push_back(i);
                        break;
                    }
                }
            }
            // Log initial values (iter -1) - BEFORE first iteration
            if (!debug_indices.empty()) {
                // Compute logNorm for debug output
                double alphaSum = 0.0;
                for (size_t i = 0; i < state.n; ++i) {
                    alphaSum += alpha[i];
                }
                double logNorm = robust_digamma(alphaSum);
                
                for (size_t idx : debug_indices) {
                    double expTheta_val = (alpha[idx] > digammaMin) ? std::exp(robust_digamma(alpha[idx]) - logNorm) : 0.0;
                    debug_stream << "-1\t" << state.names[idx] << "\t" 
                                << alpha[idx] << "\t" << logNorm << "\t"
                                << expTheta_val << "\t0.0\n";
                    debug_stream.flush(); // Ensure it's written
                }
            }
        }
    }
    
    // Store previous alpha for convergence checking
    std::vector<double> prev_alpha(state.n);
    for (size_t i = 0; i < state.n; ++i) {
        prev_alpha[i] = alpha[i];
    }
    
    // NOTE: VB/EM runs after alignment is complete, so all threads are available.
    // This avoids nested OpenMP parallelism (alignment uses threads, VB uses threads separately).
    
    // Thread-local storage for expected_counts: flat layout [num_threads * n_transcripts] for cache friendliness
    // Each thread writes to its own buffer, then we reduce deterministically
    std::vector<double> expected_counts_tls(num_threads * state.n, 0.0);
    
    // VB iterations
    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        // Compute logNorm = digamma(sum of all alphas) ONCE before E-step (Salmon's approach)
        double alphaSum = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            alphaSum += alpha[i];
        }
        double logNorm = robust_digamma(alphaSum);
        
        // Precompute expTheta vector with digammaMin guard (Salmon's approach)
        std::vector<double> expTheta(state.n, 0.0);
        for (size_t i = 0; i < state.n; ++i) {
            if (alpha[i] > digammaMin) {
                expTheta[i] = std::exp(robust_digamma(alpha[i]) - logNorm);
            } else {
                expTheta[i] = 0.0;  // Zero out very small alphas
            }
        }
        
        // E-step: compute expected counts using expTheta * aux (Salmon's exact approach)
        // Clear thread-local buffers (reuse per iteration to avoid reallocation)
        std::fill(expected_counts_tls.begin(), expected_counts_tls.end(), 0.0);
        std::memset(expected_counts.data(), 0, state.n * sizeof(double));
        
        // Check if we need to log EC details for debug transcripts
        bool log_ec_details = params.debug_trace && !params.debug_file.empty() && !debug_indices.empty();
        
        // Parallel EC loop: each thread writes to its own TLS buffer (no atomics)
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
            int thread_id = omp_get_thread_num();
            double* thread_counts = expected_counts_tls.data() + thread_id * state.n;
            const EC& ec = ecs.ecs[ec_idx];
            size_t groupSize = ec.transcript_ids.size();
            
            // Check if this EC contains any debug transcripts
            bool ec_contains_debug = false;
            if (log_ec_details) {
                for (size_t i = 0; i < groupSize; ++i) {
                    uint32_t tid = ec.transcript_ids[i];
                    for (size_t debug_idx : debug_indices) {
                        if (tid == debug_idx) {
                            ec_contains_debug = true;
                            break;
                        }
                    }
                    if (ec_contains_debug) break;
                }
            }
            
            // Single-transcript fast path (Salmon behavior: full count, no expTheta guard)
            if (groupSize == 1) {
                uint32_t tid = ec.transcript_ids[0];
                thread_counts[tid] += ec.count;  // No atomic - each thread has its own buffer
                
                // Log single-transcript EC if it's a debug transcript
                if (log_ec_details && ec_contains_debug) {
                    #pragma omp critical
                    {
                        debug_stream << "EC\t" << iter << "\t" << ec_idx << "\t" << state.names[tid] 
                                    << "\tSINGLE\t" << ec.count << "\t" << ec.count << "\t" << ec.count << "\n";
                        debug_stream.flush();
                    }
                }
                continue;
            }
            
            // Multi-transcript: compute denominator using expTheta * aux
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                if (expTheta[tid] > 0.0) {
                    // Use expTheta * aux (Salmon's VBEMOptimizer approach)
                    // If weights available, use them; otherwise fallback to 1/effLen
                    double aux = ec.has_weights() 
                        ? ec.weights[i] 
                        : (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    denom += expTheta[tid] * aux;
                }
            }
            
            // Skip EC if denom too small (minEQClassWeight guard)
            if (denom <= minEQClassWeight) {
                continue;
            }
            
            // Distribute counts proportionally
            double invDenom = ec.count / denom;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                if (expTheta[tid] > 0.0) {
                    double aux = ec.has_weights()
                        ? ec.weights[i]
                        : (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    double contribution = expTheta[tid] * aux * invDenom;
                    thread_counts[tid] += contribution;  // No atomic - each thread has its own buffer
                    
                    // Log per-transcript contribution if this is a debug transcript
                    if (log_ec_details && ec_contains_debug) {
                        bool is_debug_txp = false;
                        for (size_t debug_idx : debug_indices) {
                            if (tid == debug_idx) {
                                is_debug_txp = true;
                                break;
                            }
                        }
                        if (is_debug_txp) {
                            #pragma omp critical
                            {
                                debug_stream << "EC\t" << iter << "\t" << ec_idx << "\t" << state.names[tid]
                                            << "\t" << denom << "\t" << expTheta[tid] << "\t" << aux 
                                            << "\t" << contribution << "\n";
                                debug_stream.flush();
                            }
                        }
                    }
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
        
        // VB M-step: update alpha = prior + new expected counts
        // NOTE: Do NOT renormalize alpha - alpha_i = priorAlphas[i] + expected_counts[i] only
        for (size_t i = 0; i < state.n; ++i) {
            alpha[i] = priorAlphas[i] + expected_counts[i];
        }
        
        // Normalize abundances from alpha (for ELBO computation and next iteration)
        double total_alpha = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_alpha += alpha[i];
        }
        
        if (total_alpha > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = alpha[i] / total_alpha;
            }
        }
        
        // Debug instrumentation: log alpha/weights for selected transcripts
        // Log AFTER M-step (alpha updated) but BEFORE convergence check
        if (debug_stream.is_open()) {
            // Recompute logNorm and expTheta for debug output
            double debug_alphaSum = 0.0;
            for (size_t i = 0; i < state.n; ++i) {
                debug_alphaSum += alpha[i];
            }
            double debug_logNorm = robust_digamma(debug_alphaSum);
            
            for (size_t idx : debug_indices) {
                double expTheta_val = (alpha[idx] > digammaMin) ? std::exp(robust_digamma(alpha[idx]) - debug_logNorm) : 0.0;
                debug_stream << iter << "\t" << state.names[idx] << "\t"
                            << alpha[idx] << "\t" << debug_logNorm << "\t"
                            << expTheta_val << "\t" << expected_counts[idx] << "\n";
                debug_stream.flush();
            }
        }
        
        // Check convergence using alpha relative difference (Salmon's approach)
        bool converged = true;
        double maxRelDiff = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            if (alpha[i] > alphaCheckCutoff) {
                double relDiff = std::abs(alpha[i] - prev_alpha[i]) / alpha[i];
                maxRelDiff = std::max(maxRelDiff, relDiff);
                if (relDiff > params.tolerance) {
                    converged = false;
                }
            }
            prev_alpha[i] = alpha[i];
        }
        
        // Compute ELBO for final result (but don't use for convergence)
        double curr_elbo = compute_elbo(ecs, state.abundances.data(), state.eff_lengths.data(), params.vb_prior);
        result.final_ll = curr_elbo;
        result.iterations = iter + 1;
        
        if (converged) {
            result.converged = true;
            break;
        }
    }
    
    if (debug_stream.is_open()) {
        debug_stream.close();
    }
    
    // Copy final expected counts: use alpha - priorAlphas (after final M-step)
    // In VB, alpha[i] = priorAlphas[i] + expected_counts[i], so counts = alpha[i] - priorAlphas[i]
    for (size_t i = 0; i < state.n; ++i) {
        result.counts[i] = std::max(0.0, alpha[i] - priorAlphas[i]);
    }
    
    // Post-convergence truncation: zero out counts below minAlpha (Salmon's approach)
    for (size_t i = 0; i < state.n; ++i) {
        if (result.counts[i] <= minAlpha) {
            result.counts[i] = 0.0;
        }
    }
    
    // Zero out transcripts where alpha_i <= priorAlphas[i] + epsilon AND no unique support
    // This matches Salmon's implicit zeroing behavior
    std::vector<bool> has_unique = compute_unique_evidence(ecs);
    const double epsilon = 1e-8;
    double total_counts = 0.0;
    
    for (size_t i = 0; i < state.n; ++i) {
        if (alpha[i] <= priorAlphas[i] + epsilon && !has_unique[i]) {
            // Zero out: alpha very close to prior AND no unique evidence
            result.counts[i] = 0.0;
            state.abundances[i] = 0.0;
        } else {
            total_counts += result.counts[i];
        }
    }
    
    // Renormalize abundances for non-zero transcripts
    if (total_counts > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (result.counts[i] > 0) {
                state.abundances[i] = result.counts[i] / total_counts;
            }
        }
    }
    
    // Compute TPM the Salmon way: TPM_i = (count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6
    double total_normalized = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        if (state.eff_lengths[i] > 0 && result.counts[i] > 0) {
            total_normalized += result.counts[i] / state.eff_lengths[i];
        }
    }
    
    if (total_normalized > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (state.eff_lengths[i] > 0 && result.counts[i] > 0) {
                result.tpm[i] = (result.counts[i] / state.eff_lengths[i]) / total_normalized * 1e6;
            } else {
                result.tpm[i] = 0.0;
            }
        }
    }
    
    return result;
}
