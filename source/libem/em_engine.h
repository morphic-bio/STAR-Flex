#ifndef EM_ENGINE_H
#define EM_ENGINE_H

#include "em_types.h"

// Run EM algorithm on equivalence classes
EMResult run_em(const ECTable& ecs, TranscriptState& state, const EMParams& params);

// Initialize transcript abundances
void initialize_abundances(TranscriptState& state, const EMParams& params);

// Compute log-likelihood of current abundances given ECs (with effective-length weighting)
double compute_log_likelihood(const ECTable& ecs, const double* abundances, const double* eff_lengths);

// Zero out low-abundance transcripts (Salmon-like behavior)
// Sets counts[i] = 0 if counts[i] < threshold, then renormalizes abundances
void zero_low_abundance(std::vector<double>& counts, TranscriptState& state, double threshold);

#endif // EM_ENGINE_H
