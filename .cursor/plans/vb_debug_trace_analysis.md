# VB Debug Trace Analysis and Verification

## Overview

This document summarizes the debug trace analysis for the 9 outlier transcripts that exceed 5% tolerance in VB parity tests.

## Tasks Completed

### 1. Increased Iteration Cap/Tolerance

**Test**: Ran with `--max-iters 1000 --tolerance 1e-10`

**Results**:
- Converged at iteration 121 (vs 26 with default settings)
- **No improvement** in parity - actually slightly worse (10 vs 11 transcripts exceed 5% tolerance)
- **Conclusion**: Not stopping early - differences persist even with tighter convergence

### 2. Per-Iteration Tracing for Outliers

**Implemented**: Debug trace for all 9 outlier transcripts:
- ENST00000430101
- ENST00000464110
- ENST00000406855
- ENST00000472526
- ENST00000492516
- ENST00000403208
- ENST00000482576
- ENST00000642727
- ENST00000476077

**Trace Format**: `iter\ttranscript\talpha\tunique_count\texp_digamma\tweight\texpected_count\tabundance`

**Output File**: `/tmp/debug_all_outliers_final.tsv`

**Usage**:
```bash
EM_QUANT_DEBUG_TXPS=ENST00000430101,ENST00000464110,... ./em_quant --vb --debug-trace debug.tsv ...
```

### 3. Initialization Verification

**Verified**: Alpha seeds are exactly `prior + unique_counts` with no normalization

**Results**:
| Transcript | Unique Counts | Expected Alpha | Actual Alpha | Match |
|------------|---------------|----------------|--------------|-------|
| ENST00000430101 | 12.0 | 12.01 | 12.01 | ✓ |
| ENST00000464110 | 10.0 | 10.01 | 10.01 | ✓ |
| ENST00000406855 | 13.0 | 13.01 | 13.01 | ✓ |
| ENST00000472526 | 10.0 | 10.01 | 10.01 | ✓ |
| ENST00000492516 | 8.0 | 8.01 | 8.01 | ✓ |
| ENST00000403208 | 12.0 | 12.01 | 12.01 | ✓ |
| ENST00000482576 | 12.0 | 12.01 | 12.01 | ✓ |
| ENST00000642727 | 4.0 | 4.01 | 4.01 | ✓ |
| ENST00000476077 | 2.0 | 2.01 | 2.01 | ✓ |

**Conclusion**: ✅ Initialization is correct - `alpha_i = prior (0.01) + unique_counts_i` with no normalization

### 4. Unique Counts Verification

**Verified**: Unique counts come only from single-transcript ECs

**Implementation**: `compute_unique_counts()` correctly filters:
```cpp
if (ec.transcript_ids.size() == 1 && ec.count > 0) {
    unique_counts[ec.transcript_ids[0]] += ec.count;
}
```

**Conclusion**: ✅ Unique counts are computed correctly from single-transcript ECs only

### 5. No Extra Damping/Normalization

**Verified**: 
- ✅ Alpha values are **NOT normalized** - only abundances are normalized
- ✅ M-step: `alpha[i] = prior + expected_counts[i]` (no normalization)
- ✅ No additional smoothing beyond Dirichlet prior
- ✅ No damping factors applied

**Code Verification**:
```cpp
// VB M-step: update alpha = prior + new expected counts
// NOTE: Do NOT renormalize alpha - alpha_i = prior + expected_counts_i only
for (size_t i = 0; i < state.n; ++i) {
    alpha[i] = params.vb_prior + expected_counts[i];
}

// Normalize abundances from alpha (for ELBO computation and next iteration)
// But alpha itself is NOT normalized
double total_alpha = 0.0;
for (size_t i = 0; i < state.n; ++i) {
    total_alpha += alpha[i];
}
if (total_alpha > 0) {
    for (size_t i = 0; i < state.n; ++i) {
        state.abundances[i] = alpha[i] / total_alpha;  // Only abundances normalized
    }
}
```

**Conclusion**: ✅ No extra damping or normalization of alpha

## Trace Analysis

### Example: ENST00000476077 (Largest Outlier)

**Initialization**:
- Alpha: 2.01 (prior 0.01 + unique_counts 2.0) ✓
- exp(digamma(2.01)): 1.53605
- Weight: 0.0153605

**Convergence**:
- Converges to alpha ≈ 7.58
- Final expected count: 7.57 (alpha - prior)
- Salmon reports: 11.49
- Difference: 34.1%

**Observations**:
- VB converges smoothly to a stable solution
- Alpha values are reasonable (not near zero)
- The difference suggests Salmon may use different heuristics for highly ambiguous transcripts

### All Outliers Summary

All 9 outliers share common characteristics:
- **Low unique evidence**: <22% unique reads (most <15%)
- **High ambiguity**: Most reads shared with other transcripts
- **Stable convergence**: VB converges smoothly, not oscillating
- **Reasonable alpha values**: Not near zero or extreme

## Next Steps

### 1. Instrument Salmon (If Possible)

**Challenge**: Adding prints to Salmon's VBEMOptimizer requires:
- Access to Salmon source code
- Recompiling Salmon
- Running on the same fixture

**If successful**: Compare traces to see:
- Where divergence occurs (init vs midway)
- Whether Salmon uses different initialization
- Whether Salmon applies additional heuristics

### 2. Investigate Remaining Heuristics

If traces show divergence despite matching:
- Digamma implementation ✓
- Initialization (prior + unique_counts) ✓
- M-step (alpha = prior + expected_counts) ✓
- No extra normalization ✓

Then investigate:
- **Min alpha clamps**: Does Salmon clamp alpha to a minimum value?
- **Damping factors**: Does Salmon apply damping in early iterations?
- **Different digamma approximation**: Does Salmon use a different digamma implementation?
- **Fragment length distribution**: Even with bias disabled, does Salmon use FLD differently?

### 3. Adopt Salmon Heuristics

Based on trace comparison, adopt any remaining Salmon heuristics:
- Min alpha clamps
- Damping factors
- Different convergence criteria
- Other VB-specific optimizations

## Files

- **Debug trace**: `/tmp/debug_all_outliers_final.tsv`
- **Contains**: Per-iteration alpha, exp(digamma), weight, expected_count for all 9 outliers
- **Format**: TSV with header: `iter\ttranscript\talpha\tunique_count\texp_digamma\tweight\texpected_count\tabundance`

## Conclusions

1. ✅ **Initialization is correct**: `alpha_i = prior + unique_counts_i` with no normalization
2. ✅ **No extra damping/normalization**: Alpha values are not normalized, only abundances
3. ✅ **Unique counts correct**: Only from single-transcript ECs
4. ❌ **Tighter convergence doesn't help**: Differences persist even with 1000 iters and 1e-10 tolerance
5. ✅ **Debug tracing works**: Can trace per-iteration values for selected transcripts

**Remaining differences are likely due to**:
- Salmon's additional heuristics (min clamps, damping, etc.)
- Different digamma approximation (though ours should be accurate)
- Different handling of highly ambiguous transcripts

**Recommendation**: Compare traces with Salmon's VBEMOptimizer to identify remaining differences.
