# Trace Parity Fixes

**Date**: December 19, 2025  
**Status**: ✅ Initialization Fixed, ⚠️ EC Mismatch Remains

## Issues Identified

### 1. ✅ FIXED: Initialization Mismatch

**Problem**: 
- em_quant: alpha=5.01 (prior + unique_counts)
- Salmon: alpha=71.36 (uniformPrior = totalWeight/numActive)
- Difference: 92.98%

**Solution**: Updated em_quant uniform initialization to match Salmon's approach:
- Compute `totalWeight = sum(EC counts)` (equivalent to Salmon's `sum(projectedCounts)`)
- Compute `uniformPrior = totalWeight / numActive`
- Set `alpha[i] = uniformPrior` (per-transcript) or `uniformPrior * eff_length[i]` (per-nucleotide)

**Result**: 
- em_quant with `--uniform-init`: alpha=71.3429
- Salmon: alpha=71.3637
- Difference: **0.03%** (down from 92.98%!)

### 2. ⚠️ REMAINING: EC Mismatch

**Problem**:
- Fixture EC file: 216 ECs (`test/fixtures/salmon_eq/eq_classes.txt`)
- Salmon generated: 211 ECs (from reads)
- EC IDs don't match → "No common ECs found" in comparison

**Root Cause**: 
- em_quant uses pre-existing fixture EC file
- Salmon generates new ECs from reads during quantification
- Different EC sets = different EC IDs = can't compare per-EC contributions

**Potential Solutions**:

#### Option A: Use Salmon's `--eqclasses` flag (if supported in quant mode)
```bash
salmon quant --eqclasses test/fixtures/salmon_eq/eq_classes.txt \
    --targets test/fixtures/salmon_eq/chr22_trans.fa \
    --debugTrace /tmp/salmon_trace.txt \
    ...
```
**Status**: Need to verify if `--eqclasses` works with quant mode (currently seems to be alignment-mode only)

#### Option B: Regenerate fixture using exact same reads
- Use the exact same synthetic reads that Salmon is using
- Regenerate `eq_classes.txt` fixture
- Ensure EC IDs match

#### Option C: Match ECs by signature instead of ID
- Compare ECs by their transcript membership + weights signature
- More complex but allows comparison even with different EC IDs

## Current Status

### Initialization: ✅ FIXED
- Uniform initialization now matches Salmon (0.03% difference)
- Can compare per-iteration values successfully

### EC Comparison: ⚠️ BLOCKED
- EC IDs don't match (216 vs 211 ECs)
- Need to use same EC file or match by signature

## Next Steps

1. **Test if `--eqclasses` works with quant mode**
   - If yes: Use fixture EC file directly
   - If no: Proceed with Option B or C

2. **If `--eqclasses` doesn't work**:
   - Regenerate fixture using exact same reads Salmon uses
   - Or implement EC matching by signature

3. **Once ECs match**:
   - Compare per-EC contributions to find first divergence
   - Identify root cause (weights parsing, digamma, truncation, convergence)

## Files Modified

- `source/libem/vb_engine.cpp` - Updated uniform initialization to match Salmon's `uniformPrior` approach

## Test Commands

```bash
# Generate em_quant trace with uniform init
export EM_QUANT_DEBUG_TXPS=ENST00000215754,ENST00000215770,ENST00000645799
./tools/em_quant/em_quant --vb --uniform-init --threads 1 \
    -e test/fixtures/salmon_eq/eq_classes.txt \
    -l test/fixtures/salmon_eq/quant.sf \
    -o /tmp/em_output_uniform.tsv \
    --debug-trace /tmp/em_trace_uniform.txt

# Generate Salmon trace (currently uses different ECs)
/mnt/pikachu/salmon/build/src/salmon quant \
    -i /tmp/salmon_chr22_idx_new \
    -l A \
    -1 test/fixtures/salmon_eq/synthetic_reads1.fq \
    -2 test/fixtures/salmon_eq/synthetic_reads2.fq \
    --dumpEqWeights --threads 1 \
    --noLengthCorrection --noFragLengthDist --noEffectiveLengthCorrection \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000215754,ENST00000215770,ENST00000645799 \
    -o /tmp/salmon_trace_out

# Compare traces (per-iteration works, per-EC blocked by EC mismatch)
python3 tools/em_quant/compare_traces.py \
    /tmp/em_trace_uniform.txt /tmp/salmon_trace.txt \
    --transcript ENST00000215754
```
