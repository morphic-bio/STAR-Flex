# Debug Trace Instrumentation Summary

**Date**: December 19, 2025  
**Status**: ✅ Complete for em_quant  
**Task**: Add detailed per-iteration and per-EC logging to trace divergence between em_quant and Salmon

## What Was Implemented

### Enhanced Debug Tracing in em_quant

Added comprehensive instrumentation to `source/libem/vb_engine.cpp` to log:

1. **Per-iteration values** for selected transcripts:
   - `iter`: Iteration number (-1 for initial state)
   - `transcript`: Transcript name
   - `alpha`: Current alpha value (prior + expected_counts)
   - `logNorm`: digamma(sum of all alphas)
   - `expTheta`: exp(digamma(alpha) - logNorm)
   - `expected_count`: Expected count from E-step

2. **Per-EC details** for ECs containing debug transcripts:
   - `EC`: Marker for EC-level trace
   - `iter`: Iteration number
   - `ec_id`: EC index (0-based)
   - `transcript`: Transcript name
   - `denom`: Denominator (Σ expTheta[tid] * aux[i])
   - `expTheta`: expTheta value for this transcript
   - `aux`: Auxiliary weight (from EC weights or 1/effLen)
   - `contribution`: Per-transcript contribution to expected_counts
   - For single-transcript ECs: Shows "SINGLE" with full count

### Format

**Header**:
```
iter	transcript	alpha	logNorm	expTheta	expected_count
# EC-level trace: EC	iter	ec_id	transcript	denom	expTheta	aux	contribution
```

**Per-iteration line**:
```
0	ENST00000465752	2.88659	8.51603	0.000382601	2.87659
```

**EC-level line**:
```
EC	0	45	ENST00000465752	0.0104374	0.000309963	0.333333	0.00989909
```

**Single-transcript EC**:
```
EC	0	208	ENST00000465752	SINGLE	1	1	1
```

## Usage

### Enable Debug Tracing

1. **Set environment variable** with comma-separated transcript IDs:
   ```bash
   export EM_QUANT_DEBUG_TXPS=ENST00000465752,ENST00000464034
   ```

2. **Run em_quant with debug trace flag**:
   ```bash
   em_quant --vb --threads 1 \
       -e eq_classes.txt \
       -l quant.sf \
       -o output.tsv \
       --debug-trace /tmp/em_trace.txt
   ```

### Example Output

For transcript `ENST00000465752` (one of the outliers):
- Initial state (iter -1): alpha=1.01, expTheta=0.000309963
- Iteration 0: Shows 5 ECs containing this transcript
  - EC 45: denom=0.0104374, contribution=0.00989909
  - EC 71: denom=0.00138222, contribution=1.45763
  - EC 91: denom=0.0155012, contribution=0.299942
  - EC 97: denom=0.0014185, contribution=0.109122
  - EC 208: SINGLE-transcript EC, full count=1
- Final expected_count: 2.87659

## Next Steps: Salmon Instrumentation

To compare traces and find the first divergence, add similar logging to Salmon:

### Location
`/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp` - line 133 (VBEMUpdate_ loop)

### What to Log

For selected transcripts (same IDs as em_quant):
1. **Per-iteration** (after M-step):
   - `iter`: Iteration number
   - `transcript`: Transcript name
   - `alpha`: Current alpha value
   - `logNorm`: digamma(sum of all alphas)
   - `expTheta`: exp(digamma(alpha) - logNorm)
   - `expected_count`: Expected count from E-step

2. **Per-EC** (inside VBEMUpdate_ loop):
   - `EC`: Marker
   - `iter`: Iteration number
   - `ec_id`: EC index
   - `transcript`: Transcript name
   - `denom`: Denominator (Σ expTheta[tid] * aux[i])
   - `expTheta`: expTheta value
   - `aux`: Auxiliary weight (combinedWeights)
   - `contribution`: Per-transcript contribution

### Salmon Command

```bash
salmon quant \
    -i index \
    -l A \
    -1 reads_1.fq -2 reads_2.fq \
    --dumpEqWeights \
    --threads 1 \
    -o salmon_out \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752,ENST00000464034
```

(Note: `--debugTrace` and `--debugTranscripts` flags would need to be added to Salmon)

## Comparing Traces

Once both traces are generated:

1. **Find first divergence**:
   ```bash
   # Compare per-iteration values
   diff <(grep "^[0-9]" em_trace.txt | grep -v "^#") \
        <(grep "^[0-9]" salmon_trace.txt | grep -v "^#")
   
   # Compare EC-level details
   diff <(grep "^EC" em_trace.txt) \
        <(grep "^EC" salmon_trace.txt)
   ```

2. **Check specific iteration/EC**:
   ```bash
   # em_quant iteration 0, EC 45
   grep "^EC.*0.*45" em_trace.txt
   
   # Salmon iteration 0, EC 45
   grep "^EC.*0.*45" salmon_trace.txt
   ```

3. **Identify divergence source**:
   - **Weights parsing/order**: If `aux` values differ
   - **Digamma**: If `logNorm` or `expTheta` differ
   - **Truncation**: If values differ after many iterations
   - **Convergence**: If final values differ

## Files Modified

1. `source/libem/vb_engine.cpp`:
   - Enhanced debug header format
   - Added EC-level logging in E-step loop
   - Updated per-iteration logging format
   - Added single-transcript EC detection and logging

## Testing

✅ Code compiles successfully  
✅ Debug tracing works with environment variable  
✅ Per-iteration values logged correctly  
✅ EC-level details logged for multi-transcript ECs  
✅ Single-transcript ECs logged with "SINGLE" marker  
✅ Trace format matches requirements  
✅ Multi-transcript tracing works correctly  
✅ Helper scripts created and tested

## Helper Scripts Created

1. **`compare_traces.py`**: Compare em_quant and Salmon traces to find first divergence
2. **`extract_trace.py`**: Extract specific iterations, transcripts, or ECs from trace files
3. **`validate_trace.py`**: Validate trace file format and check for inconsistencies
4. **`README_TRACING.md`**: Comprehensive guide on using debug tracing

## Files Created/Modified

1. `source/libem/vb_engine.cpp` - Enhanced debug tracing with EC-level details
2. `tools/em_quant/compare_traces.py` - Trace comparison tool (NEW)
3. `tools/em_quant/extract_trace.py` - Trace extraction tool (NEW)
4. `tools/em_quant/validate_trace.py` - Trace validation tool (NEW)
5. `tools/em_quant/README_TRACING.md` - Tracing documentation (NEW)

## Example Trace Analysis

For transcript `ENST00000465752` (max outlier, 10.11% difference):
- Initial: alpha=1.01 (prior + unique_count=1)
- Iteration 0: expected_count=2.87659
- Final: alpha=1.77278, expected_count=1.76278
- Salmon: NumReads=1.961

The trace shows all ECs contributing to this transcript, allowing comparison with Salmon's trace to identify where divergence occurs.
