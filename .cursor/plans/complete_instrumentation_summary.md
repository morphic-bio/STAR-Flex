# Complete Debug Instrumentation Summary

**Date**: December 19, 2025  
**Status**: ✅ Complete for both em_quant and Salmon  
**Task**: Add debug tracing to both em_quant and Salmon to find divergence points

## ✅ em_quant Instrumentation (Complete)

### Implementation
- Enhanced debug tracing in `source/libem/vb_engine.cpp`
- Per-iteration logging: alpha, logNorm, expTheta, expected_count
- Per-EC logging: ec_id, denom, expTheta, aux, contribution
- Single-transcript EC detection and logging
- Multi-transcript support

### Helper Tools
- `compare_traces.py` - Compare em_quant and Salmon traces
- `extract_trace.py` - Extract specific iterations/transcripts/ECs
- `validate_trace.py` - Validate trace file format
- `README_TRACING.md` - Comprehensive documentation

## ✅ Salmon Instrumentation (Complete)

### Files Modified

1. **`/mnt/pikachu/salmon/include/SalmonOpts.hpp`**
   - Added: `bool debugTrace{false}`
   - Added: `std::string debugTraceFile`
   - Added: `std::string debugTranscripts`

2. **`/mnt/pikachu/salmon/src/ProgramOptionsGenerator.cpp`**
   - Added: `--debugTrace <file>` flag
   - Added: `--debugTranscripts <comma-separated IDs>` flag

3. **`/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`**
   - Modified single-threaded `VBEMUpdate_` (line 109): Added optional debug parameters
   - Modified multi-threaded `VBEMUpdate_` (line 251): Added EC-level logging with mutex
   - Modified `optimize` function (line 738): Added debug setup, initial state logging, per-iteration logging
   - Added includes: `<unordered_set>`, `<fstream>`, `<sstream>`, `<string>`, `<mutex>`

### Implementation Details

**Thread Safety**: Uses `std::mutex` with `std::lock_guard` for thread-safe debug output in parallel loops

**Logging Format**: Matches em_quant format exactly:
- Per-iteration: `iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count`
- EC-level: `EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution`
- Single-transcript ECs: `EC\titer\tec_id\ttranscript\tSINGLE\t1\t1\tcount`

**Initialization**: Logs initial state (iter -1) after alphas are fully initialized

## Usage

### Generate em_quant Trace

```bash
export EM_QUANT_DEBUG_TXPS=ENST00000465752,ENST00000464034
em_quant --vb --threads 1 \
    -e eq_classes.txt -l quant.sf \
    -o em_output.tsv \
    --debug-trace /tmp/em_trace.txt
```

### Generate Salmon Trace

```bash
salmon quant \
    -i /tmp/salmon_chr22_idx \
    -l A \
    -1 /tmp/sub_R1.fq -2 /tmp/sub_R2.fq \
    --dumpEqWeights \
    --threads 1 \
    --noLengthCorrection \
    --noFragLengthDist \
    --noEffectiveLengthCorrection \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752,ENST00000464034 \
    -o /tmp/salmon_out
```

### Compare Traces

```bash
cd tools/em_quant
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --transcript ENST00000465752
```

This will identify:
- First divergence in per-iteration values
- First divergence in EC-level details
- Field-by-field differences

## Next Steps

1. **Compile Salmon** with the new instrumentation
2. **Generate Salmon trace** using the same fixture as em_quant
3. **Compare traces** to find first divergence point
4. **Analyze divergence** to identify root cause:
   - **Weights parsing/order**: If `aux` values differ
   - **Digamma**: If `logNorm` or `expTheta` differ
   - **Truncation**: If values differ after many iterations
   - **Convergence**: If final values differ

## Files Created/Modified

### em_quant
- `source/libem/vb_engine.cpp` - Enhanced debug tracing
- `tools/em_quant/compare_traces.py` - Trace comparison tool
- `tools/em_quant/extract_trace.py` - Trace extraction tool
- `tools/em_quant/validate_trace.py` - Trace validation tool
- `tools/em_quant/README_TRACING.md` - Documentation

### Salmon
- `include/SalmonOpts.hpp` - Added debug parameters
- `src/ProgramOptionsGenerator.cpp` - Added CLI flags
- `src/CollapsedEMOptimizer.cpp` - Added instrumentation

## Summary

Both em_quant and Salmon now have complete debug tracing instrumentation. Once Salmon is compiled, traces can be generated and compared to identify exactly where divergence occurs, enabling precise debugging of the remaining differences.
