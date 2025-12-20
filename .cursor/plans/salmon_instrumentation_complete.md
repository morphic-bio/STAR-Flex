# Salmon Instrumentation Implementation Summary

**Date**: December 19, 2025  
**Status**: ✅ Complete  
**Task**: Add debug tracing instrumentation to Salmon for comparison with em_quant

## What Was Implemented

### 1. Added Debug Parameters to SalmonOpts

**File**: `/mnt/pikachu/salmon/include/SalmonOpts.hpp`

Added debug tracing fields:
```cpp
// Debug tracing for comparison with em_quant
bool debugTrace{false};
std::string debugTraceFile;
std::string debugTranscripts; // Comma-separated transcript IDs
```

### 2. Added Command-Line Flags

**File**: `/mnt/pikachu/salmon/src/ProgramOptionsGenerator.cpp`

Added parsing for debug flags:
```cpp
("debugTrace", po::value<std::string>(&(sopt.debugTraceFile)),
 "Enable debug trace output to file for comparison with em_quant. "
 "Logs per-iteration and per-EC values for selected transcripts.")
("debugTranscripts", po::value<std::string>(&(sopt.debugTranscripts)),
 "Comma-separated transcript IDs to trace when --debugTrace is enabled. "
 "Example: ENST00000465752,ENST00000464034")
```

### 3. Modified VBEMUpdate_ Functions

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`

#### Single-Threaded Version (Line 109)
- Added optional debug parameters with default values
- Maintains backward compatibility for bootstrapping

#### Multi-Threaded Version (Line 251)
- Added optional debug parameters
- Added EC-level logging inside parallel loop with mutex protection
- Logs: ec_id, denom, expTheta, aux, contribution
- Handles single-transcript ECs with "SINGLE" marker

### 4. Modified optimize Function

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp` (Line 738)

Added:
- Debug stream setup and header writing
- Parsing of debug transcript IDs
- Initial state logging (iter -1) after initialization
- Per-iteration logging after each VBEMUpdate_ call
- Debug stream closure at end

### 5. Added Required Includes

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`

Added:
```cpp
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <string>
#include <mutex>
```

## Implementation Details

### Thread Safety
- Used `std::mutex` with `std::lock_guard` for thread-safe debug output in parallel loops
- Static mutex ensures single lock across all threads

### Logging Format
Matches em_quant format exactly:
- Per-iteration: `iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count`
- EC-level: `EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution`
- Single-transcript ECs: `EC\titer\tec_id\ttranscript\tSINGLE\t1\t1\tcount`

### Initialization Logging
- Logs initial state (iter -1) after alphas are fully initialized
- Includes uniform/online mixing initialization

## Files Modified

1. `/mnt/pikachu/salmon/include/SalmonOpts.hpp` - Added debug parameters
2. `/mnt/pikachu/salmon/src/ProgramOptionsGenerator.cpp` - Added CLI flags
3. `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp` - Added instrumentation

## Usage

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

## Testing

After compilation, test with:

```bash
# Generate Salmon trace
salmon quant -i index -l A -1 r1.fq -2 r2.fq \
    --dumpEqWeights --threads 1 \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752 \
    -o salmon_out

# Compare with em_quant trace
cd /mnt/pikachu/STAR-Flex/tools/em_quant
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --transcript ENST00000465752
```

## Next Steps

1. **Compile Salmon** with the new instrumentation
2. **Generate Salmon trace** using the same fixture as em_quant
3. **Compare traces** to find first divergence point
4. **Analyze divergence** to identify root cause

## Notes

- Debug tracing only enabled when `--debugTrace` is provided AND `useVBEM` is true
- Single-threaded VBEMUpdate_ (used in bootstrapping) maintains backward compatibility
- All debug parameters are optional with default nullptr/empty values
- Thread-safe logging ensures correct output in parallel execution
