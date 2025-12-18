# Semiglobal Alignment Implementation - COMPLETE

**Date**: December 18, 2025  
**Status**: ✅ **COMPLETE - Perfect Parity Achieved**  
**Goal**: Match cutadapt's semiglobal alignment algorithm for 3' adapter trimming

## Executive Summary

Successfully ported cutadapt's exact semiglobal alignment algorithm from the Cython source code (`_align.pyx` v5.1). **All 9 tests pass with 0 diff lines** - perfect parity with Trim Galore/cutadapt output.

## Final Results

```
Running synthetic fixture tests...
Testing: adapter_after_quality_trim ... PASS
Testing: below_min_length ... PASS
Testing: clean_long_insert ... PASS
Testing: paired_keep_untrimmed ... PASS
Testing: quality_trim_tail ... PASS
Testing: short_insert_overlap ... PASS
Testing: short_insert_with_errors ... PASS
Synthetic fixture summary: 7 passed, 0 failed

Running integration test (nfcore_smoke)... PASS
Running integration test (nfcore_real)... PASS

Overall summary: 9 passed, 0 failed
  Synthetic fixtures: 7 passed, 0 failed
  Integration tests: 2 passed, 0 failed
```

**Diff counts:** 0 lines in both R1 and R2

## What Was Done

### 1. Downloaded and Analyzed cutadapt Source Code

**Got exact version:**
```bash
python3 -c "import cutadapt; print(cutadapt.__version__)"  # 5.1
pip download cutadapt==5.1 --no-binary :all: -d /tmp/cutadapt-src
tar -xzf cutadapt-5.1.tar.gz
```

**Key source file:** `/tmp/cutadapt-src/cutadapt-5.1/src/cutadapt/_align.pyx`

### 2. Analysis of cutadapt's Algorithm

**Source Code Analysis:**
- Examined `_align.pyx` (Cython source - the actual algorithm)
- Examined `adapters.py` (Python wrapper)

**Key Findings:**

For **3' adapter trimming** (`BackAdapter`), cutadapt uses:
```python
flags = QUERY_START | QUERY_STOP | REFERENCE_END
```

**Semiglobal Alignment Flags:**
- **QUERY_START** (2): Read prefix can be skipped (adapter can start anywhere in read)
- **QUERY_STOP** (8): Read suffix can be skipped (adapter can match partial read)
- **REFERENCE_END** (4): **Adapter suffix can be skipped** (adapter can end early)

**Critical Insight:**
- cutadapt uses **edit distance (Levenshtein)** with indels, not Hamming distance
- Error threshold is calculated as: `floor(aligned_adapter_length * max_error_rate)`
- The **aligned adapter length** (not read length) determines allowed errors
- Example: 13bp read suffix vs 12bp adapter match → `floor(12 * 0.1) = 1` error allowed

**Example from cutadapt:**
```
Read (13bp):    AGATCGGAAGAAG
Adapter (12bp): AGATCGGAAGAG-  ← adapter ends at position 12
                            ^
                         1 insertion (gap in adapter)
Result: 12bp adapter aligned, 1 error, PASS (1 ≤ 1)
```

### 2. Implementation - Exact Port from cutadapt

**File Modified:** `/mnt/pikachu/STAR-Flex/source/libtrim/adapter_trim.cpp`

**Key insight:** The previous attempts failed because they tried to approximate cutadapt's algorithm. The solution was to **port the exact algorithm** from `_align.pyx`.

**Critical details from cutadapt source:**

1. **Entry Structure** (stores more than just cost):
```cpp
struct Entry {
    int cost;    // edit distance
    int score;   // alignment score (matches*1 + mismatches*-1 + indels*-2)
    int origin;  // where alignment started (negative=ref pos, positive=query pos)
};
```

2. **Scoring Constants:**
```cpp
#define MATCH_SCORE      1
#define MISMATCH_SCORE  -1
#define INSERTION_SCORE -2
#define DELETION_SCORE  -2
```

3. **Flags for 3' adapter:**
```cpp
// QUERY_START | QUERY_STOP | REFERENCE_END = 2 | 8 | 4 = 14
int flags = FLAG_START_IN_QUERY | FLAG_STOP_IN_QUERY | FLAG_STOP_IN_REFERENCE;
```

4. **Column Initialization** (depends on flags):
```cpp
// For flags=14 (start_in_query=true, start_in_reference=false):
for (int i = 0; i <= m; i++) {
    column[i].score = i * DELETION_SCORE;  // -2 per deletion
    column[i].cost = i * DELETION_COST;    // 1 per deletion
    column[i].origin = std::max(0, min_n - i);
}
```

5. **DP Fill with origin tracking:**
```cpp
if (characters_equal) {
    cost = diag_entry.cost;
    origin = diag_entry.origin;
    score = diag_entry.score + MATCH_SCORE;  // +1
} else {
    // Choose best of mismatch, deletion, insertion
    // Each tracks origin back to start of alignment
}
```

6. **Match Selection Criteria** (from cutadapt):
```cpp
bool is_acceptable = (
    length >= min_overlap &&
    cost <= (int)(cur_effective_length * max_error_rate)
);

// Update best if:
// - First occurrence, OR
// - Overlaps previous best and has higher score, OR  
// - Longer alignment and higher score
if (is_acceptable && (
    (best.cost == m + n + 1) ||
    (origin <= best.origin + m / 2 && score > best.score) ||
    (length > best_length && score > best.score)
)) {
    // Update best
}
```

7. **Ukkonen's Trick** (optimization):
```cpp
// Only fill rows up to 'last' where cost <= k
while (last >= 0 && column[last].cost > k) {
    last--;
}
```

### 3. Testing Results

**Final Test Status:**
```
Overall summary: 9 passed, 0 failed
  Synthetic fixtures: 7 passed, 0 failed
  Integration tests: 2 passed, 0 failed
```

**Diff Counts:** 0 lines (perfect parity!)

**Test Command:**
```bash
cd /mnt/pikachu/STAR-Flex/source && make test_trim_parity
```

## What Works

✅ **All synthetic fixtures**: 7/7 pass  
✅ **nfcore_smoke integration**: PASS  
✅ **nfcore_real integration**: PASS (0 diff lines!)  
✅ **Exact cutadapt algorithm**: Ported verbatim from `_align.pyx`  
✅ **Entry structure with origin tracking**: Correctly tracks alignment start  
✅ **Scoring system**: MATCH=+1, MISMATCH=-1, INDEL=-2  
✅ **Flag-based initialization**: Correct for QUERY_START | QUERY_STOP | REFERENCE_END  
✅ **Match selection criteria**: Same tie-breaking as cutadapt  
✅ **Ukkonen's trick optimization**: Implemented for performance

## Key Lessons Learned

### Why Previous Attempts Failed

1. **Guessing instead of reading source**: Previous implementations tried to approximate cutadapt's behavior based on documentation and testing. The solution was to read the actual source code.

2. **Missing the `origin` field**: The Entry structure needs to track where the alignment started, not just cost. This is critical for determining the aligned region.

3. **Wrong scoring**: Previous attempts used `adapter_len - 2 * errors`. cutadapt uses `MATCH_SCORE=+1, MISMATCH_SCORE=-1, INDEL_SCORE=-2` accumulated during DP.

4. **Wrong initialization**: Column initialization depends on flags. For QUERY_START, the first row costs increase but origin tracking allows skipping query prefix.

5. **Missing match selection logic**: cutadapt has specific tie-breaking: prefers higher score, then considers overlap with previous best, then considers length.

### The Winning Approach

**Don't approximate - port exactly:**
```bash
# 1. Get exact source for installed version
pip download cutadapt==5.1 --no-binary :all:
tar -xzf cutadapt-5.1.tar.gz

# 2. Read the actual Cython source
cat src/cutadapt/_align.pyx

# 3. Port the algorithm verbatim to C++
```

## Testing Commands

```bash
# Rebuild everything
cd /mnt/pikachu/STAR-Flex/source/libtrim && make clean && make
cd /mnt/pikachu/STAR-Flex/tools/trimvalidate && make clean && make

# Run full test suite
cd /mnt/pikachu/STAR-Flex/source && make test_trim_parity

# Check diff counts
wc -l test/integration/trim/nfcore_real/results/diff_R*.txt

# Compare specific read with cutadapt
python3 -c "
from cutadapt.align import Aligner, EndSkip
adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
read = 'YOUR_READ'
flags = EndSkip.QUERY_START | EndSkip.QUERY_STOP | EndSkip.REFERENCE_END
aligner = Aligner(adapter, max_error_rate=0.1, flags=int(flags), min_overlap=3)
print(aligner.locate(read))
"
```

## Files Modified

1. **`/mnt/pikachu/STAR-Flex/source/libtrim/adapter_trim.cpp`**
   - Complete rewrite with semiglobal alignment
   - Added `semiglobal_align_suffix()` function
   - Modified `find_adapter_3p()` to use semiglobal alignment

2. **Build artifacts:**
   - `source/libtrim/libtrim.a` - rebuilt
   - `tools/trimvalidate/trimvalidate` - rebuilt

## Test Commands

```bash
# Rebuild
cd /mnt/pikachu/STAR-Flex/source/libtrim && make clean && make
cd /mnt/pikachu/STAR-Flex/tools/trimvalidate && make clean && make

# Run tests
cd /mnt/pikachu/STAR-Flex/source && make test_trim_parity

# Check diff counts
wc -l test/integration/trim/nfcore_real/results/diff_R1.txt \
       test/integration/trim/nfcore_real/results/diff_R2.txt

# View specific diffs
head -100 test/integration/trim/nfcore_real/results/diff_R1.txt
```

## Key Code Locations

- **Adapter trimming implementation**: `source/libtrim/adapter_trim.cpp`
- **Test fixtures**: `test/integration/trim/nfcore_real/`
- **Test script**: `tools/trimvalidate/run_parity_test.sh`
- **Debug logs**: `.cursor/debug.log` (NDJSON format)

## Conclusion

**✅ PERFECT PARITY ACHIEVED**

Successfully ported cutadapt's exact semiglobal alignment algorithm from the Cython source code (`_align.pyx` v5.1). The key was to read the actual source rather than guessing based on documentation.

**Final Results:**
- All 9 tests pass
- 0 diff lines in nfcore_real integration test
- Exact match with Trim Galore/cutadapt output

**Files Modified:**
- `/mnt/pikachu/STAR-Flex/source/libtrim/adapter_trim.cpp` - Complete rewrite with exact cutadapt algorithm

**Source Reference:**
- Downloaded from: `pip download cutadapt==5.1 --no-binary :all:`
- Key file: `src/cutadapt/_align.pyx` (Cython source)
- Algorithm: `Aligner.locate()` method (lines 298-587)
