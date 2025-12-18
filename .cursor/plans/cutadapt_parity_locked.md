# Cutadapt Parity - Locked State

**Date**: December 18, 2025  
**Status**: ✅ **LOCKED - Perfect Parity Achieved**  
**cutadapt Version**: **5.1** (explicitly synchronized)

## Summary

The adapter trimming implementation in `source/libtrim/adapter_trim.cpp` achieves **perfect parity** with cutadapt v5.1 and Trim Galore output:
- **All 9 tests pass**
- **0 diff lines** in integration tests
- **Exact algorithm port** from cutadapt `_align.pyx` v5.1

## Version Lock

**cutadapt Version**: 5.1  
**Source Reference**: `src/cutadapt/_align.pyx`, `Aligner.locate()` method (lines 298-587)  
**Download Command**: `pip download cutadapt==5.1 --no-binary :all:`

## Code Location

- **Implementation**: `source/libtrim/adapter_trim.cpp`
- **Documentation**: `docs/trimming.md`
- **Test Suite**: `tools/trimvalidate/run_parity_test.sh`
- **CI Integration**: `.travis.yml` (runs `make test_trim_parity`)

## Key Algorithm Details

- **Entry Structure**: `{cost, score, origin}` - tracks edit distance, alignment score, and alignment start position
- **Scoring**: `MATCH=+1, MISMATCH=-1, INSERTION=-2, DELETION=-2`
- **Flags**: `QUERY_START | QUERY_STOP | REFERENCE_END` (flags = 14)
- **Error Threshold**: `floor(aligned_adapter_length × max_error_rate)`
- **Match Selection**: Prefers higher score, then considers overlap with previous best, then length

## Regression Prevention

1. **CI Integration**: `.travis.yml` runs `make test_trim_parity` on every build
2. **Test Suite**: `tools/trimvalidate/run_parity_test.sh` validates all fixtures
3. **Version Documentation**: Code comments reference cutadapt v5.1 explicitly

## Upgrade Procedure

If cutadapt is upgraded in the future:

1. **Download new source**:
   ```bash
   pip download cutadapt==<NEW_VERSION> --no-binary :all: -d /tmp/cutadapt-src
   tar -xzf cutadapt-<NEW_VERSION>.tar.gz
   ```

2. **Compare algorithms**:
   - Read `src/cutadapt/_align.pyx` in new version
   - Compare with current implementation in `source/libtrim/adapter_trim.cpp`
   - Note any changes to scoring, flags, or match selection logic

3. **Update code if needed**:
   - Port any algorithm changes
   - Update version references in code comments
   - Update this document

4. **Regenerate fixtures** (if algorithm changed):
   ```bash
   # For each fixture in test/fixtures/trim/:
   cd test/fixtures/trim/<fixture>/
   trim_galore --paired --quality 20 --length 20 \
               --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
               --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
               input_R1.fastq input_R2.fastq
   mv input_R1_val_1.fq expected_R1.fastq
   mv input_R2_val_2.fq expected_R2.fastq
   ```

5. **Re-run parity tests**:
   ```bash
   cd source && make test_trim_parity
   ```

6. **Verify zero diffs**:
   ```bash
   wc -l test/integration/trim/nfcore_real/results/diff_R*.txt
   # Should show 0 lines
   ```

## Test Results (Current State)

```
Overall summary: 9 passed, 0 failed
  Synthetic fixtures: 7 passed, 0 failed
  Integration tests: 2 passed, 0 failed

Diff counts: 0 lines (perfect parity!)
```

## Files Modified for Parity Lock

- ✅ `source/libtrim/adapter_trim.cpp` - Added detailed source documentation header
- ✅ `docs/trimming.md` - Updated with cutadapt v5.1 reference
- ✅ `.travis.yml` - Added parity test to CI
- ✅ `tools/trimvalidate/run_parity_test.sh` - Added version warning
- ✅ `CHANGES.md` - Added entry documenting parity achievement
- ✅ `.cursor/plans/cutadapt_parity_locked.md` - This document

## Performance Notes

The current implementation uses:
- Space-efficient DP: single column (O(m) space) instead of full matrix (O(m×n))
- Ukkonen's trick: only fills rows where cost <= k
- Early termination: stops on exact match (cost=0, origin>=0)

Performance has not been formally benchmarked but should be comparable to cutadapt's Cython implementation.
