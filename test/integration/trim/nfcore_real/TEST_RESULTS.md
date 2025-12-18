# nfcore_real Integration Test Results

**Date**: December 18, 2025  
**Test Dataset**: GSE110004, SRR6357070  
**Source**: `/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_*.fastq.gz`

## Test Setup

### Fixture Creation
- Copied source FASTQs to `test/integration/trim/nfcore_real/`
- Gunzipped to `input_R1.fastq` and `input_R2.fastq`
- Generated expected outputs using Trim Galore 0.6.10 (cutadapt 5.1 backend):
  ```bash
  trim_galore --paired --quality 20 --length 20 \
              --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
              --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
              input_R1.fastq input_R2.fastq
  ```

### Test Configuration
- Quality cutoff: 20
- Minimum length: 20 bp
- Adapter R1: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
- Adapter R2: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
- Error rate: 0.1 (default)
- Min overlap: 1 bp (default)

## Test Results Summary

### Overall Status: **FAIL**

```
Synthetic fixture summary: 7 passed, 0 failed
Integration test (nfcore_smoke): PASS
Integration test (nfcore_real): FAIL

Overall summary: 8 passed, 1 failed
  Synthetic fixtures: 7 passed, 0 failed
  Integration tests: 1 passed, 1 failed
```

### Dataset Statistics

| File | Lines | Read Pairs |
|------|-------|------------|
| `input_R1.fastq` | 200,000 | 50,000 |
| `input_R2.fastq` | 200,000 | 50,000 |
| `expected_R1.fastq` | 198,988 | 49,747 |
| `expected_R2.fastq` | 198,988 | 49,747 |

**Note**: Trim Galore removed 253 read pairs (0.51%) because at least one mate was shorter than 20 bp after trimming.

## Failure Analysis

### Diff File Sizes
- `results/diff_R1.txt`: 16,037 lines
- `results/diff_R2.txt`: 18,802 lines
- Total differences: 34,839 lines

### Observed Differences

#### Pattern 1: Reads Longer in Our Output
Our implementation produces reads that are **longer** than Trim Galore's output in some cases.

**Example Read**: `SRR6357070.41181212`

**Input**:
```
@SRR6357070.41181212 41181212/1 kraken:taxid|4932
TCAAAAGGTTCGGTCGAAAGCAGGGAAGGTGCCACGGCGGATAGGTGGGGGCCTTTTTCAAAGCCTGCTTGATGTCTTCCCAAAGATTGTCCGTGTTTTTC
+
B/<BBF/F//B/<////</<//////////</////7////<<</</F/////<//<F//</</<<///</<///<<<7<7BFF//7/7////7////7//
```

**Expected (Trim Galore)**:
```
TCAAAAGGTTCGGTCGAAAGCAGGGAAGGTGCCACGGCGGATAGGTGGGGGCCTTTTTCAAAGCCTGCTTGATGTCTTCCCAA
```
Length: 85 bp

**Our Output**:
```
TCAAAAGGTTCGGTCGAAAGCAGGGAAGGTGCCACGGCGGATAGGTGGGGGCCTTTTTCAAAGCCTGCTTGATGTCTTCCCAAAGATTGTCCGTGTTTT
```
Length: 101 bp (16 bases longer)

**Quality Scores at End**:
- Last 20 bases: `AAAGATTGTCCGTGTTTTTC`
- Last 20 quality values: `[33, 37, 37, 14, 14, 22, 14, 22, 14, 14, 14, 14, 22, 14, 14, 14, 14, 22, 14, 14]`
- Many quality scores are below Q20 (14-22 range)

**Analysis**: Trim Galore trimmed 16 bases from the 3' end, likely due to:
1. Quality trimming (low quality tail: Q14-Q22)
2. Possible adapter contamination detection

Our implementation did not trim these bases, suggesting:
- Quality trimming algorithm may not be matching Trim Galore's behavior exactly
- Adapter detection may not be finding partial matches that Trim Galore detects

#### Pattern 2: Reads Shorter in Our Output
Some reads are **shorter** in our output than Trim Galore's.

**Example Read**: `SRR6357070.19868080`

**Expected (Trim Galore)**:
```
CTGGGTCGTAAGTGTCAGGATGTTGAATATCACCTCTTGCAACAAATCTAGCTTTATGAGTAGCGTCACGTTTCCTGTTAAAGATAACCATTGAAGATATT
```
Length: 101 bp

**Our Output**:
```
CTGGGTCGTAAGTGTCAGGATGTTGAATATCACCTCTTGCAACAAATCTAGCTTTATGAGTAGCGTCACGTTTCCTGTTAAAGATAACC
```
Length: 87 bp (14 bases shorter)

**Analysis**: Our implementation trimmed more aggressively than Trim Galore, possibly:
- Finding adapter matches that Trim Galore doesn't detect
- Quality trimming more aggressively

## Potential Root Causes

### 1. Quality Trimming Algorithm Differences
- **Issue**: Modified Mott algorithm implementation may differ from cutadapt's implementation
- **Location**: `source/libtrim/quality_trim.cpp`
- **Investigation**: Compare cumulative score calculation and reset logic

### 2. Adapter Matching Differences
- **Issue**: Adapter search algorithm may not match cutadapt's behavior exactly
- **Location**: `source/libtrim/adapter_trim.cpp`
- **Investigation**: 
  - Check if partial adapter matches at read ends are handled correctly
  - Verify error rate calculation: `floor(overlap * 0.1)`
  - Check if adapter search considers quality-trimmed sequence correctly

### 3. Order of Operations
- **Current**: Quality trim → Adapter trim
- **Question**: Does Trim Galore/cutadapt do multiple passes or different ordering?

### 4. Edge Cases in Real Data
- Real-world data may have characteristics not covered by synthetic fixtures:
  - Partial adapter contamination
  - Low-quality regions that don't trigger standard quality trimming
  - Complex quality score patterns

## Next Steps for Investigation

### 1. Compare Specific Reads
```bash
# Extract a failing read pair
cd test/integration/trim/nfcore_real
grep -A 3 "SRR6357070.41181212" input_R1.fastq > /tmp/test_read_R1.fastq
grep -A 3 "SRR6357070.41181212" input_R2.fastq > /tmp/test_read_R2.fastq

# Run Trim Galore on just this read
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            /tmp/test_read_R1.fastq /tmp/test_read_R2.fastq

# Run our trimvalidate
cd ../../../tools/trimvalidate
./trimvalidate -1 /tmp/test_read_R1.fastq -2 /tmp/test_read_R2.fastq \
               -o1 /tmp/our_R1.fastq -o2 /tmp/our_R2.fastq \
               --quality 20 --length 20

# Compare outputs
diff /tmp/test_read_R1_val_1.fq /tmp/our_R1.fastq
```

### 2. Debug Quality Trimming
- Add debug output to `quality_trim_3p()` to see:
  - Cumulative scores at each position
  - Where trimming occurs
  - Compare with cutadapt's behavior on same read

### 3. Debug Adapter Matching
- Add debug output to `find_adapter_3p()` to see:
  - All adapter match attempts
  - Best match found
  - Why matches are accepted/rejected

### 4. Check Cutadapt Documentation
- Review cutadapt's exact algorithm for:
  - Quality trimming (modified Mott)
  - Adapter matching with error rates
  - Order of operations

### 5. Compare with nfcore_smoke Success
- `nfcore_smoke` test passes - what's different?
- Is `nfcore_smoke` synthetic data that doesn't expose these edge cases?
- Do real datasets have characteristics that synthetic fixtures don't cover?

## Files Generated

### Test Results
- `test/integration/trim/nfcore_real/results/status.txt` - Test status (FAIL)
- `test/integration/trim/nfcore_real/results/diff_R1.txt` - Unified diff for R1
- `test/integration/trim/nfcore_real/results/diff_R2.txt` - Unified diff for R2
- `test/integration/trim/nfcore_real/results/README.md` - Results documentation

### Test Fixtures
- `test/integration/trim/nfcore_real/input_R1.fastq` - Input R1 (50,000 reads)
- `test/integration/trim/nfcore_real/input_R2.fastq` - Input R2 (50,000 reads)
- `test/integration/trim/nfcore_real/expected_R1.fastq` - Expected R1 (49,747 reads)
- `test/integration/trim/nfcore_real/expected_R2.fastq` - Expected R2 (49,747 reads)
- `test/integration/trim/nfcore_real/README.md` - Fixture documentation

### Test Script
- `tools/trimvalidate/run_parity_test.sh` - Updated to include nfcore_real test

## Test Command

To rerun the test:
```bash
cd /mnt/pikachu/STAR-Flex/source
make test_trim_parity
```

Or directly:
```bash
cd /mnt/pikachu/STAR-Flex/tools/trimvalidate
./run_parity_test.sh
```

## Notes

- All 7 synthetic fixtures pass, indicating core algorithms work correctly
- `nfcore_smoke` integration test passes, suggesting the issue is specific to real-world data characteristics
- The differences are systematic (not random), suggesting algorithmic differences rather than bugs
- Further investigation needed to understand exact cutadapt behavior on these edge cases
