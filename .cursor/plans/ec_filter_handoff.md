# EC Filter Module Handoff Document

## Overview

This document describes the EC (Equivalence Class) Gating and Cleanup Module implementation, which replicates Salmon's alignment filtering and EC construction pipeline with optional extended pruning modes.

## Status

**Implementation Status**: Core library complete, CLI skeleton created, test harness ready

**Compilation**: ✅ Library compiles successfully (`libem.a` builds)

**Testing**: ⚠️ Unit tests and parity tests require CLI integration (BAM reading)

## What Has Been Implemented

### Core Library (`source/libem/`)

#### 1. Alignment Filtering (`alignment_filter.h/cpp`)

**Purpose**: Filter alignments using Salmon's exact semantics before EC construction.

**Key Components**:
- `RawAlignment` - Structure matching Salmon's QuasiAlignment
- `MappingScoreInfo` - Tracks best scores per read and per transcript
- `FilterParams` - Configuration matching Salmon defaults
- `updateRefMappings()` - Per-transcript best hit tracking (keeps only best alignment per transcript)
- `filterAndCollectAlignments()` - Applies hard/soft filtering with minAlnProb gating

**Salmon Compatibility**:
- ✅ Per-transcript best hit only (matches `updateRefMappings` logic)
- ✅ Hard filter vs soft filter modes
- ✅ minAlnProb gating: `estAlnProb = exp(-scoreExp * (bestScore - currScore))`
- ✅ Decoy threshold handling
- ✅ maxReadOccs filtering

**Key Functions**:
```cpp
void updateRefMappings(uint32_t tid, int32_t score, bool is_compat, size_t idx,
                      const std::vector<bool>& is_decoy,
                      MappingScoreInfo& msi, std::vector<int32_t>& scores);

std::vector<RawAlignment> filterAndCollectAlignments(
    std::vector<RawAlignment>& alignments,
    const FilterParams& params,
    MappingScoreInfo& msi);
```

#### 2. EC Builder (`ec_builder.h/cpp`)

**Purpose**: Construct equivalence classes from filtered alignments using Salmon's auxProb computation.

**Key Components**:
- `ECBuilderParams` - Configuration (range factorization, rank ECs, etc.)
- `ReadMapping` - Intermediate representation with log-space auxProbs
- `logAdd()` - Numerically stable log-sum-exp: `log(exp(a) + exp(b))`
- `computeAuxProbs()` - Computes `auxProb = logFragProb + logFragCov + logAlignCompatProb`
- `applyRangeFactorization()` - Appends bin IDs: `rangeCount = sqrt(n) + bins`
- `applyRankEqClasses()` - Sorts transcripts by conditional probability
- `buildEquivalenceClasses()` - Main entry point for EC construction

**Salmon Compatibility**:
- ✅ Log-sum-exp normalization matching Salmon's `logAdd`
- ✅ Range factorization formula: `rangeNumber = auxProb * (sqrt(n) + bins)`
- ✅ Rank-based ECs (sort by conditional probability)
- ✅ Proper weight normalization per read

**Key Functions**:
```cpp
ReadMapping computeAuxProbs(const std::vector<RawAlignment>& alignments,
                             const ECBuilderParams& params);

void applyRangeFactorization(std::vector<uint32_t>& txp_ids,
                             const std::vector<double>& aux_probs,
                             uint32_t range_factorization_bins);

ECTable buildEquivalenceClasses(
    const std::vector<std::vector<RawAlignment>>& read_alignments,
    const ECBuilderParams& params,
    size_t num_transcripts);
```

#### 3. Extended Pruning (`extended_pruning.h/cpp`)

**Purpose**: Optional non-Salmon pruning modes (disabled by default for parity).

**Key Components**:
- `ExtendedPruningParams` - Configuration (disabled by default)
- `applyLocalPruning()` - Remove low-probability transcripts from individual ECs
- `applyGlobalPruning()` - Remove transcripts with low total weight across all ECs
- `applyExtendedPruning()` - Combined pruning function

**Important**: These are **NOT** Salmon behaviors. Must be disabled (`enable_local_pruning = false`, `enable_global_pruning = false`) when testing for Salmon parity.

**Key Functions**:
```cpp
void applyLocalPruning(EC& ec, double min_local_prob);
void applyGlobalPruning(ECTable& ecs, double min_global_abundance);
void applyExtendedPruning(ECTable& ecs, const ExtendedPruningParams& params);
```

### Test Harness (`tools/ec_filter_test/`)

#### Unit Tests

1. **`generate_fixtures.py`** - Generates test scenarios:
   - hardFilter behavior
   - minAlnProb gating
   - Per-transcript best hit
   - Decoy threshold
   - Range factorization bins
   - rankEqClasses sorting

2. **`verify_filtering.py`** - Verifies test results with assertions

#### Parity Tests

1. **`generate_parity_bam.sh`** - Creates BAM file from synthetic reads using minimap2

2. **`run_salmon_parity.sh`** - Runs Salmon and our implementation with matching parameters:
   - Uses `--dumpEqWeights` for weight output
   - Disables bias corrections (`--noLengthCorrection`, etc.)
   - Matches scoreExp, minAlnProb, rangeFactorizationBins
   - **Critical**: Uses `--no-local-pruning --no-global-pruning` flags

3. **`compare_ecs.py`** - Python harness that:
   - Parses Salmon EC format (handles gzipped files)
   - Canonicalizes ECs by (labels, weights, count)
   - Reports detailed mismatches
   - Calculates parity percentage (PASS if >=99%)

#### Build System

- **`Makefile`** - Builds `ec_filter_cli` and provides `test` and `parity-test` targets

## What Still Needs to Be Done

### 1. CLI Integration (High Priority)

**File**: `tools/ec_filter_test/ec_filter_cli.cpp`

**Current Status**: Skeleton with argument parsing only

**Required Work**:
- Integrate BAM file reading (htslib/samtools)
- Convert BAM alignments to `RawAlignment` format
- Implement transcriptome loading
- Implement Salmon-format EC output writing
- Connect filtering → EC building → pruning pipeline

**Dependencies**:
- htslib (for BAM reading) or samtools API
- FASTA parser for transcriptome

**Example Integration Pattern**:
```cpp
// Pseudocode
1. Load transcriptome (FASTA) → get transcript names/lengths
2. For each read in BAM:
   a. Extract alignments
   b. Convert to RawAlignment[] format
   c. Apply updateRefMappings() for each alignment
   d. Apply filterAndCollectAlignments()
   e. Store filtered alignments per read
3. Build ECs: buildEquivalenceClasses(read_alignments, params, num_transcripts)
4. Apply extended pruning if enabled
5. Write ECs in Salmon format (eq_classes.txt)
```

### 2. Unit Test Execution

**Status**: Fixtures generated, but CLI needs to run tests

**Required Work**:
- Complete CLI to process test fixtures
- Run `generate_fixtures.py` → process → verify with `verify_filtering.py`

### 3. Parity Test Execution

**Status**: Scripts ready, but requires:
- BAM file generation (`generate_parity_bam.sh`)
- Working CLI integration
- Salmon binary available

**Required Work**:
- Run `generate_parity_bam.sh` to create `test/fixtures/salmon_eq/alignments.bam`
- Complete CLI integration
- Run `run_salmon_parity.sh` and verify >=99% parity

## Key Implementation Details

### Salmon Parameter Defaults

All defaults match `SalmonDefaults.hpp`:
- `score_exp = 1.0`
- `min_aln_prob = 1e-5`
- `decoy_threshold = 1.0`
- `hard_filter = false` (soft filter default)
- `range_factorization_bins = 4`
- `min_score_fraction = 0.65`

### Numerical Stability

- **Log-sum-exp**: Uses `logAdd(a, b) = a + log(1 + exp(b-a))` when `a >= b`
- **Score subtraction**: Always subtracts `bestScore` before `exp()` in `estAlnProb`
- **Log-space auxProbs**: Computed in log-space, normalized via log-sum-exp, then converted to probabilities

### Thread Safety

- Filtering functions are thread-safe (read-only params, per-read state)
- EC building uses hash map (thread-safe for separate reads)
- Extended pruning operates on final EC table (single-threaded)

### EC Format

**Salmon Weighted Format** (`--dumpEqWeights`):
```
<num_transcripts>
<num_ecs>
<transcript_name_1>
...
<transcript_name_N>
<k> <idx1> ... <idxk> <w1> ... <wk> <count>
```

**Our Output**: Must match this format exactly for parity testing.

## Testing Instructions

### Build Library

```bash
cd source/libem
make clean
make -j4
```

### Build CLI (after integration)

```bash
cd tools/ec_filter_test
make
```

### Run Unit Tests (after CLI integration)

```bash
cd tools/ec_filter_test
python3 generate_fixtures.py > fixtures.json
# Process fixtures with ec_filter_cli
python3 verify_filtering.py results.json
```

### Run Parity Test (after CLI integration)

```bash
cd tools/ec_filter_test

# Step 1: Generate BAM file
./generate_parity_bam.sh

# Step 2: Run parity test
./run_salmon_parity.sh /tmp/ec_parity_test

# Step 3: Check report
cat /tmp/ec_parity_test/parity_report.txt
```

**Expected Result**: >=99% parity on EC labels, weights, and counts.

## Integration Points

### With Existing Code

The EC filter module integrates with:
- `em_types.h` - Uses existing `EC` and `ECTable` structures
- `vb_engine.cpp` / `em_engine.cpp` - Output ECs can be fed directly to VB/EM engines

### With STAR

The module is designed to run **after alignment** when all threads are available:
- No nested OpenMP (alignment uses threads, EC filtering uses threads separately)
- Thread-local storage pattern matches `vb_engine.cpp` (TLS for parallel accumulation)

## File Structure

```
source/libem/
├── alignment_filter.h          ✅ Complete
├── alignment_filter.cpp        ✅ Complete
├── ec_builder.h                ✅ Complete
├── ec_builder.cpp              ✅ Complete
├── extended_pruning.h          ✅ Complete
├── extended_pruning.cpp        ✅ Complete
└── Makefile                    ✅ Updated

tools/ec_filter_test/
├── ec_filter_cli.cpp           ⚠️  Skeleton (needs BAM integration)
├── generate_fixtures.py        ✅ Complete
├── verify_filtering.py         ✅ Complete
├── generate_parity_bam.sh      ✅ Complete
├── run_salmon_parity.sh        ✅ Complete
├── compare_ecs.py              ✅ Complete
└── Makefile                    ✅ Complete
```

## References

### Salmon Source Files

Key logic matches:
- `SalmonMappingUtils.hpp:271-390` - `filterAndCollectAlignments`
- `SalmonMappingUtils.hpp:213-269` - `updateRefMappings`
- `SalmonQuantify.cpp:506-590` - auxProb and range factorization
- `SalmonDefaults.hpp` - Default parameter values

### Plan Document

See `.cursor/plans/ec_gating_cleanup_module_489701da.plan.md` for full specification.

## Questions for Next Agent

1. **BAM Reading Library**: Which library should be used? (htslib, samtools API, or other?)
2. **FASTA Parsing**: Use existing codebase parser or add dependency?
3. **EC Output Format**: Should we match Salmon's exact format or create our own first?
4. **Error Handling**: How should we handle malformed BAM/FASTA files?
5. **Performance**: Any specific performance requirements for large BAM files?

## Next Steps

1. **Immediate**: Complete CLI integration with BAM reading
2. **Short-term**: Run unit tests and verify all scenarios pass
3. **Medium-term**: Run parity test and achieve >=99% parity
4. **Long-term**: Integrate into STAR pipeline

---

**Last Updated**: Implementation complete, awaiting CLI integration
**Contact**: See plan document for full context
