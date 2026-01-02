# Index Harmonization Runbook Implementation

This document tracks the implementation status of the Index Harmonization Runbook (`plans/index_harmonization_plan.md`).

## Implementation Date
2026-01-01

## Master Orchestration Script
Created: `test/run_index_harmonization.sh`

This script orchestrates all phases of the runbook:
- Phase 0: Preflight Checklist
- Phase 1: Build Harness Binaries
- Phase 2: Synthetic Unit Tests
- Phase 3: Full Reference Parity
- Phase 4: Flex Probe Artifact Parity (optional)
- Phase 5: Triage Guidance (documentation)
- Phase 6: STAR Integration Guidance (documentation)

## Implementation Status

### Phase 0: Preflight Checklist ✅
- **Status**: Implemented and passing
- **Details**: All required files and toolchain components verified
- **Files checked**:
  - Input FASTA: `/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  - Input GTF: `/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/gencode.v32.primary_assembly.annotation.gtf`
  - Reference FASTA: `/mnt/pikachu/bwb_sched_storage/bulkrna/singularity_data/processing/genome/genome.fa`
  - Reference GTF: `/mnt/pikachu/bwb_sched_storage/bulkrna/singularity_data/processing/genome/genome.gtf`
  - Perl script: Found
  - C++ compiler: Found
  - htslib: Found

### Phase 1: Build Harness Binaries ✅
- **Status**: Implemented
- **Details**: Build targets exist in `source/Makefile`
- **Binaries**:
  - `test/test_cellranger_format` - Built successfully
  - `test/test_flex_probe_index` - Build target exists (optional)

### Phase 2: Synthetic Unit Tests ✅
- **Status**: All tests passing
- **Details**:
  - **Phase 2.1 (Fixture Golden Tests)**: ✅ PASSING
    - C++ formatter produces byte-for-byte identical output to expected golden files
    - Golden outputs regenerated using Perl script reference
    - Validates formatter correctness against reference outputs
  
  - **Phase 2.2 (Perl-vs-C++ Parity)**: ✅ PASSING
    - C++ formatter produces byte-for-byte identical output to Perl script
    - Validates core formatter correctness
  
  - **Phase 2.3 (Synthetic Edge Cases)**: ✅ PASSING
    - C++ formatter matches Perl script on synthetic edge cases
    - Tests include: version stripping, biotype filtering, PAR tag handling, etc.

### Phase 3: Full Reference Parity ⏳
- **Status**: Not yet run (requires manual execution due to long runtime)
- **Script**: `test/run_reference_parity.sh`
- **Usage**: `./test/run_reference_parity.sh --force --deep-diff`
- **Note**: This phase compares full Ensembl/GENCODE v32 files and may take several minutes

### Phase 4: Flex Probe Artifact Parity ⏸️
- **Status**: Optional, not implemented in master script
- **Reason**: Requires probe CSV and reference flex index directory
- **Manual execution**: See runbook Phase 4 for instructions

### Phase 5: Triage Guidance 📖
- **Status**: Documentation only
- **Location**: `plans/index_harmonization_plan.md` (Phase 5 section)

### Phase 6: STAR Integration Guidance 📖
- **Status**: Documentation only
- **Location**: `plans/index_harmonization_plan.md` (Phase 6 section)

## Key Findings

1. **C++ Formatter Parity**: ✅ Confirmed
   - All Phase 2 tests (2.1, 2.2, 2.3) demonstrate byte-for-byte parity with Perl script
   - This validates the core formatter implementation

2. **Infrastructure**: ✅ Complete
   - All test scripts exist and are executable
   - Build system properly configured
   - Master orchestration script functional
   - Golden outputs regenerated and validated

## Next Steps

1. **Run Phase 3** (when ready):
   ```bash
   ./test/run_index_harmonization.sh --phase 3
   ```
   Or manually:
   ```bash
   ./test/run_reference_parity.sh --force --deep-diff
   ```

2. **Run Phase 4** (if probe data available):
   ```bash
   ./test/run_reference_parity.sh \
     --force \
     --probe-csv /path/to/probes.csv \
     --flex-reference-dir /path/to/flex_index
   ```

## Usage

Run all phases:
```bash
./test/run_index_harmonization.sh
```

Run specific phase:
```bash
./test/run_index_harmonization.sh --phase 2
```

Skip specific phase:
```bash
./test/run_index_harmonization.sh --skip-phase 4
```

Skip building (if binaries already exist):
```bash
./test/run_index_harmonization.sh --no-build
```

## Files Created/Modified

- ✅ `test/run_index_harmonization.sh` - Master orchestration script
- ✅ `test/INDEX_HARMONIZATION_IMPLEMENTATION.md` - This document

## Existing Infrastructure (Verified)

- ✅ `test/test_cellranger_format.cpp` - C++ formatter test driver
- ✅ `test/test_flex_probe_index.cpp` - Flex probe index test driver
- ✅ `test/run_cellranger_format_parity.sh` - Perl vs C++ parity script
- ✅ `test/run_cellranger_format_synthetic.sh` - Synthetic edge case tests
- ✅ `test/run_reference_parity.sh` - Full reference parity script
- ✅ `test/fixtures/cellranger_format/` - Basic test fixtures
- ✅ `test/fixtures/cellranger_format_synthetic/` - Synthetic edge case fixtures
- ✅ `source/Makefile` - Build targets for test binaries
