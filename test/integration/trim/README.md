# Integration Tests for Trimming

This directory contains integration tests for STAR-Flex's cutadapt-style trimming feature.

## Tests

### 1. FASTQ Parity Tests (`run_parity_test.sh`)

Validates that trimmed FASTQ outputs match Trim Galore/cutadapt exactly.

**Run**: `cd ../../source && make test_trim_parity`

**Status**: ✅ All 9 tests pass, 0 diff lines

### 2. End-to-End Alignment Parity Test (`run_e2e_trim_parity.sh`)

Validates that STAR alignments using built-in trimming match Trim Galore + STAR baseline.

**Run**: `cd ../../source && make test_trim_e2e`

**Status**: ✅ All comparisons identical

## Test Fixtures

- `nfcore_smoke/` - Small synthetic test dataset
- `nfcore_real/` - Real nf-core RNA-seq data (GSE110004 SRR6357070)

## Documentation

- `E2E_TEST_SUMMARY.md` - End-to-end test documentation
- `E2E_PARITY_RESULTS.md` - Test results
