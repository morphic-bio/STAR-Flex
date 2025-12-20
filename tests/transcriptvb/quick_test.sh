#!/bin/bash
# quick_test.sh - Run basic TranscriptVB validation tests
#
# Usage: ./quick_test.sh [STAR_BIN] [GENOME_DIR] [READS_1] [READS_2]
#
# If arguments not provided, uses defaults for STAR-Flex test data

set -e

# Default paths (adjust for your environment)
STAR_BIN=${1:-/mnt/pikachu/STAR-Flex/source/STAR}
GENOME_DIR=${2:-/tmp/star_vb_test/star_new_index}
READS_1=${3:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz}
READS_2=${4:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz}

OUTDIR=/tmp/transcriptvb_quick_test_$$
PASSED=0
FAILED=0

mkdir -p $OUTDIR
cd $OUTDIR

echo "========================================"
echo "TranscriptVB Quick Validation Tests"
echo "========================================"
echo "STAR: $STAR_BIN"
echo "Genome: $GENOME_DIR"
echo "Reads: $READS_1, $READS_2"
echo "Output: $OUTDIR"
echo "========================================"
echo ""

# Test 1: Basic TranscriptVB
echo "=== Test 1: Basic TranscriptVB ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 4 \
    --outFileNamePrefix basic_ \
    2>&1 | tail -5

if [ -f basic_quant.sf ]; then
    echo "✓ quant.sf created"
    
    # Check TPM sum
    TPM_SUM=$(tail -n +2 basic_quant.sf | awk '{sum+=$4} END {printf "%.0f", sum}')
    echo "  TPM sum: $TPM_SUM (expected ~1000000)"
    
    # Check for expressed transcripts
    EXPRESSED=$(tail -n +2 basic_quant.sf | awk '$5 > 0 {count++} END {print count}')
    echo "  Expressed transcripts: $EXPRESSED"
    
    # Check convergence
    if grep -q "converged: yes" basic_Log.out 2>/dev/null; then
        echo "  ✓ Quantification converged"
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ Convergence status unknown"
        FAILED=$((FAILED + 1))
    fi
else
    echo "✗ quant.sf NOT created"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 2: GC Bias Collection
echo "=== Test 2: GC Bias Collection ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBgcBias 1 \
    --runThreadN 4 \
    --outFileNamePrefix gc_ \
    2>&1 | tail -5

if grep -q "GC bias: collected" gc_Log.out 2>/dev/null; then
    GC_OBS=$(grep "GC bias: collected" gc_Log.out | sed -n 's/.*collected \([0-9]*\) fragment.*/\1/p')
    echo "✓ GC observations collected: $GC_OBS"
    PASSED=$((PASSED + 1))
else
    echo "✗ No GC observations found in log"
    FAILED=$((FAILED + 1))
fi

if grep -q "FLD-adjusted" gc_Log.out 2>/dev/null; then
    echo "✓ FLD-adjusted effective lengths"
else
    echo "  (FLD adjustment not logged)"
fi
echo ""

# Test 3: EM Mode
echo "=== Test 3: EM Mode ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBem 1 \
    --runThreadN 4 \
    --outFileNamePrefix em_ \
    2>&1 | tail -5

if [ -f em_quant.sf ]; then
    echo "✓ EM mode completed"
    PASSED=$((PASSED + 1))
else
    echo "✗ EM mode failed"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 4: Single Thread Determinism
echo "=== Test 4: Single Thread Determinism ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 1 \
    --outFileNamePrefix st1_ \
    2>&1 | tail -3

$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 1 \
    --outFileNamePrefix st2_ \
    2>&1 | tail -3

if diff -q st1_quant.sf st2_quant.sf > /dev/null 2>&1; then
    echo "✓ Deterministic results (single thread)"
    PASSED=$((PASSED + 1))
else
    echo "✗ Non-deterministic results"
    FAILED=$((FAILED + 1))
fi
echo ""

# Summary
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All tests passed!"
    echo ""
    echo "Output files in: $OUTDIR"
    exit 0
else
    echo "✗ Some tests failed"
    echo ""
    echo "Check output in: $OUTDIR"
    exit 1
fi

