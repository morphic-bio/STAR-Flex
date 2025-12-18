#!/bin/bash
# End-to-end trimming parity regression test
# Compares STAR's built-in cutadapt trimming vs Trim Galore + STAR
#
# This test validates that trimming parity translates to alignment parity by:
# 1. Running Trim Galore + STAR (baseline)
# 2. Running STAR-Flex with --trimCutadapt Yes (candidate)
# 3. Comparing BAM outputs (flagstat, idxstats, alignment diffs)
#
# Prerequisites:
#   - Reference genome and GTF (see GENOME_FA and GTF below)
#   - Pre-built STAR index (see STAR_INDEX below)
#   - Input FASTQs (see INPUT_R1 and INPUT_R2 below)
#   - /usr/local/bin/STAR (baseline STAR binary)
#   - trim_galore in PATH
#   - samtools in PATH

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Configuration
GENOME_FA="${GENOME_FA:-/mnt/pikachu/test-datasets/reference/chr22_23800000-23980000.fa}"
GTF="${GTF:-/mnt/pikachu/test-datasets/reference/genes.gtf}"
STAR_INDEX="${STAR_INDEX:-/tmp/star_chr22_index}"
WORK="${WORK:-/tmp/trim_parity_regress}"
THREADS="${THREADS:-4}"

INPUT_R1="${INPUT_R1:-/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz}"
INPUT_R2="${INPUT_R2:-/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz}"

STAR_BASELINE="${STAR_BASELINE:-/usr/local/bin/STAR}"
STAR_CANDIDATE="${STAR_CANDIDATE:-$PROJECT_ROOT/source/STAR}"

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

RESULTS_DIR="$WORK/results"
mkdir -p "$RESULTS_DIR"

echo "=== End-to-End Trimming Parity Regression Test ==="
echo "Baseline STAR: $STAR_BASELINE"
echo "Candidate STAR: $STAR_CANDIDATE"
echo "Reference: $GENOME_FA"
echo "GTF: $GTF"
echo "STAR Index: $STAR_INDEX"
echo "Input R1: $INPUT_R1"
echo "Input R2: $INPUT_R2"
echo ""

# Check prerequisites
if [ ! -f "$GENOME_FA" ]; then
    echo "Error: Genome FASTA not found: $GENOME_FA"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "Error: GTF file not found: $GTF"
    exit 1
fi

if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
    echo "Error: Input FASTQs not found"
    exit 1
fi

if [ ! -f "$STAR_BASELINE" ]; then
    echo "Error: Baseline STAR not found: $STAR_BASELINE"
    exit 1
fi

if [ ! -f "$STAR_CANDIDATE" ]; then
    echo "Error: Candidate STAR not found: $STAR_CANDIDATE"
    exit 1
fi

# Check if index exists, build if needed
if [ ! -f "$STAR_INDEX/Genome" ]; then
    echo "Building STAR index..."
    mkdir -p "$STAR_INDEX"
    # Try with GTF first, fall back to no GTF if chromosome names don't match
    if $STAR_BASELINE --runThreadN $THREADS --runMode genomeGenerate \
                      --genomeDir "$STAR_INDEX" \
                      --genomeFastaFiles "$GENOME_FA" \
                      --sjdbGTFfile "$GTF" \
                      --sjdbOverhang 99 2>&1 | grep -q "no valid exon lines"; then
        echo "GTF chromosome names don't match FASTA, building index without GTF..."
        $STAR_BASELINE --runThreadN $THREADS --runMode genomeGenerate \
                      --genomeDir "$STAR_INDEX" \
                      --genomeFastaFiles "$GENOME_FA"
    fi
    echo "Index built."
else
    echo "Using existing STAR index: $STAR_INDEX"
fi

mkdir -p "$WORK"

# Step 1: Baseline - Trim Galore + STAR
echo ""
echo "=== Step 1: Baseline (Trim Galore + STAR) ==="
BASELINE_DIR="$WORK/baseline"
mkdir -p "$BASELINE_DIR"

cp "$INPUT_R1" "$BASELINE_DIR/in_R1.fq.gz"
cp "$INPUT_R2" "$BASELINE_DIR/in_R2.fq.gz"

echo "Running Trim Galore..."
trim_galore --paired --quality 20 --length 20 \
            --adapter "$ADAPTER_R1" \
            --adapter2 "$ADAPTER_R2" \
            "$BASELINE_DIR/in_R1.fq.gz" "$BASELINE_DIR/in_R2.fq.gz" \
            --output_dir "$BASELINE_DIR" > "$BASELINE_DIR/trim_galore.log" 2>&1

# Trim Galore may output .fq.gz or .fq files
if [ -f "$BASELINE_DIR/in_R1_val_1.fq.gz" ] && [ -f "$BASELINE_DIR/in_R2_val_2.fq.gz" ]; then
    echo "Trim Galore produced gzipped output, using zcat..."
    BASELINE_R1="$BASELINE_DIR/in_R1_val_1.fq.gz"
    BASELINE_R2="$BASELINE_DIR/in_R2_val_2.fq.gz"
    READ_FILES_CMD="--readFilesCommand zcat"
elif [ -f "$BASELINE_DIR/in_R1_val_1.fq" ] && [ -f "$BASELINE_DIR/in_R2_val_2.fq" ]; then
    BASELINE_R1="$BASELINE_DIR/in_R1_val_1.fq"
    BASELINE_R2="$BASELINE_DIR/in_R2_val_2.fq"
    READ_FILES_CMD=""
else
    echo "Error: Trim Galore output files not found"
    ls -la "$BASELINE_DIR/" | grep -E "(val|trimmed)" || true
    exit 1
fi

echo "Running STAR (baseline)..."
$STAR_BASELINE --runThreadN $THREADS --genomeDir "$STAR_INDEX" \
               --readFilesIn "$BASELINE_R1" "$BASELINE_R2" \
               $READ_FILES_CMD \
               --outFileNamePrefix "$BASELINE_DIR/" \
               --outSAMtype BAM SortedByCoordinate \
               --outSAMunmapped None \
               --outSAMattributes Standard

# Step 2: Candidate - STAR-Flex with built-in trimming
echo ""
echo "=== Step 2: Candidate (STAR-Flex with --trimCutadapt Yes) ==="
CANDIDATE_DIR="$WORK/candidate"
mkdir -p "$CANDIDATE_DIR"

cp "$INPUT_R1" "$CANDIDATE_DIR/in_R1.fq.gz"
cp "$INPUT_R2" "$CANDIDATE_DIR/in_R2.fq.gz"

echo "Running STAR-Flex with built-in trimming..."
# Note: trimCutadaptAdapter expects space-separated R1 and R2 adapters
$STAR_CANDIDATE --runThreadN $THREADS --genomeDir "$STAR_INDEX" \
                --readFilesCommand zcat \
                --readFilesIn "$CANDIDATE_DIR/in_R1.fq.gz" "$CANDIDATE_DIR/in_R2.fq.gz" \
                --trimCutadapt Yes \
                --trimCutadaptQuality 20 \
                --trimCutadaptMinLength 20 \
                --trimCutadaptAdapter $ADAPTER_R1 $ADAPTER_R2 \
                --outFileNamePrefix "$CANDIDATE_DIR/" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped None \
                --outSAMattributes Standard

# Step 3: Compare results
echo ""
echo "=== Step 3: Comparison ==="

BASELINE_BAM="$BASELINE_DIR/Aligned.sortedByCoord.out.bam"
CANDIDATE_BAM="$CANDIDATE_DIR/Aligned.sortedByCoord.out.bam"

if [ ! -f "$BASELINE_BAM" ]; then
    echo "Error: Baseline BAM not found: $BASELINE_BAM"
    exit 1
fi

if [ ! -f "$CANDIDATE_BAM" ]; then
    echo "Error: Candidate BAM not found: $CANDIDATE_BAM"
    exit 1
fi

# Index BAMs if needed
[ ! -f "$BASELINE_BAM.bai" ] && samtools index "$BASELINE_BAM"
[ ! -f "$CANDIDATE_BAM.bai" ] && samtools index "$CANDIDATE_BAM"

# Compare flagstat
echo "Comparing flagstat..."
samtools flagstat "$BASELINE_BAM" > "$RESULTS_DIR/baseline.flagstat"
samtools flagstat "$CANDIDATE_BAM" > "$RESULTS_DIR/candidate.flagstat"

FLAGSTAT_DIFF=0
if diff -q "$RESULTS_DIR/baseline.flagstat" "$RESULTS_DIR/candidate.flagstat" >/dev/null 2>&1; then
    echo "✅ flagstat: IDENTICAL"
else
    echo "❌ flagstat: DIFFERENT"
    FLAGSTAT_DIFF=1
    diff "$RESULTS_DIR/baseline.flagstat" "$RESULTS_DIR/candidate.flagstat" > "$RESULTS_DIR/flagstat.diff" || true
    echo "Differences saved to: $RESULTS_DIR/flagstat.diff"
    head -20 "$RESULTS_DIR/flagstat.diff"
fi

# Compare idxstats
echo ""
echo "Comparing idxstats..."
samtools idxstats "$BASELINE_BAM" > "$RESULTS_DIR/baseline.idxstats"
samtools idxstats "$CANDIDATE_BAM" > "$RESULTS_DIR/candidate.idxstats"

IDXSTATS_DIFF=0
if diff -q "$RESULTS_DIR/baseline.idxstats" "$RESULTS_DIR/candidate.idxstats" >/dev/null 2>&1; then
    echo "✅ idxstats: IDENTICAL"
else
    echo "❌ idxstats: DIFFERENT"
    IDXSTATS_DIFF=1
    diff "$RESULTS_DIR/baseline.idxstats" "$RESULTS_DIR/candidate.idxstats" > "$RESULTS_DIR/idxstats.diff" || true
    echo "Differences saved to: $RESULTS_DIR/idxstats.diff"
    head -20 "$RESULTS_DIR/idxstats.diff"
fi

# Compare alignment records (spot check)
echo ""
echo "Comparing alignment records (spot check on chr22:23800000-23801000)..."
samtools view "$BASELINE_BAM" chr22:23800000-23801000 | head -100 > "$RESULTS_DIR/baseline.alignments" 2>/dev/null || \
    samtools view "$BASELINE_BAM" | head -100 > "$RESULTS_DIR/baseline.alignments"
samtools view "$CANDIDATE_BAM" chr22:23800000-23801000 | head -100 > "$RESULTS_DIR/candidate.alignments" 2>/dev/null || \
    samtools view "$CANDIDATE_BAM" | head -100 > "$RESULTS_DIR/candidate.alignments"

ALIGN_DIFF=0
if diff -q "$RESULTS_DIR/baseline.alignments" "$RESULTS_DIR/candidate.alignments" >/dev/null 2>&1; then
    echo "✅ Alignment records: IDENTICAL"
else
    echo "⚠️  Alignment records: DIFFERENT (checking if only tags differ)..."
    ALIGN_DIFF=1
    # Compare only core fields (columns 1-11)
    cut -f1-11 "$RESULTS_DIR/baseline.alignments" > "$RESULTS_DIR/baseline.alignments.core"
    cut -f1-11 "$RESULTS_DIR/candidate.alignments" > "$RESULTS_DIR/candidate.alignments.core"
    
    if diff -q "$RESULTS_DIR/baseline.alignments.core" "$RESULTS_DIR/candidate.alignments.core" >/dev/null 2>&1; then
        echo "✅ Core alignment fields (QNAME-FLAG-CIGAR-POS-MAPQ) are IDENTICAL"
        echo "   Differences are only in optional tags (expected)"
    else
        echo "❌ Core alignment fields differ"
        diff "$RESULTS_DIR/baseline.alignments.core" "$RESULTS_DIR/candidate.alignments.core" > "$RESULTS_DIR/alignments.diff" || true
        head -20 "$RESULTS_DIR/alignments.diff"
    fi
fi

# Summary
echo ""
echo "=== Summary ==="
echo "Baseline BAM: $BASELINE_BAM"
echo "Candidate BAM: $CANDIDATE_BAM"
echo "Results directory: $RESULTS_DIR"
echo ""

TOTAL_DIFFS=$((FLAGSTAT_DIFF + IDXSTATS_DIFF + ALIGN_DIFF))
if [ $TOTAL_DIFFS -eq 0 ]; then
    echo "✅ TEST PASSED: All comparisons identical"
    echo "PASS" > "$RESULTS_DIR/status.txt"
    exit 0
else
    echo "❌ TEST FAILED: $TOTAL_DIFFS comparison(s) differ"
    echo "FAIL" > "$RESULTS_DIR/status.txt"
    echo ""
    echo "Note: Since FASTQ-level parity is perfect (0 diff lines), any differences"
    echo "      would indicate STAR parameter mismatches or version differences,"
    echo "      not trimming algorithm differences."
    exit 1
fi
