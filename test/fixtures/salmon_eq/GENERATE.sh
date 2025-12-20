#!/bin/bash
# Script to generate Salmon equivalence class fixture
# Requires: gffread, salmon
#
# This script:
# 1. Subsets human GTF to chr22 (and optional window)
# 2. Builds transcriptome FASTA from subset GTF + chr22 FASTA slice
# 3. Builds Salmon index
# 4. Downsamples input FASTQ files (or uses existing synthetic reads)
# 5. Runs Salmon quant with bias disabled (for v1 parity testing)
# 6. Copies eq_classes.txt and quant.sf to fixture directory
#
# Usage:
#   cd test/fixtures/salmon_eq
#   ./GENERATE.sh
#
# Or with custom tool paths:
#   GFFREAD=/path/to/gffread SALMON=/path/to/salmon ./GENERATE.sh
#
# Note: If synthetic_reads1.fq and synthetic_reads2.fq exist in the current
# directory, they will be used instead of downsampling from SRR6357070.

set -e

# Allow custom tool paths via environment variables
GFFREAD="${GFFREAD:-gffread}"
SALMON="${SALMON:-salmon}"

# Check if tools are available
if ! command -v "${GFFREAD}" &> /dev/null; then
    echo "Error: gffread not found. Install with: conda install -c bioconda gffread"
    echo "       Or set GFFREAD=/path/to/gffread"
    exit 1
fi

if ! command -v "${SALMON}" &> /dev/null; then
    echo "Error: salmon not found. Install with: conda install -c bioconda salmon"
    echo "       Or set SALMON=/path/to/salmon"
    exit 1
fi

REF_FASTA="/mnt/pikachu/test-datasets/reference/chr22_23800000-23980000.fa"
# Try human GTF first, fall back to test-datasets GTF
if [ -f "/mnt/pikachu/genome.gtf.gz" ]; then
    REF_GTF_SOURCE="/mnt/pikachu/genome.gtf.gz"
elif [ -f "/mnt/pikachu/genome.gtf" ]; then
    REF_GTF_SOURCE="/mnt/pikachu/genome.gtf"
else
    REF_GTF_SOURCE="/mnt/pikachu/test-datasets/reference/genes.gtf"
fi
R1_FASTQ="/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz"
R2_FASTQ="/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz"

TMPDIR="${TMPDIR:-/tmp}"
REF_GTF="${TMPDIR}/genes_chr22.gtf"
TRANSCRIPTOME="${TMPDIR}/chr22_trans.fa"
SALMON_IDX="${TMPDIR}/salmon_chr22_idx"
SUB_R1="${TMPDIR}/sub_R1.fq"
SUB_R2="${TMPDIR}/sub_R2.fq"
SALMON_OUT="${TMPDIR}/salmon_eq"
OUTDIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Generating Salmon Equivalence Class Fixture ==="
echo "Using: gffread=$(which ${GFFREAD}), salmon=$(which ${SALMON})"
echo ""

# Step 1: Subset GTF to chr22 (and optional window)
echo "Step 1: Subsetting GTF to chr22..."
if echo "${REF_GTF_SOURCE}" | grep -q "\.gz$"; then
    # GTF is gzipped - try with coordinate filter first (23.8-23.98 Mb window)
    gunzip -c "${REF_GTF_SOURCE}" 2>/dev/null | \
        awk '$1=="chr22" && $4>=23800000 && $5<=23980000' > "${REF_GTF}" || \
    # If no records, use all chr22
    gunzip -c "${REF_GTF_SOURCE}" 2>/dev/null | \
        awk '$1=="chr22"' > "${REF_GTF}"
else
    # GTF is uncompressed - try with coordinate filter first
    awk '$1=="chr22" && $4>=23800000 && $5<=23980000' "${REF_GTF_SOURCE}" > "${REF_GTF}" 2>/dev/null || \
    # If no records, use all chr22
    awk '$1=="chr22"' "${REF_GTF_SOURCE}" > "${REF_GTF}"
fi
if [ ! -s "${REF_GTF}" ]; then
    echo "Error: No chr22 features found in GTF. Check GTF format (chromosome names should be 'chr22')."
    exit 1
fi
echo "  Created: ${REF_GTF} ($(wc -l < ${REF_GTF} | xargs) lines)"

# Step 2: Build transcriptome FASTA from subset GTF
echo ""
echo "Step 2: Building transcriptome FASTA from subset GTF..."
# Try gffread first
if "${GFFREAD}" "${REF_GTF}" -g "${REF_FASTA}" -w "${TRANSCRIPTOME}" 2>/dev/null && [ -s "${TRANSCRIPTOME}" ]; then
    echo "  Created with gffread: ${TRANSCRIPTOME} ($(wc -l < ${TRANSCRIPTOME} | xargs) lines)"
else
    echo "  gffread failed (likely coordinate issues with sliced FASTA), using Python fallback..."
    PYTHON_SCRIPT="${OUTDIR}/build_transcriptome.py"
    if [ ! -f "${PYTHON_SCRIPT}" ]; then
        echo "Error: Python script not found: ${PYTHON_SCRIPT}"
        exit 1
    fi
    python3 "${PYTHON_SCRIPT}" "${REF_GTF}" "${REF_FASTA}" "${TRANSCRIPTOME}" 23800000
    if [ ! -s "${TRANSCRIPTOME}" ]; then
        echo "Error: Failed to create transcriptome FASTA with Python script."
        exit 1
    fi
    echo "  Created with Python script: ${TRANSCRIPTOME} ($(wc -l < ${TRANSCRIPTOME} | xargs) lines)"
fi

# Step 3: Build Salmon index
echo ""
echo "Step 3: Building Salmon index..."
"${SALMON}" index -t "${TRANSCRIPTOME}" -i "${SALMON_IDX}"
echo "  Created: ${SALMON_IDX}/"

# Step 4: Prepare input reads (use synthetic reads if available, otherwise downsample)
echo ""
echo "Step 4: Preparing input reads..."
SYNTHETIC_R1="${OUTDIR}/synthetic_reads1.fq"
SYNTHETIC_R2="${OUTDIR}/synthetic_reads2.fq"

if [ -f "${SYNTHETIC_R1}" ] && [ -f "${SYNTHETIC_R2}" ]; then
    echo "  Using existing synthetic reads..."
    cp "${SYNTHETIC_R1}" "${SUB_R1}"
    cp "${SYNTHETIC_R2}" "${SUB_R2}"
    echo "  Copied: ${SUB_R1} ($(wc -l < ${SUB_R1} | xargs) lines)"
    echo "  Copied: ${SUB_R2} ($(wc -l < ${SUB_R2} | xargs) lines)"
elif [ -f "${R1_FASTQ}" ] && [ -f "${R2_FASTQ}" ]; then
    echo "  Downsampling FASTQ files (1000 read pairs = 4000 FASTQ lines)..."
    if echo "${R1_FASTQ}" | grep -q "\.gz$"; then
        zcat "${R1_FASTQ}" | head -n 4000 > "${SUB_R1}"
        zcat "${R2_FASTQ}" | head -n 4000 > "${SUB_R2}"
    else
        head -n 4000 "${R1_FASTQ}" > "${SUB_R1}"
        head -n 4000 "${R2_FASTQ}" > "${SUB_R2}"
    fi
    echo "  Created: ${SUB_R1} ($(wc -l < ${SUB_R1} | xargs) lines)"
    echo "  Created: ${SUB_R2} ($(wc -l < ${SUB_R2} | xargs) lines)"
else
    echo "Error: No input reads found."
    echo "  Expected one of:"
    echo "    - Synthetic reads: ${SYNTHETIC_R1} / ${SYNTHETIC_R2}"
    echo "    - Real reads: ${R1_FASTQ} / ${R2_FASTQ}"
    exit 1
fi

# Step 5: Run Salmon quant with bias disabled
echo ""
echo "Step 5: Running Salmon quant with --dumpEqWeights (bias disabled for parity testing)..."
# IMPORTANT: Use --dumpEqWeights to get per-EC per-transcript combinedWeights (auxs)
# Without weights, we can only achieve ~85% parity due to missing alignment quality info
# Also use --threads 1 for deterministic floating-point order
# Our EM engine uses raw transcript lengths, not FLD-adjusted effective lengths
if [ ! -f "${SUB_R1}" ] || [ ! -f "${SUB_R2}" ]; then
    echo "Error: Downsampled FASTQ files not found: ${SUB_R1} / ${SUB_R2}"
    echo "  Ensure input FASTQ files exist or synthetic reads are already in place."
    exit 1
fi
"${SALMON}" quant -i "${SALMON_IDX}" \
             -l A \
             -1 "${SUB_R1}" \
             -2 "${SUB_R2}" \
             --dumpEqWeights \
             --threads 1 \
             --noLengthCorrection \
             --noFragLengthDist \
             --noEffectiveLengthCorrection \
             -o "${SALMON_OUT}"

echo ""
echo "Step 6: Copying artifacts..."
# Handle gzipped eq_classes.txt (Salmon v1.10+)
if [ -f "${SALMON_OUT}/aux_info/eq_classes.txt.gz" ]; then
    gunzip -c "${SALMON_OUT}/aux_info/eq_classes.txt.gz" > "${OUTDIR}/eq_classes.txt"
    echo "  Extracted: eq_classes.txt (from .gz)"
elif [ -f "${SALMON_OUT}/aux_info/eq_classes.txt" ]; then
    cp "${SALMON_OUT}/aux_info/eq_classes.txt" "${OUTDIR}/"
    echo "  Copied: eq_classes.txt"
else
    echo "Error: eq_classes.txt not found in ${SALMON_OUT}/aux_info/"
    exit 1
fi

cp "${SALMON_OUT}/quant.sf" "${OUTDIR}/"
echo "  Copied: quant.sf"

# Copy transcript map if it exists
if [ -f "${SALMON_OUT}/aux_info/transcript_map.txt" ]; then
    cp "${SALMON_OUT}/aux_info/transcript_map.txt" "${OUTDIR}/"
    echo "  Copied: transcript_map.txt"
elif [ -f "${SALMON_OUT}/aux_info/transcript_to_gene_map.tsv" ]; then
    cp "${SALMON_OUT}/aux_info/transcript_to_gene_map.tsv" "${OUTDIR}/"
    echo "  Copied: transcript_to_gene_map.tsv"
fi

echo ""
echo "=== Done! Files copied to ${OUTDIR} ==="
echo ""
echo "File sizes:"
du -h "${OUTDIR}"/eq_classes.txt "${OUTDIR}"/quant.sf "${OUTDIR}"/*.txt "${OUTDIR}"/*.tsv 2>/dev/null | sort -h
echo ""
echo "Line counts:"
wc -l "${OUTDIR}"/eq_classes.txt "${OUTDIR}"/quant.sf 2>/dev/null || true
echo ""
echo "Sample equivalence classes (first 5 lines):"
head -5 "${OUTDIR}"/eq_classes.txt
echo ""
echo "Sample quant.sf (first 5 lines):"
head -5 "${OUTDIR}"/quant.sf
