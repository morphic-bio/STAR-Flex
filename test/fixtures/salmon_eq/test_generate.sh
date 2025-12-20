#!/bin/bash
# Test wrapper for GENERATE.sh using Docker for tools
# This tests the script logic using Docker containers

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Testing GENERATE.sh with Docker ==="
echo ""

# Test Step 1: GTF subsetting
echo "Test Step 1: GTF subsetting..."
REF_GTF_SOURCE="/mnt/pikachu/genome.gtf"
REF_GTF="/tmp/genes_chr22.gtf"

if echo "${REF_GTF_SOURCE}" | grep -q "\.gz$"; then
    gunzip -c "${REF_GTF_SOURCE}" 2>/dev/null | \
        awk '$1=="chr22" && $4>=23800000 && $5<=23980000' > "${REF_GTF}" || \
    gunzip -c "${REF_GTF_SOURCE}" 2>/dev/null | \
        awk '$1=="chr22"' > "${REF_GTF}"
else
    awk '$1=="chr22" && $4>=23800000 && $5<=23980000' "${REF_GTF_SOURCE}" > "${REF_GTF}" 2>/dev/null || \
    awk '$1=="chr22"' "${REF_GTF_SOURCE}" > "${REF_GTF}"
fi

if [ ! -s "${REF_GTF}" ]; then
    echo "ERROR: GTF subsetting failed"
    exit 1
fi
echo "✓ GTF subsetting: $(wc -l < ${REF_GTF} | xargs) lines"

# Test Step 2: Build transcriptome (using Docker for gffread)
echo ""
echo "Test Step 2: Building transcriptome FASTA..."
REF_FASTA="/mnt/pikachu/test-datasets/reference/chr22_23800000-23980000.fa"
TRANSCRIPTOME="/tmp/chr22_trans.fa"

docker run --rm \
    -v /mnt/pikachu:/mnt/pikachu:ro \
    -v /tmp:/tmp \
    quay.io/biocontainers/gffread:0.12.7--h8b12597_0 \
    gffread "${REF_GTF}" -g "${REF_FASTA}" -w "${TRANSCRIPTOME}"

if [ ! -s "${TRANSCRIPTOME}" ]; then
    echo "ERROR: Transcriptome creation failed"
    exit 1
fi
echo "✓ Transcriptome: $(wc -l < ${TRANSCRIPTOME} | xargs) lines"

# Test Step 3: Build Salmon index
echo ""
echo "Test Step 3: Building Salmon index..."
SALMON_IDX="/tmp/salmon_chr22_idx"

docker run --rm \
    -v /tmp:/tmp \
    combinelab/salmon:latest \
    salmon index -t "${TRANSCRIPTOME}" -i "${SALMON_IDX}"

if [ ! -d "${SALMON_IDX}" ]; then
    echo "ERROR: Salmon index creation failed"
    exit 1
fi
echo "✓ Salmon index created"

# Test Step 4: Prepare reads
echo ""
echo "Test Step 4: Preparing input reads..."
SYNTHETIC_R1="${SCRIPT_DIR}/synthetic_reads1.fq"
SYNTHETIC_R2="${SCRIPT_DIR}/synthetic_reads2.fq"
SUB_R1="/tmp/sub_R1.fq"
SUB_R2="/tmp/sub_R2.fq"

if [ -f "${SYNTHETIC_R1}" ] && [ -f "${SYNTHETIC_R2}" ]; then
    cp "${SYNTHETIC_R1}" "${SUB_R1}"
    cp "${SYNTHETIC_R2}" "${SUB_R2}"
    echo "✓ Using synthetic reads: $(wc -l < ${SUB_R1} | xargs) lines"
else
    echo "ERROR: Synthetic reads not found"
    exit 1
fi

# Test Step 5: Run Salmon quant
echo ""
echo "Test Step 5: Running Salmon quant..."
SALMON_OUT="/tmp/salmon_eq_test"

docker run --rm \
    -v /tmp:/tmp \
    combinelab/salmon:latest \
    salmon quant -i "${SALMON_IDX}" \
        -l A \
        -1 "${SUB_R1}" \
        -2 "${SUB_R2}" \
        --dumpEq \
        --threads 4 \
        --noLengthCorrection \
        --noFragLengthDist \
        --noEffectiveLengthCorrection \
        -o "${SALMON_OUT}"

# Test Step 6: Check output files
echo ""
echo "Test Step 6: Checking output files..."
if [ -f "${SALMON_OUT}/aux_info/eq_classes.txt.gz" ]; then
    gunzip -c "${SALMON_OUT}/aux_info/eq_classes.txt.gz" > /tmp/test_eq_classes.txt
    echo "✓ eq_classes.txt.gz extracted: $(wc -l < /tmp/test_eq_classes.txt | xargs) lines"
elif [ -f "${SALMON_OUT}/aux_info/eq_classes.txt" ]; then
    echo "✓ eq_classes.txt: $(wc -l < ${SALMON_OUT}/aux_info/eq_classes.txt | xargs) lines"
else
    echo "ERROR: eq_classes.txt not found"
    exit 1
fi

if [ -f "${SALMON_OUT}/quant.sf" ]; then
    echo "✓ quant.sf: $(wc -l < ${SALMON_OUT}/quant.sf | xargs) lines"
    head -5 "${SALMON_OUT}/quant.sf"
else
    echo "ERROR: quant.sf not found"
    exit 1
fi

echo ""
echo "=== All tests passed! ==="
