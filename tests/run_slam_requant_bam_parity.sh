#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

REQ_BIN="${REQ_BIN:-$ROOT_DIR/tools/slam_requant/slam_requant}"
COMPARE="${COMPARE:-$SCRIPT_DIR/slam/compare_star_outputs.py}"

# Required inputs (point to existing BAMs/annotations)
BAM="${BAM:?Set BAM=/path/to.bam}"
GTF="${GTF:?Set GTF=/path/to.gtf}"
FASTA="${FASTA:?Set FASTA=/path/to.fa}"

# Optional mask
SNPS_BED="${SNPS_BED:-}"

WORK="${WORK:-$ROOT_DIR/test/tmp_slam_requant_bam}"
STAR_REF="${STAR_REF:?Set STAR_REF=/path/to/STAR_SlamQuant.out}"
OUT_PREFIX="${OUT_PREFIX:-$WORK/requant_}"

mkdir -p "$WORK"

ARGS=(--bam "$BAM" --gtf "$GTF" --fasta "$FASTA" --out "$OUT_PREFIX")
if [[ -n "$SNPS_BED" ]]; then
  ARGS+=(--slamSnpMaskIn "$SNPS_BED")
fi

echo "=== slam_requant (BAM path) ==="
"$REQ_BIN" "${ARGS[@]}"

echo "=== Compare STAR vs slam_requant ==="
python3 "$COMPARE" \
  --reference "$STAR_REF" \
  --test "${OUT_PREFIX}SlamQuant.out"
