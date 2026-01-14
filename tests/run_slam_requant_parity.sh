#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

STAR_BIN="${STAR_BIN:-$ROOT_DIR/source/STAR}"
REQ_BIN="${REQ_BIN:-$ROOT_DIR/tools/slam_requant/slam_requant}"
FASTQ="${FASTQ:-$ROOT_DIR/test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz}"
STAR_INDEX="${STAR_INDEX:-$ROOT_DIR/test/fixtures/slam/ref/star_index}"
SNPS_BED="${SNPS_BED:-$ROOT_DIR/test/fixtures/slam/ref/snps.bed}"

WORK="${WORK:-$ROOT_DIR/test/tmp_slam_requant}"
OUT_PREFIX="${OUT_PREFIX:-$WORK/fixture_}"
DUMP_PATH="${DUMP_PATH:-$WORK/fixture.dump.bin}"

COMPARE="${COMPARE:-$SCRIPT_DIR/slam/compare_star_outputs.py}"

mkdir -p "$WORK"

if [[ ! -x "$STAR_BIN" || ! -x "$REQ_BIN" ]]; then
  echo "Missing STAR or slam_requant binary"
  exit 1
fi

echo "=== STAR: build dump + SlamQuant ==="
"$STAR_BIN" \
  --runThreadN 4 \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$FASTQ" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT_PREFIX" \
  --outSAMtype None \
  --clip3pAdapterSeq AGATCGGAAGAG \
  --clip3pAdapterMMp 0.1 \
  --slamQuantMode 1 \
  --slamSnpMaskIn "$SNPS_BED" \
  --slamDumpBinary "$DUMP_PATH" \
  > "${OUT_PREFIX}slam.log" 2>&1

echo "=== slam_requant: replay dump ==="
"$REQ_BIN" \
  --dump "$DUMP_PATH" \
  --out "$WORK/requant_" \
  --slamSnpMaskIn "$SNPS_BED"

echo "=== Compare STAR vs slam_requant ==="
python3 "$COMPARE" \
  --reference "${OUT_PREFIX}SlamQuant.out" \
  --test "$WORK/requant_SlamQuant.out"
