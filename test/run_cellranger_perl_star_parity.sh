#!/bin/bash
#
# run_cellranger_perl_star_parity.sh
# 1) Run Perl formatter to create CellRanger-style outputs and compare to fixtures.
# 2) If Perl matches, run STAR internal path and compare outputs to the same fixtures.
#
# Defaults are set for GENCODE v44 inputs and PE fixtures provided by the user.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DEFAULT_INPUT_FASTA="/mnt/pikachu/autoindex_110_44/cellranger_ref_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DEFAULT_INPUT_GTF="/mnt/pikachu/autoindex_110_44/cellranger_ref_cache/gencode.v44.primary_assembly.annotation.gtf"
DEFAULT_REFERENCE_DIR="/mnt/pikachu/bwb_sched_storage/bulkrna/singularity_data/processing/genome"
DEFAULT_PERL_SCRIPT="/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl"
DEFAULT_STAR_BIN="$REPO_ROOT/source/STAR"
DEFAULT_WORK_DIR="$REPO_ROOT/test/tmp_cr_perl_star_parity"

INPUT_FASTA="$DEFAULT_INPUT_FASTA"
INPUT_GTF="$DEFAULT_INPUT_GTF"
REFERENCE_DIR="$DEFAULT_REFERENCE_DIR"
PERL_SCRIPT="$DEFAULT_PERL_SCRIPT"
STAR_BIN="$DEFAULT_STAR_BIN"
WORK_DIR="$DEFAULT_WORK_DIR"
REUSE_STAR_DIR=""
FORCE=false
SKIP_STAR=false
SKIP_PERL=false
THREADS=8

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
  --input-fasta PATH       Raw FASTA input (default: $DEFAULT_INPUT_FASTA)
  --input-gtf PATH         Raw GTF input (default: $DEFAULT_INPUT_GTF)
  --reference-dir PATH     Fixture dir with genome.fa + genome.gtf (default: $DEFAULT_REFERENCE_DIR)
  --perl-script PATH       Perl formatter path (default: $DEFAULT_PERL_SCRIPT)
  --star-bin PATH          STAR binary path (default: $DEFAULT_STAR_BIN)
  --work-dir PATH          Working directory (default: $DEFAULT_WORK_DIR)
  --reuse-star-dir PATH    Reuse existing STAR genomeDir or cellranger_ref outputs (skip STAR run)
  --threads N              Threads for STAR (default: 8)
  --skip-star              Skip STAR step (only Perl vs fixtures)
  --skip-perl              Skip Perl step (only STAR vs fixtures)
  --force                  Re-run even if outputs exist
  -h, --help               Show help
EOF
    exit "${1:-0}"
}

log() {
    echo "[cr_perl_star_parity] $*" >&2
}

error_exit() {
    echo "[cr_perl_star_parity] ERROR: $*" >&2
    exit 1
}

compare_files() {
    local label="$1"
    local ref="$2"
    local out="$3"

    if cmp -s "$ref" "$out"; then
        log "✓ ${label} matches"
        return 0
    fi

    log "✗ ${label} differs"
    log "  ref: $ref"
    log "  out: $out"
    diff -u "$ref" "$out" | head -n 20 || true
    return 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-fasta) INPUT_FASTA="$2"; shift 2 ;;
        --input-gtf) INPUT_GTF="$2"; shift 2 ;;
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --perl-script) PERL_SCRIPT="$2"; shift 2 ;;
        --star-bin) STAR_BIN="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --reuse-star-dir) REUSE_STAR_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --skip-star) SKIP_STAR=true; shift ;;
        --skip-perl) SKIP_PERL=true; shift ;;
        --force) FORCE=true; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

# Validate inputs
[ -f "$INPUT_FASTA" ] || error_exit "Input FASTA not found: $INPUT_FASTA"
[ -f "$INPUT_GTF" ] || error_exit "Input GTF not found: $INPUT_GTF"
[ -d "$REFERENCE_DIR" ] || error_exit "Reference dir not found: $REFERENCE_DIR"
[ -f "$REFERENCE_DIR/genome.fa" ] || error_exit "Reference FASTA not found: $REFERENCE_DIR/genome.fa"
[ -f "$REFERENCE_DIR/genome.gtf" ] || error_exit "Reference GTF not found: $REFERENCE_DIR/genome.gtf"
if [ "$SKIP_PERL" = false ]; then
    [ -f "$PERL_SCRIPT" ] || error_exit "Perl script not found: $PERL_SCRIPT"
fi

if [ "$SKIP_PERL" = true ] && [ "$SKIP_STAR" = true ]; then
    error_exit "Both Perl and STAR steps are skipped; nothing to do."
fi

mkdir -p "$WORK_DIR"

# Step 1: Perl formatter outputs vs fixtures
if [ "$SKIP_PERL" = false ]; then
    PERL_OUT_DIR="$WORK_DIR/perl"
    mkdir -p "$PERL_OUT_DIR"
    PERL_FASTA="$PERL_OUT_DIR/genome.fa"
    PERL_GTF="$PERL_OUT_DIR/genes.gtf"

    if [ "$FORCE" = true ] || [ ! -s "$PERL_FASTA" ] || [ ! -s "$PERL_GTF" ]; then
        log "Running Perl formatter..."
        INPUTFA="$INPUT_FASTA" OUTPUTFA="$PERL_FASTA" \
        INPUTGTF="$INPUT_GTF" OUTPUTGTF="$PERL_GTF" \
        overwrite=1 \
        perl "$PERL_SCRIPT" 2>&1 | grep -v "^working on" || true
    else
        log "Skipping Perl formatter (outputs exist). Use --force to re-run."
    fi

    log "Comparing Perl outputs to fixtures..."
    compare_files "Perl FASTA" "$REFERENCE_DIR/genome.fa" "$PERL_FASTA" || exit 1
    compare_files "Perl GTF" "$REFERENCE_DIR/genome.gtf" "$PERL_GTF" || exit 1
else
    log "Skipping Perl step (--skip-perl set)."
fi

# Step 2: STAR internal path vs fixtures
if [ "$SKIP_STAR" = true ]; then
    log "Skipping STAR step (--skip-star set)."
    exit 0
fi

STAR_OUT_FASTA=""
STAR_OUT_GTF=""

if [ -n "$REUSE_STAR_DIR" ]; then
    if [ -d "$REUSE_STAR_DIR/cellranger_ref" ]; then
        STAR_OUT_FASTA="$REUSE_STAR_DIR/cellranger_ref/genome.fa"
        STAR_OUT_GTF="$REUSE_STAR_DIR/cellranger_ref/genes.gtf"
    else
        STAR_OUT_FASTA="$REUSE_STAR_DIR/genome.fa"
        STAR_OUT_GTF="$REUSE_STAR_DIR/genes.gtf"
    fi

    [ -f "$STAR_OUT_FASTA" ] || error_exit "STAR output FASTA not found: $STAR_OUT_FASTA"
    [ -f "$STAR_OUT_GTF" ] || error_exit "STAR output GTF not found: $STAR_OUT_GTF"
else
    [ -x "$STAR_BIN" ] || error_exit "STAR binary not found: $STAR_BIN"
    STAR_GENOME_DIR="$WORK_DIR/star_genome"
    mkdir -p "$STAR_GENOME_DIR"

    log "Running STAR internal path (this may take a long time)..."
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$STAR_GENOME_DIR" \
        --genomeFastaFiles "$INPUT_FASTA" \
        --sjdbGTFfile "$INPUT_GTF" \
        --cellrangerStyleIndex Yes \
        --runThreadN "$THREADS"

    STAR_OUT_FASTA="$STAR_GENOME_DIR/cellranger_ref/genome.fa"
    STAR_OUT_GTF="$STAR_GENOME_DIR/cellranger_ref/genes.gtf"
    [ -f "$STAR_OUT_FASTA" ] || error_exit "STAR output FASTA not found: $STAR_OUT_FASTA"
    [ -f "$STAR_OUT_GTF" ] || error_exit "STAR output GTF not found: $STAR_OUT_GTF"
fi

log "Comparing STAR outputs to fixtures..."
compare_files "STAR FASTA" "$REFERENCE_DIR/genome.fa" "$STAR_OUT_FASTA" || exit 1
compare_files "STAR GTF" "$REFERENCE_DIR/genome.gtf" "$STAR_OUT_GTF" || exit 1

log "All comparisons passed."
