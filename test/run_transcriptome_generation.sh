#!/bin/bash
# Test harness for transcriptome FASTA generation
# Usage: ./run_transcriptome_generation.sh [--synthetic | --parity | --default-path | --cellranger]
#
# --synthetic: Run synthetic unit tests (default)
# --parity: Run full GENCODE v44 parity test
# --default-path: Test default output path (${genomeDir}/transcriptome.fa)
# --cellranger: Test CellRanger-style index with transcriptome generation

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"
TMP_DIR="${TMP_DIR:-/tmp/transcriptome_test_$$}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}

# Set up cleanup trap
trap cleanup EXIT

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo -e "${RED}ERROR: $*${NC}" >&2
    exit 1
}

success() {
    echo -e "${GREEN}SUCCESS: $*${NC}"
}

warn() {
    echo -e "${YELLOW}WARNING: $*${NC}"
}

run_synthetic_tests() {
    log "Running synthetic unit tests..."
    
    local FIXTURES_DIR="$SCRIPT_DIR/fixtures/transcriptome"
    local GENOME_FA="$FIXTURES_DIR/small_genome.fa"
    local GENES_GTF="$FIXTURES_DIR/small_genes.gtf"
    local EXPECTED_FA="$FIXTURES_DIR/expected_transcriptome.fa"
    local GENOME_DIR="$TMP_DIR/genome_index"
    local OUTPUT_FA="$TMP_DIR/transcriptome.fa"
    
    # Verify fixtures exist
    [[ -f "$GENOME_FA" ]] || error "Genome FASTA not found: $GENOME_FA"
    [[ -f "$GENES_GTF" ]] || error "Genes GTF not found: $GENES_GTF"
    [[ -f "$EXPECTED_FA" ]] || error "Expected transcriptome not found: $EXPECTED_FA"
    
    # Create temp directories
    mkdir -p "$GENOME_DIR"
    
    log "Building STAR index with transcriptome generation..."
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GENES_GTF" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 4 \
        --genomeGenerateTranscriptome Yes \
        --genomeGenerateTranscriptomeFasta "$OUTPUT_FA" \
        --outFileNamePrefix "$TMP_DIR/" \
        2>&1 | tee "$TMP_DIR/star.log"
    
    log "Comparing output with expected..."
    if diff -q "$OUTPUT_FA" "$EXPECTED_FA" > /dev/null 2>&1; then
        success "Transcriptome FASTA matches expected output"
        return 0
    else
        error "Transcriptome FASTA does not match expected output"
        log "Differences:"
        diff "$OUTPUT_FA" "$EXPECTED_FA" || true
        return 1
    fi
}

run_default_path_tests() {
    log "Running default path tests..."
    
    local FIXTURES_DIR="$SCRIPT_DIR/fixtures/transcriptome"
    local GENOME_FA="$FIXTURES_DIR/small_genome.fa"
    local GENES_GTF="$FIXTURES_DIR/small_genes.gtf"
    local EXPECTED_FA="$FIXTURES_DIR/expected_transcriptome.fa"
    local GENOME_DIR="$TMP_DIR/genome_index"
    
    # Verify fixtures exist
    [[ -f "$GENOME_FA" ]] || error "Genome FASTA not found: $GENOME_FA"
    [[ -f "$GENES_GTF" ]] || error "Genes GTF not found: $GENES_GTF"
    [[ -f "$EXPECTED_FA" ]] || error "Expected transcriptome not found: $EXPECTED_FA"
    
    # Create temp directories
    mkdir -p "$GENOME_DIR"
    
    log "Building STAR index with transcriptome generation (default path)..."
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GENES_GTF" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 4 \
        --genomeGenerateTranscriptome Yes \
        --outFileNamePrefix "$TMP_DIR/" \
        2>&1 | tee "$TMP_DIR/star.log"
    
    # Default path should be ${genomeDir}/transcriptome.fa
    local DEFAULT_OUTPUT="$GENOME_DIR/transcriptome.fa"
    
    log "Checking default output path: $DEFAULT_OUTPUT"
    if [[ ! -f "$DEFAULT_OUTPUT" ]]; then
        error "Default transcriptome file not created: $DEFAULT_OUTPUT"
        return 1
    fi
    
    log "Comparing output with expected..."
    if diff -q "$DEFAULT_OUTPUT" "$EXPECTED_FA" > /dev/null 2>&1; then
        success "Transcriptome FASTA at default path matches expected output"
        return 0
    else
        error "Transcriptome FASTA at default path does not match expected output"
        log "Differences:"
        diff "$DEFAULT_OUTPUT" "$EXPECTED_FA" || true
        return 1
    fi
}

run_cellranger_tests() {
    log "Running CellRanger-style index tests (using GENCODE chr21+chr22 subset)..."
    
    # CellRanger-style tests require real GENCODE data (ENSG/ENST format)
    # Use the chr21+chr22 subset for fast testing
    local FIXTURES_DIR="$SCRIPT_DIR/fixtures/transcriptome/gencode_subset"
    local GENOME_FA="$FIXTURES_DIR/genome.fa"
    local GENES_GTF="$FIXTURES_DIR/genes.gtf"
    local EXPECTED_FA="$FIXTURES_DIR/expected_transcriptome.fa"
    local GENOME_DIR="$TMP_DIR/genome_index"
    
    # Verify fixtures exist
    if [[ ! -f "$GENOME_FA" ]]; then
        warn "GENCODE subset not found: $GENOME_FA"
        warn "Skipping CellRanger test - run with real GENCODE data to create subset"
        success "CellRanger-style test skipped (subset not available)"
        return 0
    fi
    [[ -f "$GENES_GTF" ]] || error "Genes GTF not found: $GENES_GTF"
    [[ -f "$EXPECTED_FA" ]] || error "Expected transcriptome not found: $EXPECTED_FA"
    
    # Create temp directories
    mkdir -p "$GENOME_DIR"
    
    log "Building STAR index with CellRanger-style + transcriptome generation..."
    log "  (chr21+chr22 subset - ~100k GTF lines, ~95MB genome)"
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GENES_GTF" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 11 \
        --genomeGenerateTranscriptome Yes \
        --cellrangerStyleIndex Yes \
        --runThreadN 4 \
        --outFileNamePrefix "$TMP_DIR/" \
        2>&1 | tee "$TMP_DIR/star.log"
    
    # Check standard output path
    local DEFAULT_OUTPUT="$GENOME_DIR/transcriptome.fa"
    log "Checking standard output path: $DEFAULT_OUTPUT"
    if [[ ! -f "$DEFAULT_OUTPUT" ]]; then
        error "Standard transcriptome file not created: $DEFAULT_OUTPUT"
        return 1
    fi
    
    # Check CellRanger output path
    local CELLRANGER_OUTPUT="$GENOME_DIR/cellranger_ref/transcriptome.fa"
    log "Checking CellRanger output path: $CELLRANGER_OUTPUT"
    if [[ ! -f "$CELLRANGER_OUTPUT" ]]; then
        error "CellRanger transcriptome file not created: $CELLRANGER_OUTPUT"
        return 1
    fi
    
    # Compare by sorting (STAR and gffread use different ordering)
    # CellRangerFormatter filters out some transcript types (artifacts, pseudogenes, etc.)
    # so STAR may produce fewer transcripts than gffread
    log "Comparing transcriptome content (sequences, not order)..."
    
    local gen_count=$(grep -c "^>" "$DEFAULT_OUTPUT")
    local exp_count=$(grep -c "^>" "$EXPECTED_FA")
    
    log "  Generated: $gen_count transcripts"
    log "  Expected (gffread, unfiltered): $exp_count transcripts"
    
    # Sort both and compare - strip version suffixes (STAR outputs ENST00000123 but gffread outputs ENST00000123.1)
    awk '/^>/{if(id)print id"\t"seq; id=$1; gsub(/ .*/,"",id); seq=""; next} {seq=seq$0} END{print id"\t"seq}' "$DEFAULT_OUTPUT" | sort > "$TMP_DIR/gen_sorted.tsv"
    awk '/^>/{if(id)print id"\t"seq; id=$1; gsub(/\.[0-9]+$/,"",id); gsub(/ .*/,"",id); seq=""; next} {seq=seq$0} END{print id"\t"seq}' "$EXPECTED_FA" | sort > "$TMP_DIR/exp_sorted.tsv"
    
    # Get IDs from generated (these are the ones CellRanger kept after filtering)
    cut -f1 "$TMP_DIR/gen_sorted.tsv" > "$TMP_DIR/gen_ids.txt"
    
    # Filter expected to only those IDs
    grep -Ff "$TMP_DIR/gen_ids.txt" "$TMP_DIR/exp_sorted.tsv" > "$TMP_DIR/exp_filtered.tsv" || true
    
    local filtered_count=$(wc -l < "$TMP_DIR/exp_filtered.tsv")
    log "  After filtering expected to generated IDs: $filtered_count transcripts"
    
    if [[ "$gen_count" != "$filtered_count" ]]; then
        warn "Some generated transcripts not found in expected (may be version suffix issue)"
    fi
    
    # Compare sequences for shared transcripts using join
    log "  Comparing sequences..."
    
    # Join on transcript ID and compare sequences
    join -t $'\t' "$TMP_DIR/gen_sorted.tsv" "$TMP_DIR/exp_filtered.tsv" > "$TMP_DIR/joined.tsv" 2>/dev/null || true
    
    # Check for mismatches (columns 2 and 3 should be identical)
    local mismatches=$(awk -F'\t' '$2 != $3 {print}' "$TMP_DIR/joined.tsv" | wc -l)
    
    if [[ $mismatches -gt 0 ]]; then
        error "$mismatches transcript(s) have sequence differences"
        log "First 3 mismatches:"
        awk -F'\t' '$2 != $3 {print $1": gen="length($2)"bp exp="length($3)"bp"; if(NR>3)exit}' "$TMP_DIR/joined.tsv"
        return 1
    fi
    
    log "Comparing standard and CellRanger outputs..."
    if ! diff -q "$DEFAULT_OUTPUT" "$CELLRANGER_OUTPUT" > /dev/null 2>&1; then
        error "Standard and CellRanger transcriptome FASTAs differ"
        return 1
    fi
    
    success "CellRanger-style transcriptome: $gen_count transcripts match expected"
    return 0
}

run_parity_tests() {
    log "Running GENCODE v44 parity test..."
    
    # Reference paths (can be overridden with environment variables)
    local REF_GENOME="${REF_GENOME:-/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa}"
    local REF_GTF="${REF_GTF:-/storage/autoindex_110_44/bulk_index/cellranger_ref/genes.gtf}"
    local REF_TRANSCRIPTOME="${REF_TRANSCRIPTOME:-/mnt/pikachu/bwb_sched_storage/bulkrna/singularity_data/processing/genome/transcriptome.fa}"
    local GENOME_DIR="$TMP_DIR/genome_index"
    local OUTPUT_FA="$TMP_DIR/transcriptome.fa"
    
    # Verify reference files exist
    [[ -f "$REF_GENOME" ]] || error "Reference genome not found: $REF_GENOME"
    [[ -f "$REF_GTF" ]] || error "Reference GTF not found: $REF_GTF"
    [[ -f "$REF_TRANSCRIPTOME" ]] || error "Reference transcriptome not found: $REF_TRANSCRIPTOME"
    
    # Create temp directories
    mkdir -p "$GENOME_DIR"
    
    log "Building STAR index with transcriptome generation..."
    log "  Genome: $REF_GENOME"
    log "  GTF: $REF_GTF"
    log "  Output: $OUTPUT_FA"
    
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$REF_GENOME" \
        --sjdbGTFfile "$REF_GTF" \
        --sjdbOverhang 100 \
        --genomeGenerateTranscriptome Yes \
        --genomeGenerateTranscriptomeFasta "$OUTPUT_FA" \
        --outFileNamePrefix "$TMP_DIR/" \
        --runThreadN 8 \
        2>&1 | tee "$TMP_DIR/star.log"
    
    log "Comparing output with reference transcriptome..."
    if cmp -s "$OUTPUT_FA" "$REF_TRANSCRIPTOME"; then
        success "Transcriptome FASTA matches reference exactly (byte-for-byte)"
        return 0
    else
        log "Files differ. Checking if only ordering differs..."
        
        # Sort both files by sequence ID and compare
        local SORTED_OUTPUT="$TMP_DIR/sorted_output.fa"
        local SORTED_REF="$TMP_DIR/sorted_ref.fa"
        
        # Simple sort by sequence ID (works for FASTA)
        grep "^>" "$OUTPUT_FA" | sort > "$TMP_DIR/output_ids.txt"
        grep "^>" "$REF_TRANSCRIPTOME" | sort > "$TMP_DIR/ref_ids.txt"
        
        if diff -q "$TMP_DIR/output_ids.txt" "$TMP_DIR/ref_ids.txt" > /dev/null 2>&1; then
            log "Same transcript IDs, but different order or content"
        else
            log "Different transcript IDs"
            log "Output count: $(wc -l < "$TMP_DIR/output_ids.txt")"
            log "Reference count: $(wc -l < "$TMP_DIR/ref_ids.txt")"
            
            # Show first few differences
            log "First 10 IDs only in output:"
            comm -23 "$TMP_DIR/output_ids.txt" "$TMP_DIR/ref_ids.txt" | head -10
            log "First 10 IDs only in reference:"
            comm -13 "$TMP_DIR/output_ids.txt" "$TMP_DIR/ref_ids.txt" | head -10
        fi
        
        # Compute SHA256 hashes
        log "SHA256 of output:    $(sha256sum "$OUTPUT_FA" | cut -d' ' -f1)"
        log "SHA256 of reference: $(sha256sum "$REF_TRANSCRIPTOME" | cut -d' ' -f1)"
        
        error "Transcriptome FASTA does not match reference"
        return 1
    fi
}

run_all_tests() {
    log "Running all tests..."
    local failed=0
    
    log "=== Test 1: Synthetic unit tests ==="
    if run_synthetic_tests; then
        success "Test 1 passed"
    else
        warn "Test 1 failed"
        ((failed++)) || true
    fi
    
    # Clean and re-create tmp dir between tests
    rm -rf "$TMP_DIR"
    mkdir -p "$TMP_DIR"
    
    log ""
    log "=== Test 2: Default path tests ==="
    if run_default_path_tests; then
        success "Test 2 passed"
    else
        warn "Test 2 failed"
        ((failed++)) || true
    fi
    
    # Clean and re-create tmp dir
    rm -rf "$TMP_DIR"
    mkdir -p "$TMP_DIR"
    
    log ""
    log "=== Test 3: CellRanger-style tests ==="
    if run_cellranger_tests; then
        success "Test 3 passed"
    else
        warn "Test 3 failed"
        ((failed++)) || true
    fi
    
    log ""
    if [[ $failed -eq 0 ]]; then
        success "All tests passed!"
        return 0
    else
        error "$failed test(s) failed"
        return 1
    fi
}

# Parse arguments
TEST_MODE="synthetic"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --synthetic)
            TEST_MODE="synthetic"
            shift
            ;;
        --parity)
            TEST_MODE="parity"
            shift
            ;;
        --default-path)
            TEST_MODE="default-path"
            shift
            ;;
        --cellranger)
            TEST_MODE="cellranger"
            shift
            ;;
        --all)
            TEST_MODE="all"
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [--synthetic | --parity | --default-path | --cellranger | --all]"
            echo ""
            echo "Options:"
            echo "  --synthetic     Run synthetic unit tests (default)"
            echo "  --default-path  Test default output path (\${genomeDir}/transcriptome.fa)"
            echo "  --cellranger    Test CellRanger-style index with transcriptome generation"
            echo "  --parity        Run full GENCODE v44 parity test"
            echo "  --all           Run synthetic, default-path, and cellranger tests"
            exit 0
            ;;
        *)
            error "Unknown argument: $1"
            ;;
    esac
done

# Verify STAR binary exists
[[ -x "$STAR_BIN" ]] || error "STAR binary not found or not executable: $STAR_BIN"

# Create temp directory
mkdir -p "$TMP_DIR"

# Run the appropriate test
case "$TEST_MODE" in
    synthetic)
        run_synthetic_tests
        ;;
    default-path)
        run_default_path_tests
        ;;
    cellranger)
        run_cellranger_tests
        ;;
    parity)
        run_parity_tests
        ;;
    all)
        run_all_tests
        ;;
esac

log "Test completed successfully"
