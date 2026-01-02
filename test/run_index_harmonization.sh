#!/bin/bash
#
# run_index_harmonization.sh
# Master orchestration script for Index Harmonization Runbook
# Implements all phases from plans/index_harmonization_plan.md
#
# Usage:
#   ./test/run_index_harmonization.sh [--skip-phase N] [--phase N] [--help]
#
# Options:
#   --skip-phase N    Skip a specific phase (1-6)
#   --phase N         Run only a specific phase (1-6)
#   --no-build        Skip building binaries (assume already built)
#   --help            Show this help message
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Default paths from runbook
DEFAULT_INPUT_FASTA="/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DEFAULT_INPUT_GTF="/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/gencode.v32.primary_assembly.annotation.gtf"
DEFAULT_REFERENCE_DIR="/mnt/pikachu/bwb_sched_storage/bulkrna/singularity_data/processing/genome"

# Phase control
SKIP_PHASES=()
RUN_ONLY_PHASE=""
NO_BUILD=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Master orchestration script for Index Harmonization Runbook.
Implements all phases from plans/index_harmonization_plan.md

Options:
  --skip-phase N     Skip a specific phase (1-6)
  --phase N         Run only a specific phase (1-6)
  --no-build        Skip building binaries (assume already built)
  -h, --help        Show this help message

Phases:
  Phase 0: Preflight Checklist
  Phase 1: Build Harness Binaries
  Phase 2: Synthetic Unit Tests (C++ formatter)
  Phase 3: Full Ensembl/GENCODE v32 → CellRanger Parity
  Phase 4: Optional Flex Probe Artifact Parity
  Phase 5: Triage if Diffs Remain (documentation only)
  Phase 6: Validate Inside STAR (documentation only)

Example:
  $0                    # Run all phases
  $0 --phase 2          # Run only Phase 2
  $0 --skip-phase 4    # Run all phases except Phase 4
EOF
    exit "${1:-0}"
}

log() {
    echo -e "${BLUE}[index_harmonization]${NC} $*" >&2
}

log_success() {
    echo -e "${GREEN}[index_harmonization]${NC} ✓ $*" >&2
}

log_error() {
    echo -e "${RED}[index_harmonization]${NC} ✗ $*" >&2
}

log_warning() {
    echo -e "${YELLOW}[index_harmonization]${NC} ⚠ $*" >&2
}

should_run_phase() {
    local phase="$1"
    
    # If --phase N specified, only run that phase
    if [ -n "$RUN_ONLY_PHASE" ]; then
        [ "$phase" = "$RUN_ONLY_PHASE" ]
        return
    fi
    
    # Check if phase is in skip list
    for skip in "${SKIP_PHASES[@]}"; do
        if [ "$phase" = "$skip" ]; then
            return 1
        fi
    done
    return 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-phase)
            SKIP_PHASES+=("$2")
            shift 2
            ;;
        --phase)
            RUN_ONLY_PHASE="$2"
            shift 2
            ;;
        --no-build)
            NO_BUILD=true
            shift
            ;;
        -h|--help)
            usage 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage 1
            ;;
    esac
done

# Phase 0: Preflight Checklist
run_phase_0() {
    log "=============================================="
    log "Phase 0: Preflight Checklist"
    log "=============================================="
    
    local errors=0
    
    # Check raw input files
    log "Checking raw input files..."
    if [ ! -f "$DEFAULT_INPUT_FASTA" ]; then
        log_error "Input FASTA not found: $DEFAULT_INPUT_FASTA"
        errors=$((errors + 1))
    else
        log_success "Input FASTA found"
    fi
    
    if [ ! -f "$DEFAULT_INPUT_GTF" ]; then
        log_error "Input GTF not found: $DEFAULT_INPUT_GTF"
        errors=$((errors + 1))
    else
        log_success "Input GTF found"
    fi
    
    # Check reference outputs
    log "Checking reference outputs..."
    if [ ! -d "$DEFAULT_REFERENCE_DIR" ]; then
        log_error "Reference directory not found: $DEFAULT_REFERENCE_DIR"
        errors=$((errors + 1))
    else
        if [ ! -f "$DEFAULT_REFERENCE_DIR/genome.fa" ]; then
            log_error "Reference FASTA not found: $DEFAULT_REFERENCE_DIR/genome.fa"
            errors=$((errors + 1))
        else
            log_success "Reference FASTA found"
        fi
        
        if [ ! -f "$DEFAULT_REFERENCE_DIR/genome.gtf" ]; then
            log_error "Reference GTF not found: $DEFAULT_REFERENCE_DIR/genome.gtf"
            errors=$((errors + 1))
        else
            log_success "Reference GTF found"
        fi
    fi
    
    # Check toolchain
    log "Checking toolchain..."
    if ! command -v g++ >/dev/null 2>&1 && ! command -v clang++ >/dev/null 2>&1; then
        log_error "C++ compiler not found"
        errors=$((errors + 1))
    else
        log_success "C++ compiler found"
    fi
    
    if [ ! -d "$REPO_ROOT/source/htslib" ]; then
        log_error "htslib not found: $REPO_ROOT/source/htslib"
        errors=$((errors + 1))
    else
        log_success "htslib found"
    fi
    
    # Check Perl script (optional)
    PERL_SCRIPT="/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl"
    if [ ! -f "$PERL_SCRIPT" ]; then
        log_warning "Perl script not found (optional): $PERL_SCRIPT"
        log_warning "  Phase 2.2 (Perl-vs-C++ parity) will be skipped"
    else
        log_success "Perl script found"
    fi
    
    if [ $errors -gt 0 ]; then
        log_error "Preflight checks failed with $errors error(s)"
        return 1
    fi
    
    log_success "All preflight checks passed"
    return 0
}

# Phase 1: Build Harness Binaries
run_phase_1() {
    log "=============================================="
    log "Phase 1: Build Harness Binaries"
    log "=============================================="
    
    if [ "$NO_BUILD" = true ]; then
        log "Skipping build (--no-build specified)"
        return 0
    fi
    
    cd "$REPO_ROOT/source"
    
    log "Building test_cellranger_format..."
    if make test_cellranger_format; then
        log_success "test_cellranger_format built"
    else
        log_error "Failed to build test_cellranger_format"
        return 1
    fi
    
    log "Building test_flex_probe_index..."
    if make test_flex_probe_index; then
        log_success "test_flex_probe_index built"
    else
        log_warning "Failed to build test_flex_probe_index (optional)"
    fi
    
    return 0
}

# Phase 2: Synthetic Unit Tests
run_phase_2() {
    log "=============================================="
    log "Phase 2: Synthetic Unit Tests (C++ formatter)"
    log "=============================================="
    
    local failures=0
    
    # Phase 2.1: Fixture Golden Tests
    log ""
    log "Phase 2.1: Fixture Golden Tests"
    cd "$REPO_ROOT/source"
    if make test_cellranger_golden; then
        log_success "Fixture golden tests passed"
    else
        log_error "Fixture golden tests failed"
        failures=$((failures + 1))
    fi
    
    # Phase 2.2: Perl-vs-C++ Parity on Fixtures
    log ""
    log "Phase 2.2: Perl-vs-C++ Parity on Fixtures"
    PERL_SCRIPT="/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl"
    if [ -f "$PERL_SCRIPT" ]; then
        cd "$REPO_ROOT"
        if ./test/run_cellranger_format_parity.sh; then
            log_success "Perl-vs-C++ parity test passed"
        else
            log_error "Perl-vs-C++ parity test failed"
            failures=$((failures + 1))
        fi
    else
        log_warning "Skipping Perl-vs-C++ parity (Perl script not found)"
    fi
    
    # Phase 2.3: Additional Synthetic Edge Cases
    log ""
    log "Phase 2.3: Additional Synthetic Edge Cases"
    if [ -f "$REPO_ROOT/test/run_cellranger_format_synthetic.sh" ]; then
        PERL_SCRIPT="${PERL_SCRIPT:-/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl}"
        if [ -f "$PERL_SCRIPT" ]; then
            cd "$REPO_ROOT"
            if ./test/run_cellranger_format_synthetic.sh; then
                log_success "Synthetic edge case tests passed"
            else
                log_error "Synthetic edge case tests failed"
                failures=$((failures + 1))
            fi
        else
            log_warning "Skipping synthetic edge cases (Perl script not found)"
        fi
    else
        log_warning "Synthetic test script not found: test/run_cellranger_format_synthetic.sh"
    fi
    
    if [ $failures -gt 0 ]; then
        log_error "Phase 2 completed with $failures failure(s)"
        return 1
    fi
    
    log_success "All Phase 2 tests passed"
    return 0
}

# Phase 3: Full Ensembl/GENCODE v32 → CellRanger Parity
run_phase_3() {
    log "=============================================="
    log "Phase 3: Perl → Fixtures → STAR Internal Parity"
    log "=============================================="
    
    cd "$REPO_ROOT"
    
    log "Running Perl → fixtures → STAR internal path parity..."
    log "  This may take a long time for full references..."
    
    if ./test/run_cellranger_perl_star_parity.sh --force; then
        log_success "Full reference parity test passed"
        return 0
    else
        log_error "Full reference parity test failed"
        log "  Check the output above for details on differences"
        log "  See Phase 5 in the runbook for triage guidance"
        return 1
    fi
}

# Phase 4: Optional Flex Probe Artifact Parity
run_phase_4() {
    log "=============================================="
    log "Phase 4: Optional Flex Probe Artifact Parity"
    log "=============================================="
    
    log_warning "Phase 4 requires probe CSV and reference flex index"
    log_warning "Skipping Phase 4 (not implemented in this script)"
    log_warning "To run manually:"
    log_warning "  ./test/run_reference_parity.sh \\"
    log_warning "    --force \\"
    log_warning "    --probe-csv /path/to/probes.csv \\"
    log_warning "    --flex-reference-dir /path/to/flex_index_or_flex_probe_artifacts"
    
    return 0
}

# Phase 5: Triage if Diffs Remain (documentation only)
run_phase_5() {
    log "=============================================="
    log "Phase 5: Triage if Diffs Remain"
    log "=============================================="
    
    log "Phase 5 is documentation/guidance only."
    log "See plans/index_harmonization_plan.md for triage steps."
    log ""
    log "Quick checks to localize differences:"
    log "  - Contig list diff (FASTA)"
    log "  - Gene list diff (GTF): count + missing/extra gene_ids"
    log "  - Look for mismatched header lines or GTF attribute ordering"
    
    return 0
}

# Phase 6: Validate Inside STAR (documentation only)
run_phase_6() {
    log "=============================================="
    log "Phase 6: Validate Inside STAR"
    log "=============================================="
    
    log "Phase 6 is documentation/guidance only."
    log "Once external parity is confirmed, run:"
    log "  STAR --runMode genomeGenerate --cellrangerStyleIndex Yes ..."
    log "  Verify outputs under genomeDir/cellranger_ref/ match reference"
    
    return 0
}

# Main execution
main() {
    log "=============================================="
    log "Index Harmonization Runbook Implementation"
    log "=============================================="
    log ""
    
    local phase_failures=0
    
    # Always run Phase 0 first
    if ! run_phase_0; then
        log_error "Preflight checks failed. Aborting."
        exit 1
    fi
    
    # Run phases based on configuration
    for phase in 1 2 3 4 5 6; do
        if should_run_phase "$phase"; then
            if ! "run_phase_$phase"; then
                phase_failures=$((phase_failures + 1))
                log_error "Phase $phase failed"
            fi
        else
            log "Skipping Phase $phase"
        fi
    done
    
    # Summary
    log ""
    log "=============================================="
    log "Summary"
    log "=============================================="
    
    if [ $phase_failures -eq 0 ]; then
        log_success "All phases completed successfully"
        exit 0
    else
        log_error "$phase_failures phase(s) failed"
        log "Review the output above for details"
        exit 1
    fi
}

main "$@"
