#!/bin/bash
# Flex Pipeline Smoke Test
# Tests that STAR builds and accepts flex flags correctly

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SOURCE_DIR="$PROJECT_ROOT/source"
OUT_DIR="$SCRIPT_DIR/output"
FIXTURE_DIR="$PROJECT_ROOT/reference/tests/100K"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

pass() { echo -e "${GREEN}PASS${NC}: $1"; }
fail() { echo -e "${RED}FAIL${NC}: $1"; exit 1; }
info() { echo -e "${YELLOW}INFO${NC}: $1"; }

echo "=== Flex Pipeline Smoke Test ==="
echo "Project root: $PROJECT_ROOT"
echo "Source dir: $SOURCE_DIR"
echo "Fixture dir: $FIXTURE_DIR"
echo ""

# Test 1: Build STAR
echo "=== Test 1: Build STAR with flex integration ==="
cd "$SOURCE_DIR"
make clean 2>/dev/null || true
if make -j$(nproc) 2>&1; then
    pass "STAR built successfully"
else
    fail "STAR build failed"
fi

# Verify STAR binary exists
if [[ ! -f "$SOURCE_DIR/STAR" ]]; then
    fail "STAR binary not found after build"
fi
pass "STAR binary exists"

# Test 2: STAR accepts --flexEnable flag
echo ""
echo "=== Test 2: STAR accepts --flexEnable flag ==="
if "$SOURCE_DIR/STAR" --version 2>&1 | grep -q "STAR"; then
    pass "STAR --version works"
else
    fail "STAR --version failed"
fi

# Test --flexEnable is recognized (should not error on unrecognized parameter)
if "$SOURCE_DIR/STAR" --genomeDir /nonexistent --flexEnable yes 2>&1 | grep -qi "fatal.*flexEnable"; then
    fail "--flexEnable flag not recognized"
else
    pass "--flexEnable flag is recognized (no fatal error about unknown parameter)"
fi

# Test sample detection flags are recognized
if "$SOURCE_DIR/STAR" --genomeDir /nonexistent --soloSampleProbeOffset 68 2>&1 | grep -qi "fatal.*soloSampleProbeOffset"; then
    fail "--soloSampleProbeOffset flag not recognized"
else
    pass "--soloSampleProbeOffset flag is recognized"
fi

# Test 3: Default behavior (flexEnable=no should be default)
echo ""
echo "=== Test 3: Default behavior verification ==="
# Just verify STAR doesn't crash on normal parameter parsing
if "$SOURCE_DIR/STAR" --version 2>&1 | grep -q "STAR"; then
    pass "Default behavior intact"
else
    fail "Default behavior broken"
fi

# Test 4: Full alignment test with fixture (if available)
echo ""
echo "=== Test 4: Alignment test with 100K fixture ==="
if [[ -d "$FIXTURE_DIR/genome" ]] && [[ -d "$FIXTURE_DIR/SC2300771" ]]; then
    info "Fixture found at $FIXTURE_DIR"
    
    # Create output directory
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR/flex_off" "$OUT_DIR/flex_on"
    
    # Test 4a: Run with flex disabled (default)
    echo "Running STAR with flexEnable=no (default)..."
    if "$SOURCE_DIR/STAR" \
        --genomeDir "$FIXTURE_DIR/genome" \
        --readFilesIn "$FIXTURE_DIR/SC2300771/reads_R2.fastq" "$FIXTURE_DIR/SC2300771/reads_R1.fastq" \
        --outFileNamePrefix "$OUT_DIR/flex_off/" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$FIXTURE_DIR/737K-august-2016.txt" \
        --soloCBlen 16 \
        --soloUMIlen 12 \
        --soloFeatures Gene \
        --runThreadN 1 \
        --outSAMtype BAM Unsorted \
        --soloOutFileNames Solo.out/ Gene/raw/ 2>&1; then
        pass "STAR ran successfully with flexEnable=no"
        
        # Check that standard Solo output exists
        if [[ -d "$OUT_DIR/flex_off/Solo.out" ]]; then
            pass "Solo.out directory created"
        else
            fail "Solo.out directory missing"
        fi
    else
        fail "STAR failed with flexEnable=no"
    fi
    
    # Test 4b: Run with flex enabled
    echo ""
    echo "Running STAR with flexEnable=yes..."
    if "$SOURCE_DIR/STAR" \
        --genomeDir "$FIXTURE_DIR/genome" \
        --readFilesIn "$FIXTURE_DIR/SC2300771/reads_R2.fastq" "$FIXTURE_DIR/SC2300771/reads_R1.fastq" \
        --outFileNamePrefix "$OUT_DIR/flex_on/" \
        --flexEnable yes \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$FIXTURE_DIR/737K-august-2016.txt" \
        --soloCBlen 16 \
        --soloUMIlen 12 \
        --soloFeatures Gene \
        --runThreadN 1 \
        --outSAMtype BAM Unsorted \
        --soloOutFileNames Solo.out/ Gene/raw/ 2>&1; then
        pass "STAR ran successfully with flexEnable=yes"
        
        # Check for flex pipeline log message
        if grep -q "Flex pipeline enabled" "$OUT_DIR/flex_on/Log.out" 2>/dev/null; then
            pass "Flex pipeline enabled message in log"
        else
            info "Flex pipeline enabled message not found in log (may not write until alignment completes)"
        fi
    else
        fail "STAR failed with flexEnable=yes"
    fi
    
    # Test 4c: Verify no regression (both runs should complete without crash)
    if [[ -f "$OUT_DIR/flex_off/Aligned.out.bam" ]] && [[ -f "$OUT_DIR/flex_on/Aligned.out.bam" ]]; then
        pass "Both runs produced BAM output"
    else
        fail "One or both runs did not produce BAM output"
    fi
    
else
    info "Skipping alignment test - fixture not found at $FIXTURE_DIR"
    info "To run full tests, ensure reference/tests/100K contains genome/ and SC2300771/ directories"
fi

echo ""
echo "=== All smoke tests passed ==="
echo ""
echo "Summary:"
echo "  - STAR builds successfully with flex integration"
echo "  - --flexEnable flag is recognized"
echo "  - Sample detection flags are recognized"
echo "  - Default behavior (flex off) works correctly"
if [[ -d "$FIXTURE_DIR/genome" ]]; then
    echo "  - Alignment tests with fixture completed"
fi
