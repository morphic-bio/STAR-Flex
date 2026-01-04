#!/bin/bash
# Count C/C++ lines and files in STAR-Flex vs upstream STAR

set -e

REPO_DIR="/mnt/pikachu/STAR-Flex"
UPSTREAM_URL="https://github.com/alexdobin/STAR"

# Category heuristics (edit patterns to tune classification).
FLEX_PATTERNS=(
    "*/libflex/*"
    "*/Flex*"
    "*/Flex*/*"
    "*/ProbeList*"
    "*/FlexProbe*"
)
BULK_PATTERNS=(
    "*/libem/*"
    "*/TranscriptQuant*"
    "*/Quantifications*"
    "*/ReadAlign_quantTranscriptome*"
    "*/Transcriptome_quantAlign*"
    "*/Transcriptome_geneCounts*"
    "*/TranscriptomeFasta*"
)

echo "=========================================="
echo "C/C++ Code Statistics: STAR-Flex vs Upstream"
echo "=========================================="
echo ""

# Collect source files (excluding third-party)
collect_source_files() {
    local dir="$1"
    find "$dir/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
        ! -path "*/htslib/*" \
        ! -path "*/opal/*" \
        ! -path "*third_party*" \
        -print0
}

# Match file against any glob pattern in the list
matches_any() {
    local file="$1"
    shift
    local pattern
    for pattern in "$@"; do
        if [[ "$file" == $pattern ]]; then
            return 0
        fi
    done
    return 1
}

# Count C/C++ lines and files by label
# Function to count lines and files
count_code() {
    local dir="$1"
    local label="$2"
    
    # Find all .c, .cpp, .h, .hpp files, excluding htslib and other third-party
    local files=()
    mapfile -d '' -t files < <(collect_source_files "$dir")
    local file_count=${#files[@]}
    local line_count=0
    
    if [[ "$file_count" -gt 0 ]]; then
        line_count=$(printf '%s\0' "${files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    fi
    
    echo "$label:"
    echo "  Files: $file_count"
    echo "  Lines: $line_count"
    echo ""
    
    # Return values for later use
    eval "${label//[^a-zA-Z]/_}_files=$file_count"
    eval "${label//[^a-zA-Z]/_}_lines=$line_count"
}

# Count current repo
echo "--- Current STAR-Flex Repository ---"
count_code "$REPO_DIR" "STAR_Flex"

echo "--- STAR-Flex Category Breakdown ---"
flex_files=()
bulk_files=()
shared_files=()
all_files=()
mapfile -d '' -t all_files < <(collect_source_files "$REPO_DIR")

for f in "${all_files[@]}"; do
    is_flex=0
    is_bulk=0
    if matches_any "$f" "${FLEX_PATTERNS[@]}"; then
        is_flex=1
    fi
    if matches_any "$f" "${BULK_PATTERNS[@]}"; then
        is_bulk=1
    fi

    if [[ "$is_flex" -eq 1 && "$is_bulk" -eq 0 ]]; then
        flex_files+=("$f")
    elif [[ "$is_bulk" -eq 1 && "$is_flex" -eq 0 ]]; then
        bulk_files+=("$f")
    else
        shared_files+=("$f")
    fi
done

count_category() {
    local label="$1"
    shift
    local files=("$@")
    local file_count=${#files[@]}
    local line_count=0
    if [[ "$file_count" -gt 0 ]]; then
        line_count=$(printf '%s\0' "${files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    fi
    echo "  $label:"
    echo "    Files: $file_count"
    echo "    Lines: $line_count"
}

count_category "Flex-specific" "${flex_files[@]}"
count_category "Bulk/PE-specific" "${bulk_files[@]}"
count_category "Shared/core (includes overlap)" "${shared_files[@]}"
echo ""

# Clone upstream to temp dir for comparison
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "--- Fetching Upstream STAR (shallow clone) ---"
git clone --depth 1 --quiet "$UPSTREAM_URL" "$TEMP_DIR/STAR" 2>/dev/null || {
    echo "Failed to clone upstream. Using local comparison only."
    exit 1
}
echo ""

echo "--- Upstream STAR Repository ---"
count_code "$TEMP_DIR/STAR" "Upstream_STAR"

# Calculate differences
echo "=========================================="
echo "Comparison"
echo "=========================================="

# Get the counts again for comparison
flex_files=$(find "$REPO_DIR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null | wc -l)
flex_lines=$(find "$REPO_DIR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null | xargs wc -l | tail -1 | awk '{print $1}')

upstream_files=$(find "$TEMP_DIR/STAR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null | wc -l)
upstream_lines=$(find "$TEMP_DIR/STAR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null | xargs wc -l | tail -1 | awk '{print $1}')

file_diff=$((flex_files - upstream_files))
line_diff=$((flex_lines - upstream_lines))

echo ""
echo "STAR-Flex:     $flex_files files,  $flex_lines lines"
echo "Upstream STAR: $upstream_files files,  $upstream_lines lines"
echo ""
echo "Difference:    +$file_diff files, +$line_diff lines"
echo ""

# Show new files in STAR-Flex
echo "=========================================="
echo "New/Modified Files in STAR-Flex (source/)"
echo "=========================================="
echo ""

# List files that are new or significantly different
echo "New C++ files added:"
for f in $(find "$REPO_DIR/source" -maxdepth 1 -type f \( -name "*.cpp" -o -name "*.h" \) -printf "%f\n" | sort); do
    if [[ ! -f "$TEMP_DIR/STAR/source/$f" ]]; then
        lines=$(wc -l < "$REPO_DIR/source/$f")
        echo "  + $f ($lines lines)"
    fi
done

echo ""
echo "Files in libflex/:"
if [[ -d "$REPO_DIR/source/libflex" ]]; then
    libflex_files=$(find "$REPO_DIR/source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libflex_lines=$(find "$REPO_DIR/source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libflex_files files, $libflex_lines lines"
fi

echo ""
echo "Files in solo/:"
if [[ -d "$REPO_DIR/source/solo" ]]; then
    solo_files=$(find "$REPO_DIR/source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    solo_lines=$(find "$REPO_DIR/source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $solo_files files, $solo_lines lines"
fi

echo ""
echo "Done."
