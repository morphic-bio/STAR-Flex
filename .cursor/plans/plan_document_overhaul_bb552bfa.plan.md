---
name: Plan Document Overhaul
overview: "Comprehensive update to the main plan: rename to \"Internal Spill-to-Disk Sorter\", add non-goals, promote QNAME tie-break to hard requirement, fix formatting, and restructure testing with three tiers using existing fixtures."
todos:
  - id: update-plan-naming
    content: Rename plan to Internal Spill-to-Disk Sorter, add Non-Goals section
    status: pending
  - id: update-plan-stage0
    content: Promote QNAME tie-break to hard requirement in Stage 0
    status: pending
  - id: update-plan-formatting
    content: Fix formatting issues (blank lines, Phase 4 table)
    status: pending
  - id: update-plan-testing
    content: Replace Testing Strategy with three-tier approach
    status: pending
  - id: add-code-spec-comments
    content: Add spec comments and assert guards to SamtoolsSorter.h
    status: pending
  - id: create-smoke-test
    content: Create tests/test_samtools_sorter.sh using ychrom_test fixture
    status: pending
  - id: create-unit-test
    content: Create tests/samtools_sorter_unit.cpp C++ harness (optional)
    status: pending
  - id: create-solo-fixture
    content: Create tests/solo_mini_test/ minimal fixture (optional)
    status: pending
---

# Internal Spill-to-Disk Sorter Plan Overhaul

## Summary

Comprehensive update to [`samtools_bam_sort_backend_f3c39581.plan.md`](.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md) addressing naming, requirements, formatting, and testing strategy.

---

## 1. Naming and Non-Goals

### Rename Plan

- **Old**: "Samtools BAM Sort Backend"
- **New**: "Internal Spill-to-Disk Sorter"

### Add Explicit Non-Goals Section (after Overview)

```markdown
## Non-Goals

The following are explicitly **out of scope**:

- **No htslib upgrade**: Uses existing vendored htslib BGZF APIs only
- **No samtools code vendoring**: Sorting logic implemented internally, not borrowed from samtools
- **No parity requirement with samtools sort**: Deterministic ordering defined by STAR spec, not samtools behavior
```

---

## 2. QNAME Tie-Break as Hard Requirement

### Update Stage 0 Header

```markdown
## Stage 0: Deterministic Ordering Specification (REQUIRED)

**This specification defines STAR's coordinate sorting behavior and is non-negotiable.**
```

### Add Sorting Specification Subsection (0.1)

```markdown
### 0.1 Deterministic Ordering Specification

**REQUIRED**: When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering MUST be deterministic by QNAME byte comparison.

This is a hard requirement for reproducibility and is part of STAR's coordinate sort specification going forward:

- Applies to all sources: in-memory records and spill file records
- QNAME bytes compared lexicographically (excluding trailing NUL)
- Identical input MUST produce identical output ordering
- Enforced via `assert(bamPtr != nullptr)` in heap comparator
```

Renumber existing 0.1-0.4 to 0.2-0.5.

---

## 3. Formatting Fixes

### Fix Line 54-55 (missing blank line)

```markdown
... implemented as an internal STAR component (no samtools code vendoring required).

**Key design principles:**
```

### Fix Line 66 (missing blank line after gate statement)

```markdown
**This prerequisite must be completed and verified before validation can begin.**

The internal spill-to-disk sorter must correctly implement...
```

### Fix Phase 4 Table (currently inlined, hard to read)

Convert to proper markdown table with blank lines:

```markdown
### Phase 4: Post-Sort Y/noY Split with Explicit keepBAM Semantics

**Behavior specification for samtools backend:**

| emitNoYBAM | keepBAM | Output Files |
|------------|---------|--------------|
| no | - | `Aligned.sortedByCoord.out.bam` (primary) |
| yes | no | `_Y.bam` + `_noY.bam` only (no primary) |
| yes | yes | `_Y.bam` + `_noY.bam` + `Aligned.sortedByCoord.out.bam` |

**Implementation:**
```

---

## 4. New Testing Strategy (Three Tiers)

Replace current Testing Strategy section entirely:

```markdown
## Testing Strategy

### Tier 1: Unit Tests (fast, C++ harness)

**Location**: `tests/samtools_sorter_unit.cpp` (new file)

Test `SamtoolsSorter` directly with minimal synthetic BAM records:

1. **Numeric sort order**: Records with different tid/pos sort correctly
2. **QNAME tie-break**: Records with equal keys sort by QNAME (readA < readB)
3. **Unmapped last**: `tid=-1` records sort after all mapped records
4. **Spill threshold**: Set tiny `maxRAM`, feed enough records to create >1 spill file, assert merged order matches no-spill run

### Tier 2: Functional Smoke Tests (existing fixtures)

**Primary fixture**: `tests/ychrom_test/` (self-contained: tiny genome + FASTQ)

**Wrapper script**: `tests/test_samtools_sorter.sh`

```bash
#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
GENOME_DIR="${SCRIPT_DIR}/ychrom_test/ref/star_index"
READS_R1="${SCRIPT_DIR}/ychrom_test/fastq/R1.fastq"
READS_R2="${SCRIPT_DIR}/ychrom_test/fastq/R2.fastq"
OUT_DIR="${SCRIPT_DIR}/samtools_sorter_test_output"

rm -rf "$OUT_DIR" && mkdir -p "$OUT_DIR"

# Test with spill (tiny RAM limit)
"$STAR_BIN" --runThreadN 2 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READS_R1" "$READS_R2" \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortMethod samtools \
    --limitBAMsortRAM 10000000 \
    --emitNoYBAM yes --keepBAM yes \
    --outFileNamePrefix "$OUT_DIR/test_"

# Validate
samtools quickcheck "$OUT_DIR/test_Aligned.sortedByCoord.out.bam"
samtools view -H "$OUT_DIR/test_Aligned.sortedByCoord.out.bam" | grep -q "SO:coordinate"
samtools quickcheck "$OUT_DIR/test_Aligned.sortedByCoord.out_Y.bam"
samtools quickcheck "$OUT_DIR/test_Aligned.sortedByCoord.out_noY.bam"

# Determinism check: run twice, diff output
"$STAR_BIN" ... --outFileNamePrefix "$OUT_DIR/run1_"
"$STAR_BIN" ... --outFileNamePrefix "$OUT_DIR/run2_"
diff <(samtools view "$OUT_DIR/run1_Aligned.sortedByCoord.out.bam") \
     <(samtools view "$OUT_DIR/run2_Aligned.sortedByCoord.out.bam")

echo "All tests PASSED"
```

**scRNA-seq fixture**: `tests/solo_mini_test/` (new, minimal whitelist + few reads + tiny genome)
- Avoids `/storage/*` dependencies
- Tests Solo + coordinate sorting + Y/noY split

### Tier 3: Spill Stress (targeted)

Force spills and compare to no-spill run:

```bash
# Low RAM (force spills)
./source/STAR ... --limitBAMsortRAM 10000000 --outFileNamePrefix test_spill_

# High RAM (no spills)
./source/STAR ... --limitBAMsortRAM 10000000000 --outFileNamePrefix test_normal_

# Compare: order should match exactly
diff <(samtools view test_spill_*.bam) <(samtools view test_normal_*.bam)

# Verify counts match
samtools flagstat test_spill_*.bam > spill.txt
samtools flagstat test_normal_*.bam > normal.txt
diff spill.txt normal.txt
```

### Dropped Tests

- **STAR legacy vs new backend diff**: Not required (different tie-break behavior)
- **samtools sort ordering diff**: Not a requirement (STAR defines its own spec)
```

---

## 5. Files to Modify/Create

| File | Action |
|------|--------|
| [`.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md`](.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md) | Update name, add non-goals, restructure Stage 0, fix formatting, replace testing section |
| [`source/SamtoolsSorter.h`](source/SamtoolsSorter.h) | Add spec comment header, assert guards for bamPtr |
| `tests/test_samtools_sorter.sh` | New: functional smoke test wrapper |
| `tests/samtools_sorter_unit.cpp` | New: C++ unit test harness (optional, can defer) |
| `tests/solo_mini_test/` | New: minimal Solo test fixture (optional, can defer) |

---

## 6. Implementation Order

1. Update plan document (naming, non-goals, Stage 0, formatting, testing)
2. Add spec comments and assert guards to `SamtoolsSorter.h`
3. Create `tests/test_samtools_sorter.sh` using existing `ychrom_test` fixture
4. (Optional) Create `tests/samtools_sorter_unit.cpp` C++ harness
5. (Optional) Create `tests/solo_mini_test/` fixture
