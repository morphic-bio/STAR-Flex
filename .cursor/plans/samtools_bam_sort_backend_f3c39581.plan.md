---
name: Internal Spill-to-Disk Sorter
overview: Implement an optional internal spill-to-disk coordinate-sorting backend, selectable via a new CLI parameter while preserving the legacy STAR bin sorter as the default.
todos:
  - id: add-parameter
    content: Add outBAMsortMethod parameter to Parameters.h/.cpp and parametersDefault
    status: completed
  - id: branching-logic
    content: Add backend branching at coordOneAlign call sites
    status: completed
    dependencies:
      - add-parameter
  - id: samtools-backend
    content: Implement bamSortSamtools with spill-to-disk and QNAME tie-breaker
    status: completed
    dependencies:
      - branching-logic
  - id: y-noy-split
    content: Implement post-sort Y/noY stream splitting with explicit keepBAM behavior
    status: completed
    dependencies:
      - samtools-backend
  - id: build-system
    content: Update Makefile with new sources and verify compilation
    status: completed
    dependencies:
      - samtools-backend
  - id: validation-basic
    content: Validate basic sorting, headers, and unmapped read ordering
    status: pending
    dependencies:
      - build-system
  - id: validation-split
    content: Validate Y/noY split and keepBAM behavior
    status: pending
    dependencies:
      - validation-basic
  - id: validation-stress
    content: Stress test spill-to-disk with large dataset
    status: pending
    dependencies:
      - validation-basic
  - id: validation-pipelines
    content: Test both bulk RNA-seq and scRNA-seq (Solo) pipelines
    status: pending
    dependencies:
      - validation-split
---

# Internal Spill-to-Disk Sorter

## Overview

This plan adds an internal spill-to-disk coordinate-sorting backend to STAR-Flex, selectable via `--outBAMsortMethod samtools`. The legacy STAR bin sorter remains the default. The new backend uses **raw BAM buffer sorting** with bounded-memory spill-to-disk and k-way merge, implemented as an internal STAR component (no samtools code vendoring required).

**Key design principles:**

- **Raw buffer sorting**: Parse BAM core fields directly from raw buffers without converting to `bam1_t` structures
- **Memory-bounded**: Check memory usage and spill to disk when `--limitBAMsortRAM` is exceeded
- **QNAME tie-breaking**: Deterministic ordering via byte-level QNAME comparison when numeric keys are equal
- **No BAM mutation**: Sorting metadata stored separately, never modifies BAM record content
- **Existing BGZF I/O**: Reuses STAR's existing `bgzf_write()` and `outBAMwriteHeader()` paths

## Non-Goals

The following are explicitly **out of scope**:

- **No htslib upgrade**: Uses existing vendored htslib BGZF APIs only
- **No samtools code vendoring**: Sorting logic implemented internally, not borrowed from samtools
- **No samtools ordering parity**: Deterministic ordering defined by STAR spec, not samtools behavior

---

## Stage 0: Deterministic Ordering Specification (REQUIRED)

**This specification defines STAR's coordinate sorting behavior and is non-negotiable.**

The internal spill-to-disk sorter must correctly implement deterministic ordering and handle memory pressure.

### 0.1 Deterministic Ordering Specification

**REQUIRED**: When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering MUST be deterministic by QNAME byte comparison.

This is a hard requirement for reproducibility and is part of STAR's coordinate sort specification going forward:

- Applies to all sources: in-memory records and spill file records
- QNAME bytes compared lexicographically (excluding trailing NUL)
- Identical input MUST produce identical output ordering
- Enforced via `assert(bamPtr != nullptr)` in heap comparator

**QNAME tie-breaking is part of STAR's sort spec (non-optional).**

### 0.2 Core Field Extraction

The sorter parses BAM core fields directly from raw buffers:

- `bam32[1]=tid`, `bam32[2]=pos`
- `bam32[4]=flag_nc` (flag extracted from upper 16 bits)
- `bam32[6]=mtid`, `bam32[7]=mpos`, `bam32[8]=isize`

Unmapped reads (`tid=-1` or `pos=-1`) are normalized to `INT32_MAX` to sort last.

### 0.3 SortKey Definition

```cpp
struct SortKey {
    int32_t tid;        // refID (-1 for unmapped -> INT32_MAX)
    int32_t pos;        // position (-1 for unmapped -> INT32_MAX)
    uint16_t flag;      // BAM flag
    int32_t mtid;       // mate refID
    int32_t mpos;       // mate position
    int32_t isize;      // insert size
    // QNAME comparison done separately for tie-breaking
};
```

Lexicographic comparison: `tid < pos < flag < mtid < mpos < isize`, then QNAME bytes.

### 0.4 QNAME Tie-Breaking

When all numeric SortKey fields are equal, QNAME bytes are compared directly from raw BAM buffers:

- Extract `l_qname` from core header (offset 8)
- Compare QNAME bytes (excluding trailing NUL)
- Ensures deterministic ordering without mutating BAM content
- **Required**: This tie-breaking is mandatory, not optional

### 0.5 Exit Criteria (must all pass)

- [ ] Compilation succeeds with no errors or warnings
- [ ] Core field offsets match BAM format specification
- [ ] Heap comparator correctly handles QNAME tie-breaking
- [ ] Memory accounting triggers spill-to-disk when threshold exceeded
- [ ] K-way merge produces correctly ordered output

**Gate status**: Validation can begin once exit criteria are met.

---

## Architecture

### Current STAR Bin Sorter (unchanged)

```mermaid
flowchart LR
    A[ReadAlignChunk threads] --> B["coordOneAlign()"]
    B --> C[bin temp files]
    C --> D[BAMbinSortByCoordinate]
    D --> E[bam_cat]
    E --> F[Output BAM]
```



### New Internal Spill-to-Disk Sorter

```mermaid
flowchart LR
    A2[ReadAlignChunk threads] --> S["SamtoolsSorter::addRecord()"]
    S --> |in-memory buffer| G[SortKey computation]
    G --> |memory check| T{Exceed RAM?}
    T -->|yes| SP[Sort & spill to disk]
    T -->|no| BUF[Buffer records]
    SP --> T2[(Temp spill files)]
    BUF --> FIN[finalize()]
    T2 --> FIN
    FIN --> M[k-way merge with QNAME tie-break]
    M --> P{Post-sort split}
    P --> |Y reads| Y[_Y.bam]
    P --> |noY reads| N[_noY.bam]
    P --> |keepBAM=yes| K[primary.bam]
```

**Key differences from STAR bin sorter:**

- **Direct ingestion**: Records streamed directly to `SamtoolsSorter` during alignment, bypassing bin temp files
- **Memory-bounded**: Checks `currentRAM_` against `--limitBAMsortRAM` and spills sorted chunks to disk
- **K-way merge**: Merges multiple sorted spill files + in-memory buffer using min-heap with QNAME tie-breaking
- **Raw buffer sorting**: Parses BAM core fields directly without converting to `bam1_t` structures
- **Existing BGZF I/O**: Uses STAR's existing `bgzf_write()` and `outBAMwriteHeader()` for output

---

## Data Handoff API

### Interface: `SamtoolsSorter`

```cpp
// source/SamtoolsSorter.h
class SamtoolsSorter {
public:
    // Initialize with memory limit and thread count
    SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, Parameters& P);
    
    // Thread-safe: called from coordOneAlign replacement path
    // Takes raw BAM record bytes (same format as coordOneAlign input)
    void addRecord(const char* bamData, uint32_t bamSize, bool hasY);
    
    // Called after all threads complete; performs sort + merge
    void finalize();
    
    // Iterate sorted records for output (post-sort streaming via k-way merge)
    // Returns false when no more records
    bool nextRecord(const char** bamData, uint32_t* bamSize, bool* hasY);
    
    // Cleanup temp files
    ~SamtoolsSorter();

private:
    // Thread-safe buffer with mutex
    // Spill management when buffer exceeds maxRAM
    // K-way merge with QNAME tie-breaking
};
```



### Integration Points

The coordinate-sorted output path flows through `BAMoutput::coordOneAlign()`, called from:

1. **[`source/ReadAlign_outputAlignments.cpp`](source/ReadAlign_outputAlignments.cpp)**:

- Mapped alignments (line 593)
- Unmapped mates with KeepPairs (line 627)
- Unmapped reads (line 661)

2. **[`source/ChimericAlign_chimericBAMoutput.cpp`](source/ChimericAlign_chimericBAMoutput.cpp)** (line 111)

### Branching Strategy

Modified [`source/BAMoutput.cpp`](source/BAMoutput.cpp) `coordOneAlign()` to branch:

```cpp
void BAMoutput::coordOneAlign(char *bamIn, uint bamSize, uint iRead, bool hasY) {
    if (P.outBAMsortMethod == "samtools") {
        // Route to global samtools sorter (initialized in STAR.cpp)
        if (g_samtoolsSorter != nullptr) {
            g_samtoolsSorter->addRecord(bamIn, bamSize, hasY);
        }
        return;
    }
    
    // Existing bin-based sorting logic below...
}
```

This approach:

- Requires minimal changes to call sites
- Centralizes branching in one location
- Reuses existing BAM record format (raw bytes)

---

## Implementation Steps

### Phase 1: Add Backend Selection Parameter

1. **[`source/Parameters.h`](source/Parameters.h)**: Added parameter
   ```cpp
         string outBAMsortMethod;  // "star" (default) or "samtools"
   ```




2. **[`source/Parameters.cpp`](source/Parameters.cpp)**: Parse and validate
3. **[`source/parametersDefault`](source/parametersDefault)**:
   ```javascript
         outBAMsortMethod    star
             string: star or samtools
             star     ... legacy STAR bin-based sorter (default)
             samtools ... samtools-style spill-to-disk sorter
   ```




### Phase 2: Implement SamtoolsSorter API

1. **[`source/SamtoolsSorter.h`](source/SamtoolsSorter.h)**: Interface definition

- `SortKey` struct with numeric fields
- `BAMRecord` struct with move semantics
- `BAMRecordComparator` with QNAME tie-breaking
- `SamtoolsSorter` class with spill-to-disk

2. **[`source/SamtoolsSorter.cpp`](source/SamtoolsSorter.cpp)**: Implementation

- `computeSortKey()`: Parse BAM core fields from raw buffer
- `compareQName()`: Unified QNAME byte comparison helper
- `addRecord()`: Thread-safe buffer with memory check
- `sortAndSpill()`: Sort in-memory buffer and write to temp file
- `finalize()`: Initialize k-way merge
- `nextRecord()`: Stream records from k-way merge with QNAME tie-breaking

3. **Spill-to-disk format**:

- Explicit serialization: `[bamSize:uint32][hasY:uint8][key fields][bamData:bytes]`
- Temp-only format (not a stable ABI)
- Cleanup on destruction

### Phase 3: Implement Backend Branching

1. **[`source/GlobalVariables.h`](source/GlobalVariables.h)**: Declare global sorter
   ```cpp
         extern SamtoolsSorter* g_samtoolsSorter;
   ```




2. **[`source/STAR.cpp`](source/STAR.cpp)**: Initialize before alignment threads
   ```cpp
         if (P.outBAMsortMethod == "samtools" && P.outBAMcoord) {
             g_samtoolsSorter = new SamtoolsSorter(P.limitBAMsortRAM, 
                                                    P.outBAMsortingThreadNactual,
                                                    P.outBAMsortTmpDir, P);
         }
   ```




3. **[`source/BAMoutput.cpp`](source/BAMoutput.cpp)**: Branch in `coordOneAlign()`
4. **[`source/bamSortByCoordinate.cpp`](source/bamSortByCoordinate.cpp)**: Finalization path
   ```cpp
         void bamSortByCoordinate(Parameters &P, ReadAlignChunk **RAchunk, 
                                  Genome &genome, Solo &solo) {
             if (P.outBAMsortMethod == "samtools") {
                 bamSortSamtoolsFinalize(P, genome, solo);
                 return;
             }
             // ... existing STAR bin sorter code ...
         }
   ```




### Phase 4: Post-Sort Y/noY Split with Explicit keepBAM Semantics

**Behavior specification for samtools backend:**

| emitNoYBAM | keepBAM | Output Files |
|------------|---------|--------------|
| no | - | `Aligned.sortedByCoord.out.bam` (primary) |
| yes | no | `_Y.bam` + `_noY.bam` only (no primary) |
| yes | yes | `_Y.bam` + `_noY.bam` + `Aligned.sortedByCoord.out.bam` |

**Implementation:**Uses STAR's existing `bgzf_write()` and `outBAMwriteHeader()` for consistency:

```cpp
void bamSortSamtoolsFinalize(Parameters& P, Genome& genome, Solo& solo) {
    g_samtoolsSorter->finalize();
    
    // Open output handles based on settings
    BGZF* bgzfPrimary = nullptr;
    BGZF* bgzfY = nullptr;
    BGZF* bgzfNoY = nullptr;
    
    // ... open handles ...
    
    // Stream sorted records using bgzf_write
    const char* bamData;
    uint32_t bamSize;
    bool hasY;
    while (g_samtoolsSorter->nextRecord(&bamData, &bamSize, &hasY)) {
        if (bgzfPrimary) bgzf_write(bgzfPrimary, bamData, bamSize);
        if (P.emitNoYBAMyes) {
            bgzf_write(hasY ? bgzfY : bgzfNoY, bamData, bamSize);
        }
    }
    
    // Close handles
}
```

Y-chromosome detection uses `genome.yTids` (populated in `Genome.cpp`) for consistency with existing STAR logic.

### Phase 5: Build System Integration

1. **[`source/Makefile`](source/Makefile)**: Added `SamtoolsSorter.o` to OBJECTS
2. **Compilation verified** on Linux

---

## Testing (Simplified)

**QNAME tie-breaking is part of STAR's sort spec (non-optional).**

### Tier 1: Unit Test

`make test_SamtoolsSorter` runs [`test/SamtoolsSorter_test.cpp`](test/SamtoolsSorter_test.cpp):

- QNAME tie-break: `readA` before `readB` (same numeric keys)
- Unmapped last: tid=-1 records sort after mapped
- Spill handling: tiny maxRAM forces spills, merge still correct

### Tier 2: Fixture-Based Smoke Test

[`tests/run_ychrom_spill_sorter.sh`](tests/run_ychrom_spill_sorter.sh) using `tests/ychrom_test/`:

- Spill run: `--limitBAMsortRAM 100000`
- No-spill run: `--limitBAMsortRAM 1e9`
- Asserts:
  - `samtools quickcheck` passes on `_Y.bam` and `_noY.bam`
  - `@HD` contains `SO:coordinate`
  - Y/noY counts add up
  - Determinism: `samtools view` output matches between spill and no-spill

### Tier 3: scRNA-seq Smoke Test (optional)

[`tests/run_solo_smoke.sh`](tests/run_solo_smoke.sh):

- Self-contained: builds tiny genome, whitelist, and reads
- Validates STARsolo + coordinate sorting works
- No external dependencies

---

## Rationale

### Raw Buffer Sorting

**Why parse BAM core fields directly instead of using `bam1_t`?**

- **Performance**: Avoids allocation/deallocation overhead of `bam1_t` structures
- **Memory efficiency**: Only stores raw BAM bytes + lightweight SortKey metadata
- **Simplicity**: No dependency on htslib parsing APIs beyond BGZF I/O
- **Compatibility**: Works with existing STAR BAM record format (raw bytes)

### QNAME Tie-Breaking

**Why compare QNAME bytes instead of using a hash?**

- **Deterministic ordering**: Byte-level comparison ensures identical ordering for identical inputs
- **Samtools-equivalent**: Matches samtools sort behavior exactly
- **No BAM mutation**: Comparison done on existing QNAME bytes, no aux tags added/stripped
- **Minimal overhead**: Only called when numeric keys are equal (rare case)

### Memory-Bounded Spill-to-Disk

**Why check memory and spill instead of using STAR's bin system?**

- **Fixes memory blow-ups**: STAR's bin sorter can exceed `--limitBAMsortRAM` when many bins are active
- **Predictable behavior**: Spills when threshold exceeded, regardless of bin distribution
- **Efficient merge**: K-way merge streams records without loading all spill files into memory
- **Backward compatible**: Legacy bin sorter remains default, no behavior change for existing users

---

## Rollback

- Set `--outBAMsortMethod star` to use legacy sorter
- Legacy code path remains completely unchanged
- No behavior change for users who don't specify the new parameter

---

## Notes

- **Temp file cleanup**: Handles failures gracefully (cleanup in destructor)
- **Thread-safety**: `addRecord()` uses mutex for thread-safe buffer access
- **Memory accounting**: Approximate (`bamSize + sizeof(BAMRecord) + 64` overhead per record)
- **BAM output**: Uses `bgzf_write()` + `outBAMwriteHeader()` (consistent with existing STAR code)