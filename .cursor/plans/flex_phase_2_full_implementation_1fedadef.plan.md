---
name: Flex Phase 2 Full Implementation
overview: Complete the stubbed flex implementations with minimal upstream churn, keeping all changes guarded by --flexEnable. Implement sample detection, inline hash-to-MEX, CB correction, and inline capture.
todos:
  - id: task1-sample-params
    content: Add sample/probe params to ParametersSolo.h/cpp, Parameters.cpp, parametersDefault
    status: completed
  - id: task1-sample-detector
    content: Replace SampleDetector.cpp stub with full reference implementation
    status: completed
  - id: task2-inline-bundle
    content: Add minimal InlineMatrixBundle and flex methods to SoloFeature.h
    status: completed
  - id: task2-replace-stubs
    content: Replace 4 SoloFeature stub files with reference implementations
    status: completed
  - id: task2-umi-correction
    content: Copy SoloFeature_umiCorrection.cpp and integrate into flex pipeline
    status: completed
  - id: task3-cbcorrector
    content: Wire CbCorrector with init, cleanup, and guarded use in SoloReadBarcode
    status: completed
  - id: task4-readalign
    content: Assess if ReadAlign changes needed; skip if SoloReadFeature suffices
    status: completed
  - id: task5-makefile
    content: Add only required objects; verify no excluded components
    status: completed
  - id: task6-smoke-test
    content: Update smoke test with minimal flags and specific output assertions
    status: completed
  - id: task7-cleanup
    content: Verify all allocations (CbCorrector, SampleDetector) have proper cleanup
    status: completed
---

# Flex Phase 2: Complete Implementation Plan

## Guiding Principles

1. **All changes guarded by `--flexEnable`** - zero behavior change when off
2. **Minimal upstream file diffs** - touch only: Parameters.cpp, ParametersSolo.cpp/h, parametersDefault, SoloReadBarcode_getCBandUMI.cpp, SoloFeature.h, Makefile (ReadAlign only if absolutely required)
3. **No wholesale imports** - only copy what's strictly required
4. **Keep khash** - use khash as hash storage (same as reference)
5. **Required components**: Sample detection, CB correction, UMI correction, inline hash to MEX
6. **Exclude**: Ambiguous CB resolver, TagDominance, PerSamplePipeline, BAM injection, tag-table
7. **Cleanup**: All new allocations (CbCorrector, SampleDetector) must be properly deallocated

---

## Task 1: Sample Detection Parameters and SampleDetector

**Goal**: Add sample/probe params and enable full `SampleDetector.cpp`.

### 1.1 Add Parameters to [ParametersSolo.h](source/ParametersSolo.h)

Use string intermediates for booleans:

```cpp
// Sample detection (flex pipeline) - only consumed when flexEnable
string probeListPath;                      // --soloProbeList
string sampleWhitelistPath;                // --soloSampleWhitelist
string sampleProbesPath;                   // --soloSampleProbes
uint32 sampleProbeOffset = 68;             // --soloSampleProbeOffset (default 68)
string sampleSearchNearbyStr = "yes";      // --soloSampleSearchNearby string
string sampleStrictMatchStr = "no";        // --soloSampleStrictMatch string
bool sampleSearchNearby = true;            // resolved bool (default yes)
bool sampleStrictMatch = false;            // resolved bool (default no)
```

### 1.2 Add Flag Parsing to [Parameters.cpp](source/Parameters.cpp)

```cpp
// Sample detection flags (flex)
parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloProbeList", &pSolo.probeListPath));
parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloSampleWhitelist", &pSolo.sampleWhitelistPath));
parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloSampleProbes", &pSolo.sampleProbesPath));
parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloSampleProbeOffset", &pSolo.sampleProbeOffset));
parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloSampleSearchNearby", &pSolo.sampleSearchNearbyStr));
parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloSampleStrictMatch", &pSolo.sampleStrictMatchStr));
```

### 1.3 Resolve Booleans in [ParametersSolo.cpp](source/ParametersSolo.cpp) initialize()

Resolve unconditionally to avoid uninitialized state, but only consume when flex is on:

```cpp
// Resolve sample detection booleans unconditionally (avoids uninitialized state)
// These are only consumed when flexEnable is on
sampleSearchNearby = (sampleSearchNearbyStr == "yes");
sampleStrictMatch = (sampleStrictMatchStr == "yes");
```

### 1.4 Add Defaults to [parametersDefault](source/parametersDefault)

```
soloProbeList               -
soloSampleWhitelist         -
soloSampleProbes            -
soloSampleProbeOffset       68
soloSampleSearchNearby      yes
soloSampleStrictMatch       no
```

### 1.5 Replace SampleDetector.cpp Stub

Replace [source/SampleDetector.cpp](source/SampleDetector.cpp) with full implementation from `reference/source/SampleDetector.cpp`. Only instantiated when `flexEnable` is on.

**Files**: ParametersSolo.h (+8), ParametersSolo.cpp (+3), Parameters.cpp (+6), parametersDefault (+6), SampleDetector.cpp (replace)

---

## Task 2: Inline Hash to Flex to MEX (with UMI Correction)

**Goal**: Implement InlineMatrixBundle, replace stubs, include UMI correction in flow.

### 2.1 Add to [SoloFeature.h](source/SoloFeature.h)

Keep header small - only add methods/structs that are actually invoked:

```cpp
// Inline hash matrix bundle (flex pipeline)
struct InlineMatrixBundle {
    std::vector<std::string> barcodes;
    std::vector<std::string> features;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> triplets; // row, col, count
};

// Flex pipeline methods (only when flexEnable) - add only what's invoked
InlineMatrixBundle buildInlineMatrixFromHash();
void runFlexFilterInline(const InlineMatrixBundle& bundle);
void writeMexFromInlineHashDedup(const std::string& outDir, const InlineMatrixBundle& bundle);
void umiCorrectionFromHash();  // UMI correction for flex path
```

Do NOT add loader/sink unless required by the above methods.

### 2.2 Replace Stub Files with Reference Implementations

All must check `flexEnable` at entry and return early if false:

- [SoloFeature_flexfilter.cpp](source/SoloFeature_flexfilter.cpp) - from `reference/source/SoloFeature_flexfilter.cpp`
- [SoloFeature_writeMexFromInlineHashDedup.cpp](source/SoloFeature_writeMexFromInlineHashDedup.cpp) - from reference
- [SoloFeature_materializeFromHash.cpp](source/SoloFeature_materializeFromHash.cpp) - from reference
- [SoloFeature_collapseUMI_fromHash.cpp](source/SoloFeature_collapseUMI_fromHash.cpp) - from reference

### 2.3 Add UMI Correction (Required)

Copy `reference/source/SoloFeature_umiCorrection.cpp` to `source/`. Guard execution with `flexEnable`. Integrate into flex pipeline so UMI correction happens before MEX write.

### 2.4 Assess Loader/Sink

Only copy `SoloReadInfoLoader.cpp/h` or `SoloReadInfoSink.cpp/h` if the above methods require them. Otherwise skip.

**Do NOT copy**: `SoloFeature_resolveAmbiguousCBs.cpp`, `TagDominanceCheck.*`

**Files**: SoloFeature.h (+15-20), 4 stubs replaced, SoloFeature_umiCorrection.cpp (new)

---

## Task 3: CB Correction Under flexEnable

**Goal**: Wire CbCorrector into barcode parsing when flex is enabled, with proper cleanup.

### 3.1 Add CbCorrector Pointer to [ParametersSolo.h](source/ParametersSolo.h)

```cpp
// Forward declare at top of file
class CbCorrector;

// Add inside class (near flexEnable members)
CbCorrector* cbCorrector = nullptr;
```

### 3.2 Initialize CbCorrector in [ParametersSolo.cpp](source/ParametersSolo.cpp)

At end of `initialize()`, guarded:

```cpp
// Initialize CbCorrector for flex pipeline
if (flexEnable && cbWLyes) {
    cbCorrector = new CbCorrector(cbWLstr);
    pP->inOut->logMain << "CbCorrector initialized for flex pipeline\n";
}
```

### 3.3 Add Cleanup in ParametersSolo Destructor

Ensure no memory leak:

```cpp
// In destructor or cleanup method
if (cbCorrector) {
    delete cbCorrector;
    cbCorrector = nullptr;
}
```

### 3.4 Wire into [SoloReadBarcode_getCBandUMI.cpp](source/SoloReadBarcode_getCBandUMI.cpp)

Add minimal, guarded CB correction - do not alter non-flex code paths:

```cpp
// Flex CB correction (guarded) - minimal addition
if (P.pSolo.flexEnable && P.pSolo.cbCorrector) {
    CbMatch match = P.pSolo.cbCorrector->correct(cbSeq);
    if (match.whitelistIdx != 0) {
        cbMatchInd.clear();
        cbMatchInd.push_back(match.whitelistIdx - 1);
        cbMatch = (match.hammingDist == 0) ? 0 : 1;
    }
}
```

**Files**: ParametersSolo.h (+3), ParametersSolo.cpp (+8 including cleanup), SoloReadBarcode_getCBandUMI.cpp (+8 guarded)

---

## Task 4: ReadAlign Changes (Only If Required)

**Goal**: Avoid touching ReadAlign unless inline hash absolutely needs per-read state.

### 4.1 Assessment

Check if flex pipeline can reuse existing `SoloReadFeature` data (cbMatchInd, umiSeq, etc.) instead of adding new members to ReadAlign.

**If SoloReadFeature suffices**: Skip ReadAlign.h and ReadAlign_outputAlignments.cpp changes entirely. The existing hook is sufficient.

**If per-read state is required**: Add only the minimum:

```cpp
// In ReadAlign.h - only if truly needed
uint32_t extractedCbIdxPlus1_ = 0;   // 1-based CB index, 0 = no match
uint32_t extractedUmi24_ = 0;         // 24-bit packed UMI
```

### 4.2 Preference

Prefer NOT touching ReadAlign. Keep sample detection in SoloFeature/SampleDetector, not ReadAlign.

**Files**: Potentially 0 if SoloReadFeature suffices; otherwise ReadAlign.h (+2), ReadAlign_outputAlignments.cpp (+5 guarded)

---

## Task 5: Update Makefile

### 5.1 Add Only Required Objects to [Makefile](source/Makefile)

```makefile
FLEX_OBJECTS = ... \
    SoloFeature_umiCorrection.o
    # Add loader/sink ONLY if Task 2.4 determines they're required
```

### 5.2 Verification Checklist

Before committing, verify Makefile does NOT include:

- TagDominance objects
- PerSamplePipeline objects
- BAM injection objects (BAMTagBinaryWriter, BAMunsortedAddSoloTags)
- Tag-table objects (SoloFeature_writeTagTable)

### 5.3 Build Test

Run `make clean && make -j4` and verify clean compilation.

**Files**: Makefile (+1-3 lines)

---

## Task 6: Smoke Test

### 6.1 Update [run_flex_smoke.sh](tests/flex_smoke/run_flex_smoke.sh)

Specify minimal flag set:

```bash
# Test 1: Default (flex off) - must match upstream behavior
$STAR_BIN --genomeDir "$FIXTURE/genome" \
          --readFilesIn "$FIXTURE/reads_R2.fastq" "$FIXTURE/reads_R1.fastq" \
          --outFileNamePrefix "$OUT/flex_off/" \
          --soloType CB_UMI_Simple \
          --soloCBwhitelist "$FIXTURE/whitelist.txt" \
          --soloCBlen 16 --soloUMIlen 12 \
          --soloFeatures Gene
# Assert: exit 0, no flex-specific output

# Test 2: Flex enabled
$STAR_BIN --genomeDir "$FIXTURE/genome" \
          --readFilesIn "$FIXTURE/reads_R2.fastq" "$FIXTURE/reads_R1.fastq" \
          --outFileNamePrefix "$OUT/flex_on/" \
          --flexEnable yes \
          --soloType CB_UMI_Simple \
          --soloCBwhitelist "$FIXTURE/whitelist.txt" \
          --soloCBlen 16 --soloUMIlen 12 \
          --soloFeatures Gene
# Assert: exit 0, MEX files present (if fixture supports), or no crash
```

### 6.2 Expected Outputs

| Scenario | Exit Code | Output Check |

|----------|-----------|--------------|

| `--flexEnable no` | 0 | Standard Solo.out, no flex-specific dirs |

| `--flexEnable yes` | 0 | MEX files if fixture supports, else no crash |

**Files**: run_flex_smoke.sh (update)

---

## Task 7: Cleanup Verification

### 7.1 Allocations to Track

| Allocation | Location | Cleanup |

|------------|----------|---------|

| `CbCorrector` | ParametersSolo.cpp | ParametersSolo destructor |

| `SampleDetector` | Where instantiated | Corresponding destructor |

### 7.2 Verify No Leaks

Run with valgrind or ASAN if available to confirm no memory leaks from flex allocations.

---

## File Change Summary

| File | Change Type | Lines | Notes |

|------|-------------|-------|-------|

| `ParametersSolo.h` | Modified | +11 | Sample params + CbCorrector ptr |

| `ParametersSolo.cpp` | Modified | +11 | Init + cleanup |

| `Parameters.cpp` | Modified | +6 | Flag parsing |

| `parametersDefault` | Modified | +6 | Defaults |

| `SoloFeature.h` | Modified | +15-20 | InlineMatrixBundle + methods |

| `SoloReadBarcode_getCBandUMI.cpp` | Modified | +8 | CB correction, guarded |

| `ReadAlign.h` | Maybe | +2 | Only if needed |

| `ReadAlign_outputAlignments.cpp` | Maybe | +5 | Only if needed |

| `Makefile` | Modified | +1-3 | New objects |

| `SampleDetector.cpp` | Replace | Full | From reference |

| 4 SoloFeature stubs | Replace | Full | From reference |

| `SoloFeature_umiCorrection.cpp` | New | Full | From reference |

**Upstream files with diffs**: 6-8 (all guarded by flexEnable)

---

## Out of Scope (Deferred)

- Ambiguous CB resolution (`SoloFeature_resolveAmbiguousCBs.cpp`)
- TagDominance checking
- PerSamplePipeline
- BAM tag injection (`BAMTagBinaryWriter`, `BAMunsortedAddSoloTags`)
- Tag table sidecar (`SoloFeature_writeTagTable.cpp`)
- SoloReadInfoLoader/Sink (unless required by Task 2)