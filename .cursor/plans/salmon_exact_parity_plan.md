# Salmon Exact Parity Plan

## Problem Statement

The current em_quant implementation cannot achieve exact parity with Salmon because:

1. **Missing EC weights**: The fixture (`eq_classes.txt`) was generated with `--dumpEq` which collapses the per-EC per-transcript weights (`combinedWeights`/`auxs`). Salmon's VB E-step uses `expTheta[tid] * aux`, not just `expTheta[tid] / effLen`.

2. **Salmon explicitly logs this**: When writing without weights, Salmon logs "Collapsing factorization information into simplified equivalence classes." (GZipWriter.cpp:168)

3. **Our E-step is fundamentally different**:
   - **Current**: `weight = expTheta[tid] / eff_lengths[tid]`
   - **Salmon**: `weight = expTheta[tid] * auxs[i]` where `auxs` includes alignment quality, effective length, and range-factorization info

## Key Salmon Code References

### 1. EC Weight Usage in VB E-step (CollapsedEMOptimizer.cpp:144-161)

```cpp
for (size_t i = 0; i < groupSize; ++i) {
    auto tid = txps[i];
    auto aux = auxs[i];  // <-- Per-EC per-transcript weight
    if (expTheta[tid] > 0.0) {
        double v = expTheta[tid] * aux;  // <-- NOT expTheta/effLen
        denom += v;
    }
}
// ... then distribute count proportionally
```

### 2. Single-Transcript Fast Path (CollapsedEMOptimizer.cpp:166-168)

```cpp
} else {
    salmon::utils::incLoop(alphaOut[txps.front()], count);  // Full count, no expTheta
}
```

### 3. Weight File Format (GZipWriter.cpp:217-224)

```cpp
(*equivFilePtr) << groupSize << '\t';
for (uint32_t i = 0; i < groupSize; ++i) {
    (*equivFilePtr) << txps[i] << '\t';
}
for (uint32_t i = 0; i < groupSize; ++i) {
    (*equivFilePtr) << auxs[i] << '\t';  // <-- Weights between IDs and count
}
(*equivFilePtr) << count << '\n';
```

### 4. Salmon Initialization (CollapsedEMOptimizer.cpp:785-831)

```cpp
// NOT unique counts - uses projectedCounts from online phase
alphas[i] = txp.projectedCounts;

// Mix with uniform:
alphas[i] = (alphas[i] * fracObserved) + (uniformPrior * (1.0 - fracObserved));
```

## Required Changes

### Phase 1: Fixture Regeneration

**Task**: Regenerate `test/fixtures/salmon_eq/` with weights

```bash
salmon quant \
    --dumpEqWeights \           # or -d - includes combinedWeights
    --noLengthCorrection \
    --noFragLengthDist \
    --noEffectiveLengthCorrection \
    --threads 1 \               # Deterministic floating-point order
    -i index \
    -l A \
    -1 reads_1.fq -2 reads_2.fq \
    -o salmon_out
```

**Expected format change**:
```
# Before (--dumpEq):
2   0   1   100

# After (--dumpEqWeights):
2   0   1   0.6   0.4   100
```

### Phase 2: EC Data Structure Update

**File**: `source/libem/em_types.h`

```cpp
struct EC {
    std::vector<uint32_t> transcript_ids;
    std::vector<double> weights;  // NEW: per-transcript weights (auxs)
    double count;
    
    EC() : count(0.0) {}
};
```

### Phase 3: EC Loader Update

**File**: `source/libem/ec_loader.cpp`

```cpp
// Detect format by counting fields in first EC line
// Format 1 (no weights): k idx1 idx2 ... idxk count
// Format 2 (weights):    k idx1 idx2 ... idxk w1 w2 ... wk count

// Read transcript indices
for (uint32_t j = 0; j < n_txp; ++j) {
    uint32_t idx;
    if (!(iss >> idx)) throw ...;
    ec.transcript_ids.push_back(idx);
}

// Check if weights are present
std::vector<std::string> remaining;
std::string token;
while (iss >> token) {
    remaining.push_back(token);
}

if (remaining.size() == n_txp + 1) {
    // Weighted format: n_txp weights + 1 count
    ec.weights.reserve(n_txp);
    for (uint32_t j = 0; j < n_txp; ++j) {
        ec.weights.push_back(std::stod(remaining[j]));
    }
    ec.count = std::stod(remaining[n_txp]);
} else if (remaining.size() == 1) {
    // Unweighted format: just count
    ec.count = std::stod(remaining[0]);
    // Default weights = 1/effLen (computed later)
} else {
    throw std::runtime_error("Invalid EC format");
}

// DO NOT sort transcript_ids (or sort weights together)
```

### Phase 4: VB Engine Update

**File**: `source/libem/vb_engine.cpp`

```cpp
// E-step: use weights if available
for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
    const EC& ec = ecs.ecs[ec_idx];
    size_t groupSize = ec.transcript_ids.size();
    
    // Single-transcript fast path (Salmon behavior)
    if (groupSize == 1) {
        uint32_t tid = ec.transcript_ids[0];
        #pragma omp atomic
        expected_counts[tid] += ec.count;  // Full count, no weight
        continue;
    }
    
    // Multi-transcript: use weights
    double denom = 0.0;
    for (size_t i = 0; i < groupSize; ++i) {
        uint32_t tid = ec.transcript_ids[i];
        if (expTheta[tid] > 0.0) {
            double aux = ec.weights.empty() 
                ? (1.0 / state.eff_lengths[tid])  // Fallback for unweighted
                : ec.weights[i];
            denom += expTheta[tid] * aux;
        }
    }
    
    if (denom <= minEQClassWeight) continue;
    
    for (size_t i = 0; i < groupSize; ++i) {
        uint32_t tid = ec.transcript_ids[i];
        if (expTheta[tid] > 0.0) {
            double aux = ec.weights.empty()
                ? (1.0 / state.eff_lengths[tid])
                : ec.weights[i];
            double contribution = ec.count * (expTheta[tid] * aux) / denom;
            #pragma omp atomic
            expected_counts[tid] += contribution;
        }
    }
}
```

### Phase 5: Initialization Update (Optional)

**File**: `source/libem/vb_engine.cpp`

Replace unique-counts seeding with uniform initialization:

```cpp
// Option A: Uniform initialization (simpler, matches --initUniform)
double uniform = 1.0 / state.n;
for (size_t i = 0; i < state.n; ++i) {
    alpha[i] = priorAlphas[i] + uniform * total_reads;
}

// Option B: Match Salmon exactly (requires projected counts from online phase)
// This is harder since we don't have the online phase data
```

### Phase 6: Testing

```bash
# Run both with single thread for determinism
salmon quant --threads 1 ...
./em_quant --vb --threads 1 ...

# Compare
python3 compare_quant.py salmon/quant.sf em_quant/output.tsv --tolerance 0.001
```

## Expected Impact

| Change | Impact on Parity |
|--------|------------------|
| Add EC weights | **Critical** - enables exact E-step math |
| Single-transcript fast path | Moderate - affects small ECs |
| Remove unique-counts init | Moderate - affects convergence path |
| Single-threaded comparison | Eliminates FP ordering noise |

## Implementation Order

1. **Phase 1**: Regenerate fixture with `--dumpEqWeights`
2. **Phase 2**: Update EC struct to hold weights
3. **Phase 3**: Update ec_loader to parse weights
4. **Phase 4**: Update VB E-step to use weights
5. **Phase 5**: Add single-transcript fast path
6. **Phase 6**: Test with `--threads 1`
7. **Phase 7 (optional)**: Update initialization

## Files to Modify

| File | Changes |
|------|---------|
| `test/fixtures/salmon_eq/GENERATE.sh` | Add `--dumpEqWeights --threads 1` |
| `source/libem/em_types.h` | Add `weights` vector to EC struct |
| `source/libem/ec_loader.cpp` | Parse weighted format, don't sort IDs |
| `source/libem/vb_engine.cpp` | Use weights in E-step, add single-txp fast path |
| `source/libem/em_engine.cpp` | Same changes for EM mode |
| `tools/em_quant/run_parity_test.sh` | Add `--threads 1` |

## Acceptance Criteria

- [x] EC struct updated with `weights` vector
- [x] EC loader parses weighted format correctly (auto-detects format)
- [x] VB E-step uses `expTheta * aux` (with fallback to `1/effLen` if no weights)
- [x] EM E-step uses `alpha * aux` (same fallback)
- [x] Single-transcript ECs get full count without weighting (fast path)
- [x] GENERATE.sh updated to use `--dumpEqWeights --threads 1`
- [ ] **PENDING**: Fixture regenerated with weighted format (requires salmon/gffread)
- [ ] **PENDING**: Parity test passes at **< 1% tolerance** for all transcripts

## Current Status

The code implementation is complete. The remaining blocker is regenerating the fixture with `--dumpEqWeights`, which requires `salmon` and `gffread` tools to be available.

**With unweighted fixture (current):**
- 85.7% within 5% tolerance
- 98.6% within 15% tolerance
- Max outlier: 34% (ENST00000476077, 1.2% unique reads)

**Expected with weighted fixture:**
- >99% within 1% tolerance
- All transcripts within 5% tolerance
