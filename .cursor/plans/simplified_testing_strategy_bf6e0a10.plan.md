---
name: Simplified Testing Strategy
overview: "Replace the current testing strategy in the samtools backend plan with a lighter three-tier approach: unit tests, functional smoke tests, and targeted spill stress tests."
todos:
  - id: update-testing-strategy
    content: Replace Testing Strategy section in plan with three-tier approach
    status: pending
---

# Simplified Testing Strategy

## Summary

Update the Testing Strategy section in [`samtools_bam_sort_backend_f3c39581.plan.md`](.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md) to use a three-tier approach that matches the lighter internal sorter implementation.---

## New Testing Strategy

### Tier 1: Unit/Logic Sanity (fast)

Small synthetic BAM buffers or minimal fixtures that exercise core sorting logic:

```cpp
// tests/samtools_sorter_unit.cpp or similar test harness
// Test cases:
// 1. Ordering by tid/pos (records with different coordinates)
// 2. Tie-break by QNAME (same tid/pos, different QNAME)
// 3. Unmapped last (tid=-1 records sort after mapped)
// 4. Spill threshold triggers (add records until spill occurs)
```

**Validation:**

- `nextRecord()` returns records in expected order
- Spill files created when memory threshold exceeded
- All records accounted for after merge

---

### Tier 2: Functional Output Checks (smoke)

One bulk dataset, one scRNA-seq dataset:

```bash
# Bulk RNA-seq
./source/STAR ... --outBAMsortMethod samtools --outFileNamePrefix test_bulk_

# scRNA-seq (Solo)
./source/STAR ... --outBAMsortMethod samtools --soloType CB_UMI_Simple ...
```

**Validation:**

- `samtools quickcheck` passes
- `@HD` header contains `SO:coordinate`
- Unmapped reads at end (spot-check: `samtools view | tail -100 | cut -f3`)
- Deterministic output: run twice with same input, diff outputs match

---

### Tier 3: Spill Stress (targeted)

Force spills with small RAM limit, compare to high RAM run:

```bash
Chec# Low RAM (force spills)
./source/STAR ... --outBAMsortMethod samtools \
    --limitBAMsortRAM 100000000 \
    --outFileNamePrefix test_spill_

# High RAM (no spills)
./source/STAR ... --outBAMsortMethod samtools \
    --limitBAMsortRAM 10000000000 \
    --outFileNamePrefix test_normal_

# Compare outputs
diff <(samtools view test_spill_Aligned.sortedByCoord.out.bam) \
     <(samtools view test_normal_Aligned.sortedByCoord.out.bam)

# Verify counts match
samtools flagstat test_spill_Aligned.sortedByCoord.out.bam > spill_flagstat.txt
samtools flagstat test_normal_Aligned.sortedByCoord.out.bam > normal_flagstat.txt
diff spill_flagstat.txt normal_flagstat.txt
```

---

## Dropped Tests

- **STAR vs legacy sorter diff**: Not required (different implementations may have different tie-break behavior)
- **samtools ordering diff**: Not a requirement (we define our own deterministic spec)

---

## Changes to Plan Document

Replace the current "Testing Strategy" section (lines ~280-380) with the three-tier approach above. Keep it concise and actionable.---