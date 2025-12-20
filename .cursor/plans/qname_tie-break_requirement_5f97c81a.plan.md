---
name: QNAME Tie-Break Requirement
overview: Elevate QNAME tie-breaking from implementation detail to a hard requirement in STAR's coordinate sorting specification, with corresponding updates to plan documentation and code comments.
todos:
  - id: qname-spec-plan
    content: Add QNAME tie-breaking spec statement to plan document
    status: pending
  - id: qname-spec-code
    content: Add spec comment header and update comparator comments in SamtoolsSorter.h
    status: pending
  - id: qname-safety
    content: Add assert guards for bamPtr validity in HeapLess comparator
    status: pending
---

# QNAME Tie-Breaking as Hard Requirement

## Summary

Update plan documentation and code comments to establish QNAME byte comparison as a non-negotiable requirement for coordinate sorting when numeric keys are equal. Add safety guards for `bamPtr` validity.

---

## 1. Plan Document Updates

Update [`source/.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md`](.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md):

### Add Spec Statement (in Stage 0 or new "Sorting Specification" section)

```markdown
### Sorting Specification (REQUIRED)

**QNAME tie-breaking is mandatory for coordinate sorting:**

When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering MUST be deterministic by QNAME byte comparison. This is a non-negotiable requirement for reproducibility and is part of STAR's coordinate sort specification going forward.

- Applies to all sources: in-memory records and spill file records
- QNAME bytes compared lexicographically (excluding trailing NUL)
- Ensures identical input produces identical output ordering
```

### Update Rationale Section

Add explicit statement that QNAME tie-breaking is required (not optional) for reproducibility.

---

## 2. Code Comment Updates

### [`source/SamtoolsSorter.h`](source/SamtoolsSorter.h)

**Add spec comment at top of file (after includes):**

```cpp
// =============================================================================
// SORTING SPECIFICATION (REQUIRED):
// When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering
// MUST be deterministic by QNAME byte comparison. This is a non-negotiable
// requirement for reproducibility and is part of STAR's coordinate sort spec.
// =============================================================================
```

**Update `SamtoolsSorterHelpers::compareQName` comment:**

```cpp
// REQUIRED: Compare QNAME bytes from raw BAM buffers for tie-breaking.
// Called when all numeric SortKey fields are equal.
// This is mandatory for deterministic ordering - not optional.
```

**Update `HeapLess` comparator comment:**

```cpp
// Comparator for min-heap with REQUIRED QNAME tie-breaking.
// INVARIANT: bamPtr must be non-null for any entry in the heap.
// When numeric keys are equal, QNAME byte comparison determines order.
```

---

## 3. Safety Guard for `bamPtr`

### [`source/SamtoolsSorter.h`](source/SamtoolsSorter.h) - HeapLess comparator

Add assert/guard before QNAME comparison:

```cpp
struct HeapLess {
    bool operator()(const HeapEntry& a, const HeapEntry& b) const {
        // Compare numeric SortKey first
        if (a.key < b.key) return false;
        if (b.key < a.key) return true;
        // Numeric keys equal - QNAME tie-break is REQUIRED
        // INVARIANT: bamPtr must be valid for all heap entries
        assert(a.bamPtr != nullptr && "HeapEntry::bamPtr must not be null");
        assert(b.bamPtr != nullptr && "HeapEntry::bamPtr must not be null");
        return SamtoolsSorterHelpers::compareQName(a.bamPtr, b.bamPtr) > 0;
    }
};
```

Need to add `#include <cassert>` at top of file.

---

## 4. Files to Modify

| File | Changes |
|------|---------|
| [`.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md`](.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md) | Add sorting spec section, update rationale |
| [`source/SamtoolsSorter.h`](source/SamtoolsSorter.h) | Add spec comment header, update comparator comments, add assert guards |

---

## 5. Verification

After changes:
- Compilation succeeds (assert available)
- Comments clearly state QNAME tie-breaking is required
- Plan and code are in sync on this requirement
