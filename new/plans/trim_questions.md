# Questions for Review: Cutadapt-Parity Trimming Plan

**STATUS: ALL RESOLVED** - See decisions below each question.

These questions arose from analysis of the trim plan and STAR-Flex codebase.

---

## 1. Quality Trimming Algorithm

The plan specifies "quality Q20" but doesn't clarify the algorithm.

- **Cutadapt** uses a modified Mott algorithm (dynamic programming to find optimal cutpoint).
- **Trim Galore** wraps cutadapt, so it uses the same algorithm.
- **BWA** and others use a simpler "slide until quality >= threshold" approach.

**Question**: Should we exactly replicate cutadapt's modified Mott algorithm, or is a simpler threshold-based approach acceptable for this milestone?

**RESOLVED**: Replicate cutadapt's modified Mott algorithm so we get the same Q20 behavior as Trim Galore. Simpler heuristics would break parity.

---

## 2. 5' Quality Trimming

The plan includes `quality_trim_head` test case for "Low-quality bases at 5'".

However, Trim Galore defaults do NOT perform 5' quality trimming - only 3'. To get 5' trimming in cutadapt, you must explicitly pass `-q 20,20` (two values).

**Question**: Should 5' quality trimming be included in the default behavior (deviating from Trim Galore), or should it be an optional parameter?

**RESOLVED**: Trim Galore defaults only trim 3'. We'll match that and keep 5' trimming off unless a user explicitly enables it later. (Note: Remove the `quality_trim_head` fixture from scope.)

---

## 3. Adapter Sequences for R1 vs R2

Illumina TruSeq uses different adapter sequences:
- **Read 1**: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
- **Read 2**: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`

Trim Galore auto-detects these. Both share a common prefix `AGATCGGAAGAGC` (13bp).

**Question**: Should we:
- (a) Use separate full adapter sequences for R1/R2 (Trim Galore behavior)?
- (b) Use only the common prefix for both mates (simpler)?
- (c) Make this configurable?

**RESOLVED**: Use the full Illumina TruSeq R1 and R2 adapters separately (as cutadapt does) while allowing future configurability.

---

## 4. Interaction with Existing STAR Clip Parameters

STAR already has adapter clipping parameters:
- `--clip3pAdapterSeq`
- `--clip3pAdapterMMp` (default 0.1)
- `--clip3pNbases`
- `--clipAdapterType` (Hamming, CellRanger4, None)

**Question**: When `--trimWithCutadapt` (proposed flag) is enabled:
- (a) Should it completely replace existing clipping (bypass `ClipMate`)?
- (b) Should it run in addition to existing clipping?
- (c) Should existing clip parameters be errors if `--trimWithCutadapt` is set?

**RESOLVED**: Enabling the new cutadapt-style trimming should replace ClipMate for those mates; we'll bypass ClipMate and document that clip parameters are ignored when this mode is active.

---

## 5. ClipMate Bypass vs Integration

The plan mentions "Update ClipMate accounting so clipped bases propagate to CIGAR/soft-clip output (or bypass ClipMate when the new trimmer runs)."

Current clipping flow in `readLoad.cpp`:
1. Read sequence (line 31)
2. Convert to numeric (line 50)
3. Clip 5' and 3' via ClipMate (lines 60-61)
4. Read quality (line 64)

The plan requires reading quality BEFORE trimming (for quality-based trim).

**Question**: Which approach is preferred:
- (a) Bypass ClipMate entirely when new trimmer is active, record clipped info separately?
- (b) Refactor ClipMate to accept quality input and integrate new logic?
- (c) Run new trimmer in readLoad.cpp before ClipMate, let ClipMate handle CIGAR propagation?

**RESOLVED**: Bypass ClipMate when the new trimmer runs. The trimmer returns clipped lengths and a stats struct that STAR records directly, so we don't have to retrofit ClipMate.

---

## 6. readFilter Handling for Dropped Pairs

When a read pair fails min-length after trimming, Trim Galore discards BOTH mates.

STAR's `readFilter` field currently handles read filtering.

**Question**: For dropped pairs, should we:
- (a) Set `readFilter` to skip both mates (no output)?
- (b) Set a new filter flag specific to length filtering?
- (c) Zero out Lread to signal "don't map"?

**RESOLVED**: When post-trim length <20 on either mate, drop both by marking readFilter='Y' (fail QC) and skip mapping; no new flag needed.

---

## 7. Statistics and Logging

Trim Galore/cutadapt output detailed trimming statistics (reads processed, bases trimmed, reads too short, etc.).

**Question**: Should trimming statistics:
- (a) Be integrated into STAR's existing `Log.final.out`?
- (b) Go to a separate trimming log file?
- (c) Both?

**RESOLVED**: Add a trimming section to Log.final.out (reads processed, trimmed, too short) so users see the info alongside other STAR metrics; no separate log required.

---

## 8. CLI Tool Scope

The plan calls for a standalone CLI tool `cutadapt_lite` for validation.

**Question**: 
- Should this tool be in `tools/` directory (like `flexfilter`)?
- Should it link against `libflex.a` or be a separate `libtrim.a`?
- Should it produce identical FASTQ output format to cutadapt (including read name annotations)?

**RESOLVED**: Put the validator under tools/ as a standalone utility (renamed to `trimvalidate`) linked against the new trimming library (libtrim.a). Its FASTQ output should match cutadapt (same read names/order) but we don't need to mimic cutadapt's read-name annotations.

---

## 9. Thread Safety

STAR processes reads in parallel chunks across multiple threads.

**Question**: The trimming library should be thread-safe. Should we:
- (a) Make functions pure (no global state) - trivially thread-safe?
- (b) Use thread-local storage for any buffers?
- (c) Accept external buffers as parameters?

**RESOLVED**: The trimming API will be pure/stateless; callers pass pointers/buffers, so it's naturally thread-safe (no globals). Any scratch buffers can be stack-local or provided by the caller.

---

## 10. Performance Expectations

The existing `localSearch` function is O(n*m) for read length n and adapter length m.

**Question**: Are there performance requirements? Should we:
- (a) Match existing `localSearch` performance (sufficient)?
- (b) Optimize with SIMD (like CR4's Opal usage)?
- (c) Use bit-parallel algorithms for Hamming distance?

**RESOLVED**: Matching existing localSearch/Hamming performance is sufficient for now. No SIMD/bit tricks needed until we prove it's a bottleneck.

---

## Summary

All questions resolved. See the full implementation plan in `.cursor/plans/cutadapt_parity_trimming_78c03f1c.plan.md`.
