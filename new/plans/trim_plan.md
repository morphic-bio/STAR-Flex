## Cutadapt-Parity Trimming Plan

### Scope
- Match Trim Galore defaults (TruSeq 3′ adapter, `-e 0.1`, `-O 1`, quality Q20, min length 20, paired-end drop-if-short).
- No general cutadapt features (linked adapters, wildcards, indels) in the initial milestone.
- Only paired-end support; single-end can come later if needed.
- Adapter list: Illumina TruSeq, with STAR’s existing polyA handling (`--clipAdapterType CellRanger4`) for tail clipping.

### Library design
1. **Language**: C++11 to match STAR’s build; expose a C-style header (plain structs, `char*`) for easy integration and unit testing.
2. **API**:
   - `TrimResult trim_read(char* seq, char* qual, size_t len, const TrimParams&)`: trims in place, returns new length + flags.
   - `TrimResult trim_pair(Read& r1, Read& r2, const TrimParams&)`: orchestrates paired-end policy (drop pair if either fails min length).
3. **Pipeline**: quality trimming → adapter search/trimming → min-length filter.
4. **Adapter matching**: exact sequence/Hamming search with min overlap 1 and max error 0.1 (i.e. allow floor(len*0.1) mismatches); no indel support initially.
5. **Extensibility hooks**: parameter struct with slots reserved for future features (e.g., multiple adapters, indel mode, polyA detection, SE handling).

### Standalone harness & tests
1. Build a CLI tool (`cutadapt_lite`) that wraps the library; accepts FASTQ pairs, outputs trimmed FASTQ (for local validation against Trim Galore/cutadapt).
2. Fixture matrix (FASTQ pairs generated synthetically; expected outputs from Trim Galore defaults):
   | Case | Description | Expected behavior |
   |------|-------------|-------------------|
   | `clean_long_insert` | Insert length >> read length; no adapters, high quality | No trimming, reads unchanged |
   | `short_insert_overlap` | Insert shorter than read; adapter read-through, high quality | Trim to overlap length (adapters removed) |
   | `short_insert_with_errors` | Overlap + adapters but with ≤10% mismatches | Trim correctly despite mismatches |
   | `quality_trim_tail` | Low-quality tail dropping below Q20 before adapters | Quality trimming removes tail; adapters untouched |
   | `quality_trim_head` | Low-quality bases at 5′ | Quality trimming removes leading bases |
   | `adapter_after_quality_trim` | Low-quality tail exposes adapter | Quality trimming first, then adapter trimming |
   | `below_min_length` | After trims, read length <20 bp | Drop both mates |
   | `polyA_tail` | PolyA tail only (STAR handles via CR4) | No trimming from new library; STAR handles via existing option |
   | `paired_keep_untrimmed` | One mate clean, other needs adapter trimming | Trim only affected mate, keep both if lengths ≥20 |
3. Automate comparisons: run Trim Galore (cutadapt) on each fixture, capture trimmed FASTQs, and diff against our tool; integrate into CI.

### STAR integration
1. Refactor `source/readLoad.cpp` so qualities are read before trimming.
2. Call the new trimming library when a new CLI flag (e.g., `--trimWithCutadapt`) is enabled.
3. Update `ClipMate` accounting so clipped bases propagate to CIGAR/soft-clip output (or bypass ClipMate when the new trimmer runs).
4. Ensure paired-end drop logic is respected (set `readFilter`/flags to skip mapping if pair is discarded).
5. Wire parameters into `Parameters.cpp`/`parametersDefault` (adapter string, Q cutoff, min length, enable flag).

### Milestones
1. **Library MVP**: Adapter + Q trimming functions with unit tests.
2. **CLI validator**: Standalone tool + fixtures verifying parity with Trim Galore defaults.
3. **STAR hook**: Integrate library behind a feature flag; add regression tests using existing STAR test harness.
4. **Extras (future)**: Additional adapters, SE mode, indel matching, multiple quality thresholds, optional unpaired output.
