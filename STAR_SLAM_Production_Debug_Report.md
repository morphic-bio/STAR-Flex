# STAR-SLAM Production Sample Debug Analysis Report

**Date:** January 9, 2026  
**Sample:** WDHD1-0h-3  
**Analysis Type:** Debug instrumentation for discordant genes  
**STAR Version:** 2.7.11b (compiled 2026-01-09)

---

## Executive Summary

The debug run was rerun with the same trimming flags used in the GEDI pipeline (EndToEnd + 3' adapter clipping). The recomputed top-20 discordant genes remain **highly divergent**, so the discrepancies are **not explained by trimming** alone. The divergence is **bidirectional**: STAR NTR is higher than GEDI for 12 genes and lower for 8. A focused single‑mapping subset (avg weight >= 0.95) still shows large NTR deltas, indicating that multi‑mapping is **not** the sole driver.

**Key Findings (trimmed, recomputed list):**
- **SNP mask filtering:** Not a factor (0 drops)
- **Strandness filtering:** Not a factor (0 drops)
- **Multi-mapping:** Present in some genes, but several high-delta genes are single‑mapping
- **Conversion counts:** Large deltas remain; STAR is often higher than GEDI
- **Position bias:** Some genes have conversions concentrated at read starts, others at read ends

---

## Weight Mode Comparison Experiment

**Date:** January 9, 2026  
**Purpose:** Compare `--slamWeightMode Alignments` (default) vs `--slamWeightMode Uniform` to determine which weighting scheme produces results closer to GEDI, using the same flags as the trimmed baseline run (no extra `outFilter*` overrides).

### Test Configuration

Both runs used identical parameters except for `--slamWeightMode`:

```bash
--slamQuantMode 1
--slamSnpDetect 1
--slamConvRate 0.02024
--slamStrandness Sense
--alignEndsType EndToEnd
--clip3pAdapterSeq AGATCGGAAGAG
--clip3pAdapterMMp 0.1
--slamWeightMode <Alignments|Uniform>
```

**Genome Index:** `/mnt/pikachu/STAR-Flex/test/fixtures/slam/ref/star_index`  
**Input FASTQ:** `/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz`  
**GEDI Reference:** `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_rerun.tsv.gz`

### Weight Mode Definitions

- **Alignments (default):** `weight = 1.0 / nTr` - Multi-mapped reads are down-weighted by the number of alignments
- **Uniform:** `weight = 1.0` - All alignments receive full weight regardless of multi-mapping status

### Correlation Results

| Weight Mode | Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-------------|-----------|--------|---------|-------------|--------------|--------------|---------------|
| **Alignments** | >=20 | readcount | 4255 | **0.920467** | **0.924726** | **0.905483** | **0.931271** |
| **Alignments** | >=50 | readcount | 2166 | **0.926568** | **0.929443** | **0.930352** | **0.936203** |
| **Alignments** | >=100 | readcount | 1149 | **0.906414** | **0.910138** | **0.824392** | **0.917799** |
| **Uniform** | >=20 | readcount | 4255 | 0.914377 | 0.919544 | 0.881119 | 0.926008 |
| **Uniform** | >=50 | readcount | 2166 | 0.919508 | 0.922726 | 0.903878 | 0.929250 |
| **Uniform** | >=100 | readcount | 1149 | 0.895682 | 0.900634 | 0.771093 | 0.908097 |

### Key Observations

1. **Alignments mode is closer to GEDI:** Across all thresholds, Alignments mode shows higher correlations than Uniform mode:
   - NTR Pearson: **0.920** (Alignments) vs **0.914** (Uniform) at threshold >=20
   - NTR Spearman: **0.925** (Alignments) vs **0.920** (Uniform) at threshold >=20
   - k/nT Pearson: **0.905** (Alignments) vs **0.881** (Uniform) at threshold >=20
   - k/nT Spearman: **0.931** (Alignments) vs **0.926** (Uniform) at threshold >=20

2. **Uniform mode over-counts multi-mapped reads:** Uniform mode shows higher ReadCount values for multi-mapped genes (e.g., RPL17: 30,383 vs 4,849 in GEDI), indicating that giving full weight to all alignments inflates counts.

3. **Correlation stays high but dips at >=100:** Both modes remain high at >=50, with a modest drop at >=100 (fewer genes), consistent with residual discordance in a subset of highly expressed genes.

4. **ReadCount discrepancies:** 
   - Alignments mode: 648 genes with ReadCount mismatches
   - Uniform mode: 2,836 genes with ReadCount mismatches (more mismatches)

5. **Conversion discrepancies:**
   - Alignments mode: 916 genes with Conversion mismatches
   - Uniform mode: 1,527 genes with Conversion mismatches (more mismatches)

6. **Alignment mode matches baseline:** The Alignments-mode correlations now match the trimmed baseline comparison (see “Full Comparison Stats”), confirming the earlier drop was driven by mismatched filter flags.

### Interpretation

**Alignments mode (default) is closer to GEDI** than Uniform mode. The down-weighting of multi-mapped reads (`weight = 1.0 / nTr`) produces more accurate quantification compared to giving full weight to all alignments. Uniform mode over-counts multi-mapped reads, leading to inflated ReadCount and Conversion values, particularly for highly multi-mapped genes like ribosomal proteins.

**Recommendation:** Continue using `--slamWeightMode Alignments` (default) as it provides better parity with GEDI.

### Output Files

- **Alignments mode output:** `/storage/SLAM-Seq-prod-compare-20260109/star_weight_alignments_correct/WDHD1_0h3_alignments_SlamQuant.out`
- **Uniform mode output:** `/storage/SLAM-Seq-prod-compare-20260109/star_weight_uniform_correct/WDHD1_0h3_uniform_SlamQuant.out`
- **Comparison results:** 
  - `/storage/SLAM-Seq-prod-compare-20260109/compare_weight_alignments_correct.txt`
  - `/storage/SLAM-Seq-prod-compare-20260109/compare_weight_uniform_correct.txt`

---

## Test Configuration (trimmed debug run)

```bash
--slamQuantMode 1
--slamSnpDetect 1
--slamConvRate 0.02024
--slamStrandness Sense
--alignEndsType EndToEnd
--clip3pAdapterSeq AGATCGGAAGAG
--clip3pAdapterMMp 0.1
--slamDebugGeneList /mnt/pikachu/STAR-Flex/top20_genes.txt
--slamDebugOutPrefix /storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed/WDHD1_0h3
--slamDebugMaxReads 2000
```

**Genome Index:** `/mnt/pikachu/STAR-Flex/test/fixtures/slam/ref/star_index`  
**Input FASTQ:** `WDHD1-0h-3_S201_R1_001.fastq.gz`

---

## GEDI Error-Model Parameters (lock for parity)

From `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_rerun.param`:

- `err=0.001`
- `errlm=TA*1.434`
- `errlm2=TA*1.434`
- `conv=0.02024` (from `WDHD1_0h3_Sense_rerun.rates.tsv`, `single_new`)

Action: for future GEDI reruns, pass `-err 0.001 -errlm 'TA*1.434' -errlm2 'TA*1.434' -conv 0.02024` to remove hidden variability. Note: locking `-conv` disables estimation; we should remove it after the remaining divergences are resolved so conversion is learned from data again.

---

## Recomputed Top-20 Discordant Genes (trimmed)

List recomputed from trimmed STAR output vs GEDI, using |NTR delta| with GEDI readcount >= 50.

| Gene | Symbol | GEDI NTR | STAR NTR | ΔNTR | GEDI Conv | STAR Conv | AvgWeight | IntronicFrac |
|------|--------|----------|----------|------|-----------|-----------|-----------|-------------|
| ENSG00000172053 | QARS1 | 0.087200 | 0.999788 | +0.912588 | 18.000 | 151.000 | 0.997 | 0.003 |
| ENSG00000233954 | UQCRHL | 0.036400 | 0.922599 | +0.886199 | 3.000 | 32.000 | 0.970 | 0.000 |
| ENSG00000254726 | MEX3A | 0.201400 | 0.999962 | +0.798562 | 3.000 | 11.000 | 0.998 | 0.000 |
| ENSG00000197061 | H4C3 | 0.215900 | 0.999997 | +0.784097 | 269.500 | 1195.500 | 0.998 | 0.000 |
| ENSG00000229833 | PET100 | 0.219600 | 0.999999 | +0.780399 | 40.000 | 198.200 | 0.507 | 0.000 |
| ENSG00000140525 | FANCI | 1.000000 | 0.237848 | -0.762152 | 11.000 | 3.000 | 0.814 | 0.132 |
| ENSG00000102230 | PCYT1B | 1.000000 | 0.253975 | -0.746025 | 12.000 | 4.000 | 0.747 | 0.193 |
| ENSG00000196365 | LONP1 | 0.254400 | 0.999998 | +0.745598 | 14.000 | 57.000 | 0.994 | 0.002 |
| ENSG00000170540 | ARL6IP1 | 1.000000 | 0.300216 | -0.699784 | 11.000 | 3.000 | 0.987 | 0.000 |
| ENSG00000169592 | INO80E | 0.731900 | 0.049424 | -0.682476 | 13.000 | 2.000 | 0.991 | 0.000 |
| ENSG00000147403 | RPL10 | 0.211800 | 0.873707 | +0.661907 | 165.500 | 610.067 | 0.413 | 0.000 |
| ENSG00000181392 | SYNE4 | 0.338500 | 0.999999 | +0.661499 | 9.000 | 7.000 | 0.529 | 0.471 |
| ENSG00000265681 | RPL17 | 0.001300 | 0.644677 | +0.643377 | 26.780 | 725.945 | 0.160 | 0.003 |
| ENSG00000173598 | NUDT4 | 0.263400 | 0.903623 | +0.640223 | 3.533 | 11.533 | 0.532 | 0.137 |
| ENSG00000100084 | HIRA | 0.000000 | 0.622129 | +0.622129 | 0.000 | 9.000 | 0.457 | 0.525 |
| ENSG00000161634 | DCD | 0.726400 | 0.133520 | -0.592880 | 33.080 | 10.583 | 0.759 | 0.000 |
| ENSG00000075618 | FSCN1 | 0.170200 | 0.760784 | +0.590584 | 37.000 | 164.500 | 0.944 | 0.000 |
| ENSG00000106211 | HSPB1 | 0.013300 | 0.598491 | +0.585191 | 10.000 | 76.000 | 1.000 | 0.000 |
| ENSG00000084072 | PPIE | 0.955200 | 0.377594 | -0.577606 | 16.000 | 7.000 | 0.946 | 0.054 |
| ENSG00000198356 | GET3 | 0.010400 | 0.587852 | +0.577452 | 2.000 | 17.000 | 0.990 | 0.005 |

---

## Single‑Mapping Discordant Genes (avg weight >= 0.95)

These genes remain highly discordant even when multi‑mapping is minimal.

| Gene | Symbol | GEDI NTR | STAR NTR | ΔNTR | GEDI Conv | STAR Conv | AvgWeight | IntronicFrac | ConvFrac<=5 | MedianPos |
|------|--------|----------|----------|------|-----------|-----------|-----------|-------------|-------------|-----------|
| ENSG00000172053 | QARS1 | 0.087200 | 0.999788 | +0.912588 | 18.000 | 151.000 | 0.997 | 0.003 | 0.033 | 46.0 |
| ENSG00000233954 | UQCRHL | 0.036400 | 0.922599 | +0.886199 | 3.000 | 32.000 | 0.970 | 0.000 | 0.000 | 45.0 |
| ENSG00000254726 | MEX3A | 0.201400 | 0.999962 | +0.798562 | 3.000 | 11.000 | 0.998 | 0.000 | 0.273 | 13.0 |
| ENSG00000197061 | H4C3 | 0.215900 | 0.999997 | +0.784097 | 269.500 | 1195.500 | 0.998 | 0.000 | 0.927 | 2.0 |
| ENSG00000196365 | LONP1 | 0.254400 | 0.999998 | +0.745598 | 14.000 | 57.000 | 0.994 | 0.002 | 0.000 | 46.0 |
| ENSG00000170540 | ARL6IP1 | 1.000000 | 0.300216 | -0.699784 | 11.000 | 3.000 | 0.987 | 0.000 | 0.000 | 46.0 |
| ENSG00000169592 | INO80E | 0.731900 | 0.049424 | -0.682476 | 13.000 | 2.000 | 0.991 | 0.000 | 0.923 | 2.0 |
| ENSG00000106211 | HSPB1 | 0.013300 | 0.598491 | +0.585191 | 10.000 | 76.000 | 1.000 | 0.000 | 0.961 | 1.0 |
| ENSG00000198356 | GET3 | 0.010400 | 0.587852 | +0.577452 | 2.000 | 17.000 | 0.990 | 0.005 | 1.000 | 0.0 |

**Position bias signals:**
- Strong **early‑position** conversion enrichment: H4C3, HSPB1, GET3, INO80E.
- Strong **late‑position** conversion enrichment: QARS1, LONP1, UQCRHL, ARL6IP1.

This suggests systematic positional differences may be contributing to the discordance, independent of multi‑mapping.

---

## Mismatch Definition Comparison (GEDI vs STAR)

Used `tests/slam/compare_mismatch_summaries.py` on trimmed outputs.

Output: `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_mismatch_summaries.txt`

**Result:** Large mismatches in **Intronic / IntronicSense** categories for both mismatches.tsv and mismatchdetails.tsv. This aligns with our prior decision to retain STAR’s intronic definition rather than GEDI’s. Exonic categories are closer but still show differences.

---

## Full Comparison Stats (trimmed STAR vs GEDI)

Reference: `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_rerun.tsv.gz`  
Test: `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed/SlamQuant.out`  
Output: `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_debug.txt`

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20 | readcount | 4255 | 0.920467 | 0.924726 | 0.905483 | 0.931271 |
| >=50 | readcount | 2166 | 0.926568 | 0.929443 | 0.930352 | 0.936203 |
| >=100 | readcount | 1149 | 0.906414 | 0.910138 | 0.824392 | 0.917799 |

---

## Fixed-Parameter GEDI Rerun (err/conv locked)

Reference: `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_fixed.tsv.gz`  
Output: `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_fixed.txt`

**Result:** The correlation metrics and mismatch counts are unchanged from the trimmed run above. Locking `err`/`errlm`/`errlm2`/`conv` does not explain the divergence.

---

## Read-Level Divergence Trace (GEDI debug vs STAR debug)

We rebuilt the GEDI jar with debug instrumentation and added a `ReadLoc` key to STAR’s debug output (chromosome coordinates, end‑exclusive) to enable per‑read matching.

**STAR debug output:** `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed_loc/WDHD1_0h3.reads.tsv`  
**GEDI debug output:** `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_fixed.debug.tsv`  
**Comparison output:** `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_debug_reads_top20.tsv`

Summary from `tests/slam/compare_debug_reads.py` (top20 genes):

- STAR keys: 591
- GEDI keys: 769
- STAR-only (pass weight > 0): 53
- GEDI-only (consistent weight > 0): 166
- Mismatches (|delta| >= 0.1): 649

**Pattern:** The largest deltas are dominated by **weight differences at the same read locations**, not by STAR drop reasons (top mismatches are all `PASS`). GEDI weights are often higher for the same locus, which points to **assignment/weighting differences** (e.g., distinct counts or transcript compatibility) rather than filtering differences.

---

## Diagnostic Statistics (trimmed debug run)

From `SlamQuant.out.diagnostics`:

- Reads processed: **1,499,512**
- Reads dropped SNP mask: **0**
- Reads dropped strandness: **0**
- Reads with nAlignWithGene=0: **1,602,153**
- Reads with sumWeight < 1.0: **1,817,144**

---

## Files Generated

**Trimmed Debug Outputs:**
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed/WDHD1_0h3.gene.tsv`
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed/WDHD1_0h3.reads.tsv`
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed/SlamQuant.out.diagnostics`

**Read-Location Debug Outputs:**
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_trimmed_loc/WDHD1_0h3.reads.tsv`
- `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_fixed.debug.tsv`

**Single‑Mapping Debug Outputs:**
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_singlemap/WDHD1_0h3.gene.tsv`
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_singlemap/WDHD1_0h3.reads.tsv`
- `/storage/SLAM-Seq-prod-compare-20260109/debug/WDHD1_0h3_singlemap/SlamQuant.out.diagnostics`

**Comparison Outputs:**
- `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_debug.txt`
- `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_mismatch_summaries.txt`
- `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_Sense_trimmed_fixed.txt`
- `/storage/SLAM-Seq-prod-compare-20260109/compare/compare_debug_reads_top20.tsv`

**Weight Mode Outputs (corrected flags):**
- `/storage/SLAM-Seq-prod-compare-20260109/star_weight_alignments_correct/WDHD1_0h3_alignments_SlamQuant.out`
- `/storage/SLAM-Seq-prod-compare-20260109/star_weight_uniform_correct/WDHD1_0h3_uniform_SlamQuant.out`
- `/storage/SLAM-Seq-prod-compare-20260109/compare_weight_alignments_correct.txt`
- `/storage/SLAM-Seq-prod-compare-20260109/compare_weight_uniform_correct.txt`

---

## Conclusion

Trimming was necessary for consistent instrumentation, but it does **not** resolve the STAR vs GEDI divergence. The single‑mapping subset remains highly discordant, which points to **conversion counting differences and/or positional filtering differences** as the remaining likely sources. Intronic mismatch definitions also diverge, but those are expected given STAR’s standard handling.

**Recommended next steps:**
1) Inspect conversion position filters: test whether ignoring first/last N positions reduces NTR deltas for the single‑mapping genes.
2) Compare STAR vs GEDI conversion counts on the same BAM reads for a few discordant genes (manual spot checks).
3) Decide whether to add a compatibility mode for conversion position filtering, if justified.

---

## Remaining Hypotheses (from instrumentation)

1) **Read‑assignment/compatibility differences:** Read‑level tracing shows STAR `PASS` at the same loci where GEDI assigns higher weights, suggesting different transcript compatibility logic or counting granularity.
2) **Positional filtering bias:** Strong early/late conversion enrichment in discordant genes implies position‑dependent handling (implicit trimming, quality filters, or conversion masking).
3) **Intronic definition mismatch:** GEDI’s intronic/IntronicSense categories diverge from STAR’s standard annotation logic; this affects mismatch summaries even when exonic looks closer.
4) **Conversion counting details:** Differences in how conversions are attributed within overlapping or multi‑part alignments could still be driving residual deltas.

**Report Updated:** January 9, 2026  
**Analysis Tool:** STAR-SLAM Debug Instrumentation  
**Reference:** GRAND-SLAM (GEDI) v1.0
