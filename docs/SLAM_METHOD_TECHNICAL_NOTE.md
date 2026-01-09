# SLAM Method Technical Note: Metrics and Benchmarks

This note captures the key SLAM metrics used for interpretation and testing, and
the published inter-method correlation ranges that provide a sanity check for
cross-tool comparisons. It is intended to inform STAR-SLAM implementation and
parity testing against GRAND-SLAM.

## Core GRAND-SLAM outputs (GEDI wiki)

Source: https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM

- Per-gene read count (fractional for multimappers).
- New-to-total RNA ratio (NTR): reported as MAP and mean, with 0.05/0.95
  quantiles and Beta posterior parameters (alpha, beta).
- Mismatch summaries for exonic/intronic and sense/antisense, single vs double,
  and position-specific profiles to diagnose trimming artifacts.
- Binomial mismatch rates (p_c) for new vs old RNA used by the EM model.

Preprocessing guidance emphasizes:
- Mapped reads must include `MD` (and preferably `NH`) tags.
- Avoid soft-clipping; use end-to-end mapping and trim 5' and 3' to remove
  read-end artifacts detected via mismatch-position plots.

## Published inter-method comparison benchmarks

The GRAND-SLAM paper itself does not publish explicit "acceptable" correlation
thresholds. A relevant external benchmark is the comparison study:

- Boileau et al. 2021, "A comparison of metabolic labeling and statistical
  methods to infer genome-wide dynamics of RNA turnover." (Briefings in
  Bioinformatics)
- Analysis scripts and tables are in:
  https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling

This study compares decay-rate estimates across protocols and methods (pulseR
and GRAND-SLAM). Pearson correlations are used, along with residual standard
error and RMSD, and confidence-interval overlap.

Using the published `paper/tables/tbl.xlsx` from that repo, Pearson r between
pulseR and GRAND-SLAM decay-rate estimates across 11,603 genes is:

- Timepoints 0-1-2h:
  - SLAM: r = 0.7722
  - TLS:  r = 0.8076
  - TUC:  r = 0.8360
- Timepoints 0-2-4h:
  - SLAM: r = 0.8629
  - TLS:  r = 0.8818
  - TUC:  r = 0.8776
- Timepoints 0-4-8h:
  - SLAM: r = 0.7810
  - TLS:  r = 0.8138
  - TUC:  r = 0.8104
- Timepoints 0-1-2-4-8h:
  - SLAM: r = 0.8856
  - TLS:  r = 0.8972
  - TUC:  r = 0.8938

Note: these correlations are for decay-rate estimates (not NTR). They are
useful as a sanity range for cross-method agreement but are not a direct
acceptance threshold for STAR-SLAM vs GRAND-SLAM parity.

## Implications for STAR-SLAM testing

- For direct parity, compare NTR and conversion-related metrics on a common
  gene set with adequate coverage (e.g., nT >= 20-50).
- Expect high correlation for well-covered genes when preprocessing is matched
  (annotation, SNP mask, trimming, mapping mode).
- Systematic deviations are often concentrated in mitochondrial and ribosomal
  genes or low-complexity regions; document these explicitly.

## Fixture parity benchmarks (STAR-SLAM vs GRAND-SLAM)

These values come from the 100k-read human fixture and the GRAND-SLAM reference
`fixture_ref_human.tsv.gz`, using `tests/slam/compare_fixture.py` with
read-count thresholds (>=20, >=50, >=100). They provide a concrete baseline for
expected parity in this repository.

### SNP mask from BED (slamSnpBed)

This is the most direct parity mode because the reference fixture was generated
using the same SNP mask. Use a BED only when a sample-specific VCF/BED is known;
otherwise prefer internal SNP detection for production data.

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20      | readcount | 384 | 0.998865 | 0.981026 | 0.999865 | 0.999489 |
| >=50      | readcount | 76  | 0.996069 | 0.986681 | 0.998412 | 0.996712 |
| >=100     | readcount | 23  | 0.994379 | 0.994066 | 0.996847 | 0.994071 |

### Internal SNP detection (slamSnpDetect 1)

This mode is expected to diverge slightly from the BED-based reference (different
SNP calling), but still provides strong agreement. This is the recommended
default when no external SNP set is available.

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20      | readcount | 384 | 0.989764 | 0.967831 | 0.991470 | 0.990457 |
| >=50      | readcount | 76  | 0.996314 | 0.980295 | 0.998116 | 0.994169 |
| >=100     | readcount | 23  | 0.993684 | 0.991098 | 0.996284 | 0.988142 |

Note: k/nT correlations consistently reach the ~0.999 level under the BED mask,
which reflects very tight agreement in the raw conversion signal. NTR remains
the primary biological endpoint, with correlation values interpreted in the
context of coverage thresholds and SNP handling mode.
