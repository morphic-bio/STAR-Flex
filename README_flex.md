# STAR-Flex: Inline Hash Pipeline for 10x Flex Samples

This document describes the flex pipeline integration in this STAR fork.

## Overview

STAR-Flex extends STAR with an inline hash-based pipeline optimized for 10x Genomics Flex (Fixed RNA Profiling) samples using probes for transcript detection and RTL tags for multiplexing. We generate a hybrid reference with the regular genome and with synthetic chromosomes for each of the probes. This allows to use the STAR alignment routines to quantify probe alignment and use the genomic hits to confirm the match and detect off-probe noise. However, the rest of the workflow diverges from the standard STAR solo workflow, largely due to the presence of RTL tags for multiplexing samples. Because these are on the same mate as the probe and not the cell barcode, STAR's barcode and UMI correction, and UMI deduping routines could not be used. Furthermore, the noise characteristics of Flex are different that the native STAR's multimapping ad emptyDrops functions could not be used. A fast inline path was created to handle Flex processing after STAR alignment.

1. **Sample tag detection** during alignment identifies multiplexed sample barcodes
2. **Inline hash capture** stores CB/UMI/gene tuples directly in memory
3. **Cell Barcode (CB) correction** applies 1MM pseudocount-based correction (Cell Ranger compatible)
4. **UMI correction** uses clique-based 1MM deduplication
5. **Cell filtering** via OrdMag (simple EmptyDrops) or full EmptyDrops per sample
6. **Tag occupancy filtering** via Monte Carlo estimation of the expected distribution of samples per cell barcode
7. **MEX output** produces raw and per-sample filtered matrices

When `--flex no` (default), STAR behavior is identical to upstream.

For detailed technical documentation of the data flow and algorithms, see [docs/flex_methodology.md](docs/flex_methodology.md).

## Quick Start

```bash
STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist sample_whitelist.tsv \
  --soloProbeList probe_list.txt \
  --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexOutputPrefix output/per_sample \
  --soloMultiMappers Rescue \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloFeatures Gene \
  --outFileNamePrefix output/
```

## Required Inputs

| Input | Description |
|-------|-------------|
| Flex reference genome | Hybrid genome with probe pseudo-chromosomes (see [Building References](#building-references)) |
| CB whitelist | 10x barcode whitelist (e.g., `737K-fixed-rna-profiling.txt`) |
| Sample whitelist | TSV mapping sample tag sequences to labels |
| Probe list | Gene list from probe set |
| Sample probe barcodes | 10x probe barcode sequences file |

## Parameters

### Master Switch

| Flag | Default | Description |
|------|---------|-------------|
| `--flex` | `no` | Enable flex pipeline (`yes`/`no`) |

### Sample Detection

| Flag | Default | Description |
|------|---------|-------------|
| `--soloSampleWhitelist` | - | Path to sample tag whitelist TSV |
| `--soloProbeList` | auto | Path to probe gene list (auto-detects from genome index if not specified) |
| `--soloSampleProbes` | - | Path to 10x sample probe barcodes |
| `--soloSampleProbeOffset` | 0 | Offset in read for sample probe sequence |
| `--soloSampleSearchNearby` | `yes` | Search nearby positions for sample tag |
| `--soloSampleStrictMatch` | `no` | Require strict match for sample tag |

### FlexFilter (Cell Calling)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloFlexExpectedCellsPerTag` | 0 | Expected cells per sample tag |
| `--soloFlexExpectedCellsTotal` | 0 | Total expected cells (alternative to per-tag) |
| `--soloFlexAllowedTags` | - | Optional: restrict to specific sample tags |
| `--soloFlexOutputPrefix` | - | Output prefix for per-sample MEX |

### EmptyDrops Parameters (Advanced)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloFlexEdNiters` | 10000 | Monte Carlo simulation iterations |
| `--soloFlexEdFdrThreshold` | 0 (disabled) | FDR threshold for cell calling; if set (>0), FDR gate is used |
| `--soloFlexEdPvalueThreshold` | 0.05 | Raw p-value threshold when FDR gate is disabled (default behavior) |
| `--soloFlexEdLower` | 100 | Lower UMI bound for ambient profile |

## Output Structure

```
output/
├── Solo.out/Gene/raw/          # Raw MEX (all barcodes)
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── matrix.mtx
└── per_sample/                  # Per-sample filtered MEX (labels from whitelist)
    ├── SampleA/Gene/filtered/
    ├── SampleB/Gene/filtered/
    └── flexfilter_summary.tsv   # Cell calling statistics
```

## Building References

The flex pipeline requires a hybrid reference genome that includes pseudo-chromosomes for probe sequences. Scripts are provided in `scripts/` to build these references:

### Integrated Index Generation (Recommended)

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/flex_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf \
  --sjdbOverhang 100 \
  --flexGeneProbeSet /path/to/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv \
  --runThreadN 8
```

#### Required Inputs

| Input | Description |
|-------|-------------|
| `--genomeFastaFiles` | Base genome FASTA file |
| `--sjdbGTFfile` | Gene annotation GTF file (can be gzipped) |
| `--flexGeneProbeSet` | 10x Flex probe CSV file (50bp gene probes) |

#### Flex Index Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--flexGeneProbeSet` | - | Path to 50bp gene probe CSV file |
| `--flexGeneProbeLength` | 50 | Expected probe length (fails if mismatch) |

#### Output Structure

```
flex_index/
├── probe_gene_list.txt           # Unique gene IDs with probes (auto-detected for --soloProbeList)
├── flex_probe_artifacts/         # Probe processing artifacts
│   ├── filtered_probe_set.csv    # Probes matching GTF genes
│   ├── probes_only.fa            # Probe-only FASTA
│   ├── probes_only.gtf           # Probe-only GTF entries
│   ├── genome.filtered.fa        # Hybrid FASTA (used for indexing)
│   ├── genes.filtered.gtf        # Hybrid GTF (used for indexing)
│   ├── probe_genes_exons.bed     # Probe coordinates
│   ├── probe_list.txt            # Unique gene IDs
│   └── metadata/
│       └── reference_manifest.json
├── Genome                        # Standard STAR index files
├── SA
├── SAindex
└── ... (other STAR index files)
```

#### Probe Filtering Rules

The integrated preprocessor applies these filters:
1. **50bp A/C/G/T only** - Fails if any probe has invalid length or characters
2. **Skip DEPRECATED** - Excludes probes marked as deprecated
3. **Gene match** - Keeps only probes whose gene_id exists in the target GTF
4. **Deterministic ordering** - Stable sort by gene_id then probe_id

### Alternative: Shell Scripts

For custom workflows or debugging, standalone shell scripts are available:

```bash
# Filter probes and build hybrid reference
./scripts/filter_probes_to_gtf.sh \
  --probe-set /path/to/probes.csv \
  --gtf /path/to/genes.gtf.gz \
  --base-fasta /path/to/genome.fa \
  --output-dir ./probe_artifacts
```

The legacy `build_filtered_reference.sh` and `make_filtered_star_index.sh` scripts are also available. See [scripts/README.md](scripts/README.md) for details.

### Using the Flex Index

After building, use the index with the probe gene list:

```bash
STAR \
  --genomeDir /path/to/flex_index \
  --flex yes \
  ... # other flex parameters
  # --soloProbeList is auto-detected from probe_gene_list.txt in the index directory
```

## Standalone FlexFilter Tool

A standalone tool `run_flexfilter_mex` is available for offline MEX processing. This allows re-running the OrdMag/EmptyDrops cell calling pipeline on existing composite MEX files without re-running STAR alignment.

**Use cases:**
- Parameter tuning (adjust expected cells, EmptyDrops thresholds)
- Reprocessing with different filtering settings
- Integration with non-STAR pipelines (any tool producing composite CB+TAG MEX)
- Batch reprocessing of archived STAR outputs

### Building

The tool is optional and not built by the default `make STAR` target:

```bash
cd source
make flexfilter
```

This produces `tools/flexfilter/run_flexfilter_mex`.

### Input Requirements

The tool expects a composite MEX directory containing:
- `matrix.mtx` - Matrix Market sparse matrix (or `InlineHashDedup_matrix.mtx`)
- `barcodes.tsv` - Composite barcodes in CB16+TAG8 format (24 characters)
- `features.tsv` - Gene IDs (tab-separated)

The composite barcode format concatenates the 16bp cell barcode with the 8bp sample tag:
```
AAACCCAAGAAACACTACGTACGT  # CB16 (AAACCCAAGAAACACT) + TAG8 (ACGTACGT)
```

### Basic Usage

```bash
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /path/to/Solo.out/Gene/raw \
  --total-expected 12000 \
  --output-prefix /path/to/filtered_output
```

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--mex-dir` | Path to composite MEX directory (required) |
| `--total-expected` | Total expected cells across all samples (required) |
| `--output-prefix` | Output directory prefix (required) |
| `--sample-whitelist` | TSV file mapping sample names to tag sequences |
| `--ed-lower-bound` | Lower UMI bound for EmptyDrops (default: 500) |
| `--ed-fdr` | FDR threshold for EmptyDrops (default: 0.01) |
| `--disable-occupancy` | Skip occupancy post-filter (for testing) |

### Output Structure

```
output_prefix/
├── SampleA/Gene/filtered/
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── EmptyDrops/
│       └── emptydrops_results.tsv
├── SampleB/Gene/filtered/
│   └── ...
└── flexfilter_summary.tsv
```

### Example Workflow: Reprocess with Different Expected Cells

```bash
# Original STAR run produced Solo.out/Gene/raw/
# Reprocess with higher cell expectation
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /storage/run1/Solo.out/Gene/raw \
  --total-expected 20000 \
  --output-prefix /storage/run1/refiltered_20k

# Or with explicit sample whitelist
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /storage/run1/Solo.out/Gene/raw \
  --sample-whitelist samples.tsv \
  --total-expected 15000 \
  --output-prefix /storage/run1/refiltered_explicit
```

Sample whitelist format (`samples.tsv`):
```
Sample_A	ACGTACGT
Sample_B	TGCATGCA
Sample_C	GGCCGGCC
```
Labels in the first column are used verbatim for per-sample directories, and the order in the whitelist is preserved.

### Testing

```bash
# Requires tests/gold_standard/ fixtures
./tools/flexfilter/test_smoke.sh

# Validate output format
./tools/flexfilter/validate_output.py /path/to/output
```

See [tools/flexfilter/README.md](tools/flexfilter/README.md) for complete CLI reference and advanced options.

## Building STAR-Flex

Standard STAR build process:

```bash
cd source
make -j8
```

The flex objects are automatically included in the build.

## Testing

See [docs/TESTING_flex.md](docs/TESTING_flex.md) for detailed testing instructions.

Quick test:

```bash
./tests/run_flex_multisample_test.sh
```

Gold standard comparison files are bundled in `tests/gold_standard/`.

## Code Organization

```
source/
├── libflex/                    # Core flex filtering library
│   ├── FlexFilter.cpp/h        # Main filter orchestration
│   ├── EmptyDropsMultinomial.cpp/h  # Full EmptyDrops
│   ├── OrdMagStage.cpp/h       # Simple EmptyDrops (OrdMag)
│   └── OccupancyGuard.cpp/h    # Occupancy-based filtering
├── solo/
│   └── CbCorrector.cpp/h       # CB correction with pseudocounts
├── SampleDetector.cpp/h        # Sample tag detection
├── InlineCBCorrection.cpp/h    # Inline hash CB correction
├── UMICorrector.cpp/h          # Clique-based UMI correction
├── MexWriter.cpp/h             # MEX matrix output
├── GeneResolver.cpp/h          # Probe-to-gene mapping
├── SoloFeature_flexfilter.cpp  # FlexFilter integration
├── SoloFeature_writeMexFromInlineHashDedup.cpp
└── UmiCodec.h                  # UMI encoding/decoding helpers

tools/flexfilter/               # Standalone FlexFilter CLI
├── run_flexfilter_mex.cpp      # Main CLI wrapper
├── Makefile                    # Build configuration
├── README.md                   # CLI documentation
├── test_smoke.sh               # Smoke test script
└── validate_output.py          # Output validation script
```

## Compatibility

- Baseline: STAR 2.7.11b
- When `--flex no` (default), behavior is identical to upstream STAR
- Upstream `README.md` and `CHANGES.md` are not modified
