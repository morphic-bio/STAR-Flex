# STAR-Flex Parameter Reference

This document lists **STAR-Flex-only parameters** that are not present in upstream STAR. For all other parameters, see upstream `README.md` and `parametersDefault`.

## Bulk RNA-seq Features

### Cutadapt-Style Trimming

**Note**: Trimming is a general-purpose feature for bulk RNA-seq (and other workflows), not specific to the Flex pipeline.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--trimCutadapt` | `-` | Enable cutadapt-style trimming (`Yes`/`-`). When enabled, all `clip*` parameters are ignored and ClipMate clipping is bypassed. |
| `--trimCutadaptQuality` | `20` | Quality threshold for 3' trimming (Phred scale, integer) |
| `--trimCutadaptMinLength` | `20` | Minimum read length after trimming. Pairs with either mate shorter than this are dropped (integer) |
| `--trimCutadaptAdapter` | `-` | Custom adapter sequences for R1 and R2 (space-separated). Default (`-`) uses TruSeq adapters: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` (R1) and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` (R2) |

**Note**: Trimming achieves perfect parity with Trim Galore/cutadapt v5.1. See [docs/trimming.md](docs/trimming.md) for algorithm details and usage examples.

### Y-Chromosome BAM Split

**Note**: This feature was developed for **Morphic requirements for KOLF cell lines**. It is **not connected to the Flex pipeline** and is a general-purpose feature usable with any STAR workflow.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--emitNoYBAM` | `no` | Enable Y-chromosome BAM splitting (`yes`/`no`). When enabled, emits `<out>_Y.bam` and `<out>_noY.bam`. Primary BAM is suppressed unless `--keepBAM yes` is specified. |
| `--keepBAM` | `no` | Keep primary BAM output when `--emitNoYBAM yes` is enabled (`yes`/`no`) |
| `--noYOutput` | - | Optional: override default path for noY BAM output |
| `--YOutput` | - | Optional: override default path for Y BAM output |

See [docs/Y_CHROMOSOME_BAM_SPLIT.md](Y_CHROMOSOME_BAM_SPLIT.md) for details.

## Flex-Specific Features

### Flex Pipeline

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--flex` | `no` | Enable flex pipeline (`yes`/`no`) |

**Note**: Flex pipeline is specific to 10x Genomics Flex (Fixed RNA Profiling) samples. See [README_flex.md](../README_flex.md) for complete flex parameter documentation.

## Upstream Parameters

For all other parameters (genome generation, alignment, solo, etc.), refer to upstream STAR documentation:
- `README.md` - User guide
- `parametersDefault` - Complete parameter list with defaults
- `doc/STARmanual.pdf` - Full manual
