# Salmon Equivalence Class Fixture Status

## Current Status: ✅ Complete (Bias-Disabled for v1 Parity)

The Salmon equivalence class fixture has been **regenerated with bias correction disabled** for v1 parity testing with `em_quant`.

## What Was Created

✅ **Transcriptome FASTA** (`chr22_trans.fa`): 70 transcripts from human chr22 slice (23.8-23.98 Mb), ~105KB.

✅ **Quantification File** (`quant.sf`): Created with `--noLengthCorrection --noFragLengthDist --noEffectiveLengthCorrection` flags.
- EffectiveLength = 100 for all transcripts (bias disabled)
- 70 transcripts, ~3.3KB

✅ **Equivalence Classes** (`eq_classes.txt`): 188 ECs from ~5000 synthetic read pairs, ~4.2KB.

## Parity Test Results

With `em_quant` v1 (20% tolerance, 3.0 near-zero threshold):
- ✓ All 55 significant transcripts within tolerance
- ℹ 15 transcripts have near-zero differences (Salmon=0, em_quant<3.0) - expected behavior due to different regularization

## Generation Process

1. **Transcriptome Creation**: Used a Python script to extract transcripts directly from the GTF and FASTA (avoiding `gffread` coordinate issues with slice boundaries).

2. **Synthetic Read Generation**: Generated synthetic paired-end reads from the transcriptome using `art_illumina` (10x coverage, 100bp reads, Illumina HiSeq 2500 error profile).

3. **Salmon Quantification**: Ran Salmon quant with `--dumpEq` on synthetic reads, achieving 99.96% mapping rate and successfully generating equivalence classes.

## Generation Method

The transcriptome was created using a Python script that:
1. Reads the chr22 FASTA slice
2. Filters the human GTF to chr22 features within the 23.8-23.98 Mb window
3. Extracts exon sequences and concatenates them per transcript
4. Handles reverse-complement for minus-strand transcripts
5. Writes a FASTA file with transcript sequences

## Regeneration

To regenerate the fixture with bias disabled (required for v1 parity testing):

1. **Generate synthetic reads** (if not already present):
   ```bash
   art_illumina -ss HS25 -i chr22_trans.fa -l 100 -f 10 -m 200 -s 10 -o synthetic_reads
   ```

2. **Build Salmon index**:
   ```bash
   docker run --rm -v $(pwd):/data combinelab/salmon:latest \
       salmon index -t /data/chr22_trans.fa -i /data/salmon_idx
   ```

3. **Run Salmon quant with bias disabled**:
   ```bash
   docker run --rm -v $(pwd):/data combinelab/salmon:latest \
       salmon quant -i /data/salmon_idx -l A \
       -1 /data/synthetic_reads1.fq -2 /data/synthetic_reads2.fq \
       --dumpEq --threads 4 \
       --noLengthCorrection --noFragLengthDist --noEffectiveLengthCorrection \
       -o /data/salmon_output
   ```

4. **Copy artifacts**:
   ```bash
   gunzip -c salmon_output/aux_info/eq_classes.txt.gz > eq_classes.txt  # Salmon v1.10+ uses gzip
   cp salmon_output/quant.sf .
   ```

## Tools Status

The generation process uses:
- **Python 3**: For transcript extraction (built-in)
- **Salmon**: For indexing and quantification (via Docker: `continuumio/miniconda3`)

## Files Included

- `chr22_trans.fa`: Transcriptome FASTA (71 transcripts)
- `salmon_idx/`: Salmon index directory
- `eq_classes.txt`: Equivalence classes file (generated from synthetic reads)
- `quant.sf`: Quantification results
- `synthetic_reads1.fq` / `synthetic_reads2.fq`: Synthetic paired-end reads (optional, can be regenerated)
