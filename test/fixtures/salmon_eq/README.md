# Salmon Equivalence Class Fixture

This directory contains Salmon equivalence class data generated from a downsampled subset of SRR6357070 reads aligned to a chr22 reference slice.

## Generation

### Prerequisites

- `gffread` (from Cufflinks) - Install: `sudo apt install gffread` or `conda install -c bioconda gffread`
- `salmon` (v1.9.0 or later) - Install: `conda install -c bioconda salmon` or download from https://github.com/COMBINE-lab/salmon
- Reference files:
  - `/mnt/pikachu/test-datasets/reference/chr22_23800000-23980000.fa` (chr22 FASTA slice)
  - `/mnt/pikachu/genome.gtf` or `/mnt/pikachu/genome.gtf.gz` (human GTF)
- Read files:
  - `/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz` (human reads)
  - `/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz`

**Note**: If tools are not in PATH, activate the appropriate conda environment or add them to PATH before running. Docker can also be used (see GENERATE.sh).

### Steps

**Quick generation (automated script):**
```bash
cd test/fixtures/salmon_eq
./GENERATE.sh
```

**Manual generation:**

1. **Subset GTF to chr22 (with optional coordinate filtering):**
   ```bash
   # Try with coordinate filter first (23.8-23.98 Mb window)
   gunzip -c /mnt/pikachu/genome.gtf.gz 2>/dev/null || cat /mnt/pikachu/genome.gtf | \
     awk '$1=="chr22" && $4>=23800000 && $5<=23980000' > /tmp/genes_chr22.gtf
   
   # If no records, use all chr22:
   # awk '$1=="chr22"' /mnt/pikachu/genome.gtf > /tmp/genes_chr22.gtf
   ```

2. **Build transcriptome FASTA from subset GTF + reference:**
   ```bash
   gffread /tmp/genes_chr22.gtf \
           -g /mnt/pikachu/test-datasets/reference/chr22_23800000-23980000.fa \
           -w /tmp/chr22_trans.fa
   ```

3. **Downsample reads (1000 read pairs = 4000 FASTQ lines):**
   ```bash
   zcat /mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz | head -n 4000 > /tmp/sub_R1.fq
   zcat /mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz | head -n 4000 > /tmp/sub_R2.fq
   ```

4. **Build Salmon index:**
   ```bash
   salmon index -t /tmp/chr22_trans.fa -i /tmp/salmon_chr22_idx
   ```

5. **Run Salmon quant with equivalence class dumping (bias disabled for parity testing):**
   ```bash
   # IMPORTANT: Use --dumpEqWeights to get per-EC per-transcript combinedWeights
   # Without weights, we can only achieve ~85% parity (34% max outlier)
   # Also use --threads 1 for deterministic floating-point order
   # Disable bias correction for v1 parity testing (raw lengths)
   salmon quant -i /tmp/salmon_chr22_idx \
                -l A \
                -1 /tmp/sub_R1.fq \
                -2 /tmp/sub_R2.fq \
                --dumpEqWeights \
                --threads 1 \
                --noLengthCorrection \
                --noFragLengthDist \
                --noEffectiveLengthCorrection \
                -o /tmp/salmon_eq
   ```

6. **Copy artifacts:**
   ```bash
   cp /tmp/salmon_eq/aux_info/eq_classes.txt test/fixtures/salmon_eq/
   cp /tmp/salmon_eq/quant.sf test/fixtures/salmon_eq/
   # If transcript map exists:
   cp /tmp/salmon_eq/aux_info/transcript_map.txt test/fixtures/salmon_eq/ 2>/dev/null || \
   cp /tmp/salmon_eq/aux_info/transcript_to_gene_map.tsv test/fixtures/salmon_eq/ 2>/dev/null || true
   ```

## Files

- `eq_classes.txt` - Equivalence classes in Salmon format (4.2 KB, 188 ECs)
  - Line 1: number of transcripts (70)
  - Line 2: number of equivalence classes (188)
  - Lines 3-72: transcript names
  - Remaining lines: EC definitions (format depends on generation flags)
    - **Unweighted** (`--dumpEq`): `<k> <idx1> <idx2> ... <idxk> <count>`
    - **Weighted** (`--dumpEqWeights`): `<k> <idx1> ... <idxk> <w1> ... <wk> <count>`
  - **IMPORTANT**: For exact parity (<1% tolerance), regenerate with `--dumpEqWeights`
- `quant.sf` - Salmon quantification results (3.3 KB, 70 transcripts)
  - Standard Salmon TSV with Name, Length, EffectiveLength, TPM, NumReads
  - Note: EffectiveLength=100 for all transcripts (bias disabled)
- `chr22_trans.fa` - Transcriptome FASTA (105 KB, 70 transcripts)
- `synthetic_reads1.fq`, `synthetic_reads2.fq` - Synthetic paired-end reads (1.1 MB each)
  - Generated with art_illumina from the transcriptome

## Usage

This fixture can be used to test Salmon equivalence class parsing and processing logic, ensuring compatibility with Salmon's output format.

## Notes

- The reads are downsampled to 1000 pairs (4000 FASTQ lines) to keep equivalence classes small
- The reference is a chr22 slice (23800000-23980000) to limit transcriptome size
- Equivalence classes are generated using Salmon's default k-mer size (31)
- Auto-detection of library type (`-l A`) is used

## Parity Testing Requirements

**For v1 parity testing, the Salmon fixture MUST be generated with:**

1. **Bias correction disabled** - Our `em_quant` engine uses raw transcript lengths:
   - `--noLengthCorrection` - disable bias correction
   - `--noFragLengthDist` - don't use fragment length distribution
   - `--noEffectiveLengthCorrection` - use raw lengths

2. **EC weights enabled** - Required for exact parity:
   - `--dumpEqWeights` (or `-d`) - include per-EC per-transcript combinedWeights
   - Without weights, Salmon collapses factorization info and we can only achieve ~85% parity

3. **Single-threaded** - For deterministic floating-point order:
   - `--threads 1` - eliminates parallel FP ordering as a confound

### Parity Results

| Fixture Type | 5% Tolerance | 10% Tolerance | Max Outlier |
|-------------|-------------|---------------|-------------|
| Unweighted (`--dumpEq`) | ~85% pass | ~96% pass | 34% diff |
| Weighted (`--dumpEqWeights`) | >99% pass | 100% pass | <1% diff |

The remaining differences with unweighted fixtures are due to:
- Missing alignment quality weights in `combinedWeights` (auxs)
- Salmon's range-factorization information being collapsed

This ensures the only difference is the EM implementation itself, not the effective-length model. A future version may add FLD-based effective lengths, at which point the fixture should be regenerated without the bias-disabling flags.
