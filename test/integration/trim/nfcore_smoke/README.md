# nf-core RNA-seq Smoke Test Integration Fixture

This directory contains a small real-world dataset for integration testing of the trimming functionality.

## Source

This dataset is derived from nf-core/rnaseq test data. The original test data can be found at:
- Repository: https://github.com/nf-core/test-datasets
- Path: `datasets/rnaseq/testdata/`

## Downloading Original Test Data

To regenerate this fixture with the original nf-core test data:

1. Clone the nf-core test-datasets repository:
   ```bash
   git clone --depth 1 https://github.com/nf-core/test-datasets.git
   cd test-datasets/datasets/rnaseq/testdata
   ```

2. Identify the smallest paired-end FASTQ files (typically named `*_R1.fastq.gz` and `*_R2.fastq.gz`)

3. Copy or symlink the files to this directory:
   ```bash
   cp <smallest_R1>.fastq.gz test/integration/trim/nfcore_smoke/input_R1.fastq.gz
   cp <smallest_R2>.fastq.gz test/integration/trim/nfcore_smoke/input_R2.fastq.gz
   gunzip test/integration/trim/nfcore_smoke/input_R*.fastq.gz
   ```

4. Generate expected outputs using Trim Galore:
   ```bash
   cd test/integration/trim/nfcore_smoke
   trim_galore --paired --quality 20 --length 20 \
               --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
               --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
               input_R1.fastq input_R2.fastq

   mv input_R1_val_1.fq expected_R1.fastq
   mv input_R2_val_2.fq expected_R2.fastq
   ```

## Current Dataset

**Note**: The current `input_R*.fastq` files in this directory are **synthetic** - they are derived from (but not identical to) the nf-core/rnaseq test dataset. They are designed to mimic real-world characteristics:
- Realistic read lengths (50-100 bp)
- Mixed quality scores
- Some adapter contamination
- Some low-quality tails

These synthetic files can be used for immediate testing without requiring external downloads. For testing with the actual nf-core dataset, follow the "Downloading Original Test Data" instructions above.

## Running Integration Tests

To test trimming on this dataset:

```bash
cd tools/trimvalidate
./trimvalidate -1 ../../test/integration/trim/nfcore_smoke/input_R1.fastq \
               -2 ../../test/integration/trim/nfcore_smoke/input_R2.fastq \
               -o1 /tmp/out_R1.fastq \
               -o2 /tmp/out_R2.fastq \
               --quality 20 --length 20

# Compare with expected output
diff /tmp/out_R1.fastq ../../test/integration/trim/nfcore_smoke/expected_R1.fastq
diff /tmp/out_R2.fastq ../../test/integration/trim/nfcore_smoke/expected_R2.fastq
```

## Notes

- This is a "smoke test" dataset - small enough for quick validation but realistic enough to catch integration issues
- For comprehensive testing, use the full nf-core test dataset
- Expected outputs should be regenerated whenever Trim Galore behavior changes or when updating to match newer versions
