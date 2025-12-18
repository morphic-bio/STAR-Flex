# nf-core RNA-seq Real Dataset Integration Fixture

This directory contains a real-world dataset for integration testing of the trimming functionality.

## Source

This dataset is derived from GSE110004, SRR6357070:
- GEO Accession: GSE110004
- SRA Run: SRR6357070
- Source FASTQs: `/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz` and `SRR6357070_2.fastq.gz`

## Regenerating Expected Outputs

To regenerate the expected outputs with Trim Galore (cutadapt backend):

```bash
cd /mnt/pikachu/STAR-Flex/test/integration/trim/nfcore_real

# Ensure input files are present (gunzipped)
gunzip -c /mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz > input_R1.fastq
gunzip -c /mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz > input_R2.fastq

# Run Trim Galore
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

# Rename outputs
mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq

# Clean up intermediate files
rm -f *_report.txt *_trimmed.fq
```

## Running Integration Tests

To test trimming on this dataset:

```bash
cd tools/trimvalidate
./trimvalidate -1 ../../test/integration/trim/nfcore_real/input_R1.fastq \
               -2 ../../test/integration/trim/nfcore_real/input_R2.fastq \
               -o1 /tmp/out_R1.fastq \
               -o2 /tmp/out_R2.fastq \
               --quality 20 --length 20

# Compare with expected output
diff /tmp/out_R1.fastq ../../test/integration/trim/nfcore_real/expected_R1.fastq
diff /tmp/out_R2.fastq ../../test/integration/trim/nfcore_real/expected_R2.fastq
```

## Notes

- This is a real dataset from a published study (GSE110004)
- Contains 50,000 read pairs
- Expected outputs generated with Trim Galore 0.6.10 (cutadapt 5.1 backend)
- Test results are stored in `results/` directory after running the parity test suite
