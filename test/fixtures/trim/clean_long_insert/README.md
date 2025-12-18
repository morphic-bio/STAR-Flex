# clean_long_insert

## Description

Insert length >> read length; no adapters, high quality.

## Expected Behavior

No trimming expected, reads unchanged.

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/clean_long_insert
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
