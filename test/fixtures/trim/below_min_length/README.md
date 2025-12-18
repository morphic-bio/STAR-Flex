# below_min_length

## Description

After trimming, reads are exactly at the minimum length threshold (20bp). Tests edge case where reads are just barely kept.

## Expected Behavior

Reads trimmed to exactly 20bp and kept (at minimum threshold).

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/below_min_length
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
