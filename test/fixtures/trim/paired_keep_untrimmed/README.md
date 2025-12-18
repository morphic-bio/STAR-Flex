# paired_keep_untrimmed

## Description

One mate needs trimming, other is clean.

## Expected Behavior

Trim affected mate, keep both if lengths ≥ 20.

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/paired_keep_untrimmed
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
