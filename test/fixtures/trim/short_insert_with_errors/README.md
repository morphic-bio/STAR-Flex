# short_insert_with_errors

## Description

Overlap + adapters but with ≤10% mismatches.

## Expected Behavior

Trim correctly despite mismatches.

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/short_insert_with_errors
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
