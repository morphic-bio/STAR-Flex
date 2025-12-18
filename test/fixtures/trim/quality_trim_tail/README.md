# quality_trim_tail

## Description

Low-quality tail (Q < 20) at 3' end before any adapter.

## Expected Behavior

Quality trimming removes tail.

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/quality_trim_tail
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
