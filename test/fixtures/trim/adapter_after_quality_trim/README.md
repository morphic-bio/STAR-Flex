# adapter_after_quality_trim

## Description

Low-quality tail that, when trimmed, exposes adapter.

## Expected Behavior

Quality trim first, then adapter trim.

## Regenerating Expected Outputs

```bash
cd test/fixtures/trim/adapter_after_quality_trim
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
