# EM Unit Test Fixture

This directory contains a handcrafted synthetic equivalence class fixture for unit testing the EM algorithm.

## Format

- 3 transcripts: T1, T2, T3
- 3 equivalence classes:
  - EC1: 100 reads map uniquely to T1
  - EC2: 50 reads map uniquely to T2
  - EC3: 30 reads map ambiguously to T1 or T2

## Expected Results

With uniform initialization:
- T1: ~115 reads (100 unique + ~15 from ambiguous EC3)
- T2: ~65 reads (50 unique + ~15 from ambiguous EC3)
- T3: 0 reads

## Usage

```bash
cd ../../../tools/em_quant
./em_quant -e ../../test/fixtures/em_unit/eq_classes.txt -o /tmp/em_unit_output.tsv
```
