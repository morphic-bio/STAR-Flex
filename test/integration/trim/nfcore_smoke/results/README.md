# Integration Test Results

This directory contains the results of running the nfcore_smoke integration test.

## Files

- `status.txt` - Test status (PASS/FAIL)
- `diff_R1.txt` - Unified diff between expected and actual R1 output
- `diff_R2.txt` - Unified diff between expected and actual R2 output

## Running the Test

The integration test is automatically run as part of the parity test suite:

```bash
# From project root
cd source && make test_trim_parity

# Or directly
cd tools/trimvalidate && make test
```

## Regenerating Results

Results are automatically regenerated each time the parity tests run. To manually regenerate:

```bash
cd tools/trimvalidate
./trimvalidate -1 ../../test/integration/trim/nfcore_smoke/input_R1.fastq \
               -2 ../../test/integration/trim/nfcore_smoke/input_R2.fastq \
               -o1 /tmp/out_R1.fastq -o2 /tmp/out_R2.fastq \
               --quality 20 --length 20

diff -u ../../test/integration/trim/nfcore_smoke/expected_R1.fastq /tmp/out_R1.fastq > ../../test/integration/trim/nfcore_smoke/results/diff_R1.txt
diff -u ../../test/integration/trim/nfcore_smoke/expected_R2.fastq /tmp/out_R2.fastq > ../../test/integration/trim/nfcore_smoke/results/diff_R2.txt
```

## Notes

- Empty diff files indicate perfect parity with Trim Galore output
- These results are stored alongside the parity test logs for verification
- The integration test uses synthetic data derived from nf-core test datasets
- For testing with actual nf-core data, see the parent directory README.md
