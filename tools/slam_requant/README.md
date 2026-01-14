# slam_requant

Re-quantify SLAM outputs from a STAR dump without re-alignment.

## Build
```
make
```

## Usage (minimal)
```
./slam_requant \
  --dump dump.bin \
  --out out_prefix \
  --slamSnpMaskIn mask.bed.gz
```

## Usage (BAM + GTF + FASTA)
```
./slam_requant \
  --bam input.bam \
  --gtf genes.gtf \
  --fasta genome.fa \
  --out out_prefix \
  --slamSnpMaskIn mask.bed.gz
```

## Key options
- `--dump <path>`: STAR dump created with `--slamDumpBinary`.
- `--bam <path>`: BAM input (alternative to `--dump`).
- `--gtf <path>`: GTF annotation (required with `--bam`).
- `--fasta <path>`: Reference FASTA (required with `--bam`).
- `--out <prefix>`: output prefix for `SlamQuant.out` and diagnostics.
- `--slamSnpMaskIn <bed.gz>`: apply SNP mask during replay (optional).
- `--snpMaskFromBam`: build a simple SNP mask from the BAM if no bed is provided.
- `--snpMaskOut <path>`: optional output BED for the auto-built mask.
- `--trim5p/--trim3p`: manual trims (optional).
- `--autoTrim variance`: compute trims from dump using variance method.
- `--trimScope first|per-file`: auto-trim scope (default: first).
- `--strandness none|sense|antisense`: drop opposite strand reads (default: none).
- `--slamQcReport <prefix>`: write QC JSON + HTML (optional).

## Outputs
- `<out>.SlamQuant.out`
- `<out>.SlamQuant.out.diagnostics`
- `<out>.SlamQuant.out.transitions.tsv`
- `<out>.SlamQuant.out.mismatches.tsv`
- `<out>.SlamQuant.out.mismatchdetails.tsv`
- QC JSON/HTML if `--slamQcReport` is provided
