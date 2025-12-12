# Flex helper scripts

This folder holds the minimal helpers to build a flex-capable STAR reference:

- `make_flex_reference.sh`: Create probe-only pseudo-chromosome FASTA/GTF from the 10x Flex probe CSV (one contig per probe). Outputs `<prefix>.fa` and `<prefix>.gtf`.
- `make_flex_star_index.sh`: End-to-end wrapper that (1) runs `make_flex_reference.sh`, (2) concatenates the base genome FASTA/GTF with the probe FASTA/GTF, and (3) runs `STAR --runMode genomeGenerate` on the combined reference. Outputs combined FASTA/GTF and the STAR index under `star_index/`.
- `filter_probes.sh`: Drop probes with `DEPRECATED` in the probe_id and emit a manifest.
- `build_filtered_reference.sh`: Full pipeline to filter probes, generate probe-only FASTA/GTF, build a hybrid FASTA/GTF (base genome + filtered probes), and emit BED/probe lists + manifests under `filtered_reference/`.
- `make_filtered_star_index.sh`: Run STAR genomeGenerate on the filtered hybrid reference produced by `build_filtered_reference.sh`.
- `mk_probe_bed.py`: Helper used by `build_filtered_reference.sh` to generate exon BED entries for probe genes.

Usage examples:

```bash
# Build probe-only FASTA/GTF
./scripts/make_flex_reference.sh probes.csv work/probes_only

# Build combined reference + STAR index
./scripts/make_flex_star_index.sh \
  --probe-csv probes.csv \
  --base-fasta /path/to/genome.fa \
  --base-gtf /path/to/genes.gtf \
  --out-dir work/flex_reference \
  --threads 24

# Build filtered hybrid reference (drops DEPRECATED probes) and STAR index
./scripts/build_filtered_reference.sh \
  --probe-set probes.csv \
  --base-fasta /path/to/genome.fa \
  --base-gtf /path/to/genes.gtf \
  --work-dir work
./scripts/make_filtered_star_index.sh \
  --filtered-reference work/filtered_reference \
  --output-dir work/star_index_filtered \
  --threads 24
```

Notes:
- Requires an available `STAR` binary in `PATH` (or pass `--star-bin`).
- `build_filtered_reference.sh` depends on Python 3 with pandas (`mk_probe_bed.py`) and uses `filter_probes.sh`/`make_flex_reference.sh` internally.
- Defaults are conservative; adjust `--threads`, `--sa-index-n-bases`, `--chr-bin-n-bits`, and `--genome-limit` as needed for your hardware/reference.
- Scripts are self-contained and do not depend on the legacy reference tooling.
