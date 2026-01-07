# AutoIndex + CellRanger-Style References

STAR-Flex adds an index-time workflow to:
- Download a pinned reference (or use provided inputs)
- Verify download integrity with `cksum` (when possible)
- Optionally format FASTA/GTF using CellRanger-style rules
- Build (or rebuild) the STAR genome index in `--genomeDir`

This is implemented inside `STAR --runMode genomeGenerate`.

## Quick Start (Download + Format + Index)

Build a CellRanger-style GRCh38 reference (default `--cellrangerRefRelease 2024-A`, i.e. Ensembl 110 + GENCODE 44):

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --autoIndex Yes \
  --cellrangerStyleIndex Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

Rebuild from scratch (re-download + re-format + re-index):

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --autoIndex Yes \
  --forceAllIndex Yes \
  --cellrangerStyleIndex Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

## Choosing a Default Release

AutoIndex supports a small set of named “CellRanger-style” releases for default URLs:

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --autoIndex Yes \
  --cellrangerRefRelease 2020-A \
  --cellrangerStyleIndex Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

Supported values are documented in `source/parametersDefault` under `cellrangerRefRelease`.

## Using Custom URLs

Provide both `--faUrl` and `--gtfUrl`:

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --faUrl  ftp://.../genome.fa.gz \
  --gtfUrl ftp://.../genes.gtf.gz \
  --cellrangerStyleIndex Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

If a URL is not in the trusted table and has no known checksum, STAR-Flex requires `--allUntrustedUrl Yes` (this disables integrity checking for those URLs).

## Using Local Files (No Download)

`--cellrangerStyleIndex` is independent of downloads; it can format local inputs too:

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf \
  --cellrangerStyleIndex Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

## Output Layout

When `--cellrangerStyleIndex Yes` is enabled, formatted inputs are written to:
- `${genomeDir}/cellranger_ref/genome.fa`
- `${genomeDir}/cellranger_ref/genes.gtf`

If `--genomeGenerateTranscriptome Yes` is also set, STAR-Flex writes:
- `${genomeDir}/transcriptome.fa`
- `${genomeDir}/cellranger_ref/transcriptome.fa`

Downloads are cached (default) under:
- `${genomeDir}/cellranger_ref_cache/`

Override with `--cellrangerStyleCacheDir /path/to/cache`.

## Integrity, Overwrites, and “Force” Flags

- `--forceIndex Yes`: rebuild index files even if the index exists (keeps cached downloads and formatted files).
- `--forceAllIndex Yes`: delete cached downloads and formatted files, then re-download/re-format/rebuild.
- `--autoCksumUpdate Yes`: for trusted URLs missing embedded checksums, attempt to retrieve checksums from CHECKSUMS files and cache them.
- `--replaceUnverifiableFiles Yes`: if a final decompressed file exists but cannot be verified, re-download/re-decompress and overwrite it.
  - Checksum mismatches are always hard errors (independent of this flag).

## Download-Only Mode

`--cellrangerStyleDownloadOnly Yes` exits after downloads complete (no formatting, no index build). This is useful if you want to stage the cache first:

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --faUrl  ftp://.../genome.fa.gz \
  --gtfUrl ftp://.../genes.gtf.gz \
  --cellrangerStyleCacheDir /path/to/cache \
  --cellrangerStyleDownloadOnly Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 1
```

If you use `--autoIndex Yes` to populate default URLs and the index already exists, STAR-Flex may skip the run; add `--forceIndex Yes` (or `--forceAllIndex Yes`) to force a download-only pass.

## Tests and Parity Checks

- Perl vs C++ formatter (byte-for-byte): `./test/run_cellranger_format_parity.sh`
- Perl fixtures vs STAR internal CellRanger path: `./test/run_cellranger_perl_star_parity.sh`
- AutoIndex smoke test (force flags + URL resolution): `./test/run_autoindex_smoke.sh`
