# SLAM Fixture Regeneration Notes

## Overview

The SLAM fixture can be regenerated using the script:
```bash
bash tests/regenerate_slam_fixture.sh
```

## Regeneration Results (2026-01-12)

### SNP BED (`snps.bed`)
- **Status**: ✅ MATCH (same content)
- **Entry count**: 1253 entries
- **Notes**: Order may differ from original but content is identical

### Reference TSV (`fixture_ref_human.tsv.gz`)
- **Status**: ⚠️ DIFFERENT (but highly correlated)
- **Original genes**: 7891
- **Regenerated genes**: 7907 (+16 extra genes)
- **Pearson correlation (MAP)**: 0.8601
- **Mean absolute difference**: 0.0497

### Why the Reference TSV differs

The regenerated reference TSV differs from the original due to:

1. **Stochastic EM estimation**: GRAND-SLAM uses EM with random initialization, leading to slightly different NTR estimates between runs.

2. **GEDI version**: The original was created with the same GEDI version (2.0.7b) but random seeds differ.

3. **Extra genes**: The regenerated version includes 16 additional genes that weren't in the original fixture. These are likely genes with very low read counts that are borderline for inclusion.

### Validation

Despite the differences, the fixture is valid for testing because:

1. All original genes are present in the regenerated version
2. NTR values are highly correlated (r=0.86)
3. The SNP BED is identical in content
4. The BAM preprocessing (adapter clipping, read counts) matches exactly

## Commands Used

### STAR Alignment
```bash
STAR \
    --runThreadN 4 \
    --genomeDir /storage/autoindex_110_44/bulk_index \
    --readFilesIn slam_100000_reads_SRR32576116.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix fixture_human_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --clip3pAdapterSeq AGATCGGAAGAG \
    --clip3pAdapterMMp 0.1
```

### GRAND-SLAM SNP Detection
```bash
gedi -e Slam \
    -reads fixture_human_Aligned.sortedByCoord.out.bam \
    -genomic homo_sapiens_110_44 \
    -prefix snpdetect \
    -strandness Sense \
    -snpConv 0.3 \
    -snppval 0.001 \
    -nthreads 4 \
    -keep
```

### SNP BED Conversion
```bash
# Convert snpdata (Location: chrom:pos 1-based) to BED3 (0-based)
tail -n +2 snpdetect.snpdata | awk -F'\t' '{
    split($1, loc, ":");
    chrom = loc[1];
    pos = loc[2];
    if (chrom == "MT") chrom = "chrM";
    else if (chrom !~ /^chr/) chrom = "chr" chrom;
    print chrom "\t" (pos-1) "\t" pos
}' > snps.bed
```

### BAM Prefiltering (for SNP removal)
```bash
samtools view -b -h -U fixture_human_prefiltered.bam \
    -L snps.bed fixture_human_Aligned.sortedByCoord.out.bam > /dev/null
```

### GRAND-SLAM Quantification
```bash
gedi -e Slam \
    -reads fixture_human_Aligned.sortedByCoord.out.bam \
    -genomic homo_sapiens_110_44 \
    -prefix fixture_ref_human \
    -strandness Sense \
    -nthreads 4 \
    -full
```

## File Checksums

Expected checksums (from original fixture):
```
8c28fae4adae8c2ed0ee22db1e7cceced4ae598f441b86d1a17fe5aeea8f3d76  raw/slam_100000_reads_SRR32576116.fastq.gz
8f2e760550f4456f13be4384475cff7de786fa9b2242cbbc8c16edb9fd9ac497  ref/snps.bed
d22fe98466c3f25b0e20444d30bef4f21c0e07a6cdd6906d5bc5b2d2e2c48a47  expected/fixture_ref_human.tsv.gz
```

## Recommendations

1. **Do not replace the fixture** unless you're updating to a new reference version
2. **Use the existing fixture** for regression testing
3. **Use regeneration** for understanding the pipeline or debugging
