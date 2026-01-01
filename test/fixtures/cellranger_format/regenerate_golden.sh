#!/bin/bash
# Regenerate golden output files from Perl script
# Run this when Perl script rules change

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PERL_SCRIPT="/mnt/pikachu/Bulk-RNA-seq/star_salmon_aws/widgets/star_salmon_aws/10x_format_fa_gtf/Dockerfiles/format-fa-gtf.pl"

if [ ! -f "$PERL_SCRIPT" ]; then
    echo "ERROR: Perl script not found: $PERL_SCRIPT"
    exit 1
fi

INPUTFA="$SCRIPT_DIR/test_genome.fa"
INPUTGTF="$SCRIPT_DIR/test_genes.gtf"
OUTPUTFA="$SCRIPT_DIR/expected_genome.fa"
OUTPUTGTF="$SCRIPT_DIR/expected_genes.gtf"

echo "Regenerating golden outputs..."
echo "Input FASTA: $INPUTFA"
echo "Input GTF: $INPUTGTF"
echo "Output FASTA: $OUTPUTFA"
echo "Output GTF: $OUTPUTGTF"

INPUTFA="$INPUTFA" OUTPUTFA="$OUTPUTFA" \
INPUTGTF="$INPUTGTF" OUTPUTGTF="$OUTPUTGTF" \
overwrite=1 \
perl "$PERL_SCRIPT"

echo "Golden outputs regenerated successfully!"

