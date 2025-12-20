#!/usr/bin/env python3
"""
Build transcriptome FASTA from GTF and sliced reference FASTA.

This script avoids gffread's coordinate issues with sliced FASTA by:
1. Reading the sliced FASTA directly
2. Extracting exon sequences for transcripts fully within the slice
3. Handling coordinate translation relative to slice start
4. Reverse complementing minus-strand transcripts
"""

import sys
import re
from collections import defaultdict

def parse_gtf(gtf_path):
    """Parse GTF file and return transcript exons."""
    transcripts = defaultdict(list)  # transcript_id -> [(start, end, strand), ...]
    
    gtf_file = open(gtf_path, 'r')
    for line in gtf_file:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        chrom, source, feature, start_str, end_str, score, strand, frame, attrs = fields[:9]
        
        if feature != 'exon':
            continue
        
        # Parse attributes
        attrs_dict = {}
        for match in re.finditer(r'(\S+)\s+"([^"]+)";', attrs):
            attrs_dict[match.group(1)] = match.group(2)
        
        transcript_id = attrs_dict.get('transcript_id')
        if not transcript_id:
            continue
        
        start = int(start_str)
        end = int(end_str)
        
        transcripts[transcript_id].append({
            'start': start,
            'end': end,
            'strand': strand,
            'chrom': chrom
        })
    
    gtf_file.close()
    return transcripts

def read_fasta(fasta_path):
    """Read FASTA file and return sequence dictionary."""
    sequences = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr is not None:
                    sequences[current_chr] = ''.join(current_seq)
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_chr is not None:
            sequences[current_chr] = ''.join(current_seq)
    
    return sequences

def reverse_complement(seq):
    """Return reverse complement of sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(comp.get(base, base) for base in reversed(seq))

def build_transcriptome(gtf_path, fasta_path, output_path, slice_start=23800000):
    """Build transcriptome FASTA from GTF and sliced FASTA."""
    print(f"Parsing GTF: {gtf_path}")
    transcripts = parse_gtf(gtf_path)
    print(f"  Found {len(transcripts)} transcripts")
    
    print(f"Reading FASTA: {fasta_path}")
    sequences = read_fasta(fasta_path)
    print(f"  Found {len(sequences)} chromosomes")
    
    # Get chromosome name from FASTA (should be chr22 or similar)
    chr_name = list(sequences.keys())[0] if sequences else None
    if not chr_name:
        print("Error: No sequences found in FASTA")
        return False
    
    print(f"  Using chromosome: {chr_name}")
    
    # Get chromosome sequence
    chr_seq = sequences[chr_name]
    chr_len = len(chr_seq)
    print(f"  Chromosome length: {chr_len} bp")
    
    # Build transcript sequences
    output_file = open(output_path, 'w')
    n_success = 0
    n_failed = 0
    
    for transcript_id, exons in sorted(transcripts.items()):
        # Sort exons by start position
        exons_sorted = sorted(exons, key=lambda x: x['start'])
        
        # Get strand (should be consistent across exons)
        strand = exons_sorted[0]['strand']
        
        # Check if all exons are within the slice
        # Coordinates in GTF are 1-based, slice coordinates are relative to slice_start
        all_within = True
        transcript_seq_parts = []
        
        for exon in exons_sorted:
            # Translate coordinates: GTF coordinates are absolute, need to convert to slice-relative
            exon_start_abs = exon['start']
            exon_end_abs = exon['end']
            
            # Convert to slice-relative coordinates (0-based for Python)
            exon_start_rel = exon_start_abs - slice_start - 1  # -1 for 0-based
            exon_end_rel = exon_end_abs - slice_start  # end is inclusive in GTF
            
            # Check bounds
            if exon_start_rel < 0 or exon_end_rel > chr_len:
                all_within = False
                break
            
            # Extract sequence
            exon_seq = chr_seq[exon_start_rel:exon_end_rel]
            transcript_seq_parts.append(exon_seq)
        
        if not all_within:
            n_failed += 1
            continue
        
        # Concatenate exons
        transcript_seq = ''.join(transcript_seq_parts)
        
        # Reverse complement if minus strand
        if strand == '-':
            transcript_seq = reverse_complement(transcript_seq)
        
        # Write FASTA entry
        output_file.write(f">{transcript_id}\n")
        # Write sequence in 60-char lines
        for i in range(0, len(transcript_seq), 60):
            output_file.write(transcript_seq[i:i+60] + '\n')
        
        n_success += 1
    
    output_file.close()
    
    print(f"\nResults:")
    print(f"  Successfully extracted: {n_success} transcripts")
    print(f"  Failed (out of bounds): {n_failed} transcripts")
    print(f"  Output written to: {output_path}")
    
    return n_success > 0

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: build_transcriptome.py <gtf_file> <fasta_file> <output_fasta> [slice_start]")
        print("  slice_start: Start coordinate of FASTA slice (default: 23800000)")
        sys.exit(1)
    
    gtf_path = sys.argv[1]
    fasta_path = sys.argv[2]
    output_path = sys.argv[3]
    slice_start = int(sys.argv[4]) if len(sys.argv) > 4 else 23800000
    
    success = build_transcriptome(gtf_path, fasta_path, output_path, slice_start)
    sys.exit(0 if success else 1)
