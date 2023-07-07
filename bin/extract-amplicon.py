#!/usr/bin/env python3

import os
import sys
import argparse
import pysam
from typing import List, Tuple

def parse_command_line_args() -> Tuple[str, str, str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Extract mapped reads within amplicon coordinates.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("bed_file", help="Path to the BED file containing amplicon coordinates.")
    parser.add_argument("amplicon_name", help="Name of the amplicon of interest.")
    parser.add_argument("fwd_primer_suffix", help="Extension on the names of the forward primers (e.g., _LEFT).")
    parser.add_argument("rev_primer_suffix", help="Extension on the names of the forward primers (e.g., _RIGHT).")
    args = parser.parse_args()
    return args.bam_file, args.bed_file, args.amplicon_name, args.fwd_primer_suffix, args.rev_primer_suffix


def resolve_symlink(path: str) -> str:
    """Resolve symbolic links in the provided path."""
    return os.path.realpath(path)


def parse_bed_file(bed_file: str, amplicon_name: str, fwd_primer: str, rev_primer: str) -> Tuple[str, int, int]:
    """Parse the BED file and find amplicon coordinates by name."""
    with open(bed_file, 'r') as bed:
        ref = None
        amplicon_start = None
        amplicon_end = None
        
        for line in bed:
            chrom, start, end, name, num, strand = line.strip().split('\t')
            if amplicon_name in name and fwd_primer in name:
                ref = chrom
                amplicon_start = int(end)
                continue
            elif amplicon_name in name and rev_primer in name:
                amplicon_end = int(start)
            if ref is not None and amplicon_start is not None and amplicon_end is not None:
                break
            
        if ref is not None and amplicon_start is not None and amplicon_end is not None:
            return ref, amplicon_start, amplicon_end
        else:
            raise ValueError(f"Amplicon '{amplicon_name}' not found in the BED file.")


def extract_mapped_reads(bam_file: str, chrom: str, start: int, end: int) -> List[pysam.AlignedSegment]:
    """Extract mapped reads within the specified amplicon coordinates."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    assert chrom in bam.references, f"Contig '{chrom}' not found in the BAM file."
    
    mapped_reads = []
    for read in bam.fetch(chrom, start, end):
        if not read.is_unmapped and read.reference_start >= start and read.reference_end <= end:
            # Filter out reads that align outside the amplicon coordinates or that have 
            # "junction" annotations
            if "junction" not in read.query_name:
                mapped_reads.append(read)
    
    bam.close()
    return mapped_reads


def write_mapped_reads_to_bam(mapped_reads: List[pysam.AlignedSegment], output_bam: str):
    """Write the extracted mapped reads to a new BAM file."""
    header = mapped_reads[0].header.to_dict()  # Get the header as a dictionary
    new_header = pysam.AlignmentHeader.from_dict(header)  # Create a new AlignmentHeader object
    with pysam.AlignmentFile(output_bam, "wb", header=new_header) as output:
        for read in mapped_reads:
            output.write(read)


def main():
    # Parse command line arguments
    bam_file, bed_file, amplicon_name, fwd_primer, rev_primer = parse_command_line_args()

    # Resolve symbolic links in file paths
    bam_file = resolve_symlink(bam_file)
    bed_file = resolve_symlink(bed_file)

    # Parse BED file to find amplicon coordinates
    chrom, start, end = parse_bed_file(bed_file, amplicon_name, fwd_primer, rev_primer)

    # Extract mapped reads within amplicon coordinates
    mapped_reads = extract_mapped_reads(bam_file, chrom, start, end)

    # Write mapped reads to a new BAM file
    output_bam = f"{amplicon_name}_mapped_reads.bam"
    write_mapped_reads_to_bam(mapped_reads, output_bam)
    print(f"Extracted mapped reads written to {output_bam}")


if __name__ == "__main__":
    main()
