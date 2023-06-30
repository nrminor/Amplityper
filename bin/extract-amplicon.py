#!/usr/bin/env python3

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
    args = parser.parse_args()
    return args.bam_file, args.bed_file, args.amplicon_name


def resolve_symlink(path: str) -> str:
    """Resolve symbolic links in the provided path."""
    return os.path.realpath(path)


def parse_bed_file(bed_file: str, amplicon_name: str) -> Tuple[str, int, int]:
    """Parse the BED file and find amplicon coordinates by name."""
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end, name = line.strip().split('\t')
            if name.contains(amplicon_name):
                return chrom, int(start), int(end)
    raise ValueError(f"Amplicon '{amplicon_name}' not found in the BED file.")


def extract_mapped_reads(bam_file: str, chrom: str, start: int, end: int) -> List[pysam.AlignedSegment]:
    """Extract mapped reads within the specified amplicon coordinates."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    mapped_reads = [read for read in bam.fetch(chrom, start, end) if not read.is_unmapped]
    bam.close()
    return mapped_reads


def write_mapped_reads_to_bam(mapped_reads: List[pysam.AlignedSegment], output_bam: str):
    """Write the extracted mapped reads to a new BAM file."""
    header = pysam.AlignmentHeader.from_template(mapped_reads[0].header)
    with pysam.AlignmentFile(output_bam, "wb", header=header) as output:
        for read in mapped_reads:
            output.write(read)


def main():
    # Parse command line arguments
    bam_file, bed_file, amplicon_name = parse_command_line_args()

    # Resolve symbolic links in file paths
    bam_file = resolve_symlink(bam_file)
    bed_file = resolve_symlink(bed_file)

    # Parse BED file to find amplicon coordinates
    chrom, start, end = parse_bed_file(bed_file, amplicon_name)

    # Extract mapped reads within amplicon coordinates
    mapped_reads = extract_mapped_reads(bam_file, chrom, start, end)

    # Write mapped reads to a new BAM file
    output_bam = f"{amplicon_name}_mapped_reads.bam"
    write_mapped_reads_to_bam(mapped_reads, output_bam)
    print(f"Extracted mapped reads written to {output_bam}")


if __name__ == "__main__":
    main()
