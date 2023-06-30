#!/usr/bin/env python3

import sys
import gzip
from typing import List

### POSSIBLE OPTIMIZATIONS ###
### ---------------------- ###
#   - read FASTQ in batches rather than all at once
#   - remove condition() function
#   - multiprocess FASTQ if performance is an issue.
### ---------------------- ###

def parse_bed_file(bed_file_path: str, amplicon_name: str) -> int:
    amplicon_length = None
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end, name = line.strip().split('\t')
            if name.contains(amplicon_name):
                amplicon_length = int(end) - int(start)
                break
    return amplicon_length

def filter_fastq(fastq_file_path: str, amplicon_length: int) -> List[str]:
    def condition(i: int, line: str) -> bool:
        return (i + 2) % 4 == 0 and len(line.strip()) == amplicon_length

    with gzip.open(fastq_file_path, 'rt') as fastq_file:
        lines = fastq_file.readlines()
        filtered_reads = [line for i, line in enumerate(lines) if condition(i, line)]
    return filtered_reads

def resolve_symlink(path: str) -> str:
    if os.path.islink(path):
        return os.path.realpath(path)
    return path

def main():
    # Get command line arguments
    fastq_file_path = resolve_symlink(sys.argv[1])
    bed_file_path = resolve_symlink(sys.argv[2])
    amplicon_name = sys.argv[3]

    # Parse BED file to determine amplicon length
    amplicon_length = parse_bed_file(bed_file_path, amplicon_name)
    if amplicon_length is None:
        print(f"Amplicon '{amplicon_name}' not found in the BED file.")
        return

    # Filter FASTQ reads based on amplicon length
    filtered_reads = filter_fastq(fastq_file_path, amplicon_length)

    # Output filtered reads to stdout
    for read in filtered_reads:
        print(read, end='')

if __name__ == '__main__':
    main()
