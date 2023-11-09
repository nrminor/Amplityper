# READ-ZAP: Read Extraction And De novo assembly for Zero-error Amplicon Phasing
[![Docker CI](https://github.com/nrminor/READ-ZAP/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/nrminor/READ-ZAP/actions/workflows/docker-image.yaml)

## Overview
Many library preparations use PCR amplification to provide ample template for the sequencer. But this amplification also has bioinformatic benefits: By filtering out reads that don't have a given amplicon's forward and reverse primer, we can amass a dataset comprised exclusively of complete and therefore less error-prone amplicons. READ-ZAP uses this approach to supply clean, accurate reads to a stringent de novo assembly algorithm. The resulting contigs, derived from a "stack" of identical amplicon sequences, are each likely to represent a real viral lineage that is present in a given infection. READ-ZAP's approach helps overcome the fact that viral sequencing, far from sequencing a single "individual," is really a population sample, with many co-occurring lineages competing and evolving at once.

**Please note:** READ-ZAP currently only supports Illumina paired-end short reads, though support for PacBio and Oxford Nanopore long reads will be added in the future.

## Quick Start

At this stage, READ-ZAP has three core software requirements: Nextflow to manage the data flow, Docker to provide software (though Apptainer also works), and Geneious Prime for de novo assembly. If Nextflow, Docker, and Geneious Prime are already installed on your system, no `git clone` commands or other setup is required. Simply run the workflow and reproduce our findings with:

```
nextflow run nrminor/READ-ZAP -latest \
--fastq_dir /path/to/FASTQ/files.fastq.gz \
--primer_bed /path/to/primers.bed \
--desired_amplicon name_of_amplicon
```

With this command, Nextflow will automatically pull the workflow bundle of files from this GitHub repository and run it.
