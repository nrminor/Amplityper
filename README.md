# Amplityper: PCR Amplicon-driven Haplotype Phasing
[![Docker CI](https://github.com/nrminor/Amplityper/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/nrminor/Amplityper/actions/workflows/docker-image.yaml) [![Open Source Starter Files](https://github.com/nrminor/Amplityper/actions/workflows/open-source-starter.yaml/badge.svg)](https://github.com/nrminor/Amplityper/actions/workflows/open-source-starter.yaml) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

## Overview
Many library preparations use PCR amplification to provide ample template for the sequencer. But this amplification also has bioinformatic benefits: By filtering out reads that don't have a given amplicon's forward and reverse primer, we can amass a dataset comprised exclusively of complete and therefore less error-prone amplicons. Amplityper uses this approach to supply clean, accurate reads to a stringent read deduplication process. The resulting haplotypes, derived from a "stack" of identical and complete amplicon sequences, are each likely to represent a real viral lineage that is present in a given infection sample. Amplityper's approach helps overcome the fact that viral sequencing, far from sequencing a single "individual," is really a population sample, with many co-occurring lineages competing and evolving at once.

**Please note:** Amplityper currently only supports Illumina paired-end short reads, though support for PacBio and Oxford Nanopore long reads will be added in the future.

## Quick Start

Amplityper has just two core software requirements: Nextflow to manage the data flow and Docker to provide software (though Apptainer also works). If Nextflow and Docker are already installed on your system, no `git clone` commands or other setup is required. Simply run the workflow and reproduce our findings with:

```
nextflow run nrminor/Amplityper -latest \
--fastq_dir /path/to/FASTQ/files.fastq.gz \
--primer_bed /path/to/primers.bed \
--desired_amplicon name_of_amplicon
```

With this command, Nextflow will automatically pull the workflow bundle of files from this GitHub repository and run it.

## Developer setup
Amplityper also allows users to run locally on the user's system. The upside of removing Docker and running in locally installed software is that it can make the workflow faster, as there is no longer any Docker overhead for pulling images and running containers. The downside is that configuring the development environment is more complex. Users will need to use the configuration files in this repo to install the Python environment (using `pyproject.toml` or `requirements.txt`), the R environment (using `renv.lock`), and the Rust environment (using `Cargo.toml`). More detailed instructions will be provided in the future, alongside a `justfile` that can set up the whole development environment for you.
