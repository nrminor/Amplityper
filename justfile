default:
    just --list

# Install the Rust toolchain and the crates used by Amplityper
rust:
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    cargo install tidyvcf
    cargo install fastqc-rs
    cargo install nanoq
    cargo install scidataflow
    cargo install nushell
alias rs := rust

# Install Intel-chip MacOS packages available via homebrew
homebrew:
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" && \
    -brew install \
    wget \
    curl \
    zstd \
    pigz \
    unzip \
    java \
    samtools \
    bcftools \
    bedtools \
    seqkit \
    csvtk \
    vsearch \
    fastp \
    r
    -brew install --cask docker
alias brew := homebrew
alias hb := homebrew

# Install Ubuntu packages via apt-get
ubuntu-apt-get:
    apt update && \
    apt install software-properties-common -y && \
    apt update && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt update && \
    apt-get update && \
    apt-get install -y \
    build-essential \
    gcc \
    make \
    wget \
    curl \
    gzip \
    pigz \
    zstd \
    unzip \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    autotools-dev \
    autoconf \
    default-jre \
    perl \
    python3.12 \
    python3-pip \
    cython \
    git \
    bash \
    pkg-config \
    libssl-dev \
    python3.12-distutils \
    r-base

# Packages to build from source localled
local-builds:
    touch ~/.zprofile
    -mkdir ~/bioinformatics
    cd ~/bioinformatics
    wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz
    -tar -xzf bbmap.tar.gz
    rm bbmap.tar.gz
    echo "export PATH=$PATH:~/bioinformatics/bbmap" >> ~/.zprofile
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    -unzip snpEff_latest_core.zip
    rm snpEff_latest_core.zip
    chmod +x snpEff/exec/*
    echo "export PATH=$PATH:~/bioinformatics/snpEff/exec" >> ~/.zprofile
alias lb := local-builds

# R libraries
r-packages:
    Rscript -e "install.packages('tidyverse',  repos='http://cran.us.r-project.org', clean = TRUE)"
    Rscript -e "install.packages('arrow',  repos='http://cran.us.r-project.org', clean = TRUE)"
alias r := r-packages

# Install Python packages for x86 (Intel chips)
python-x86:
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh
    mamba install -y -c conda-forge python==3.12.0 poetry==1.7.1
    python3 -m pip install -r requirements.txt
    poetry install
    @echo "It is recommended to run `poetry shell` in this directory before using Amplityper."
alias py86 := python-x86

# Install Python packages for ARM (Apple chips)
python-arm:
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh
    mamba install -y -c conda-forge python==3.12.0 poetry==1.7.1
    python3 -m pip install -r requirements.txt
    poetry install
    @echo "It is recommended to run `poetry shell` in this directory before using Amplityper."
alias pyarm := python-arm

[macos, confirm]
# install all packages for Intel chip Macs
macos-x86:
    just homebrew-x86
    just local-builds
    just rust
    just r-packages
    just python-x86

[macos, confirm]
# install all packages for Apple chip Macs
macos-arm:
    just homebrew-arm
    just local-builds
    just rust
    just r-packages
    just python-arm

[linux, confirm]
# Install all packages for Linux ubuntu systems
ubuntu:
    just ubuntu-apt-get
    just rust
    just r-packages
    just python-x86
    just local-builds
