# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Avoid prompts from APT during build
ARG DEBIAN_FRONTEND=noninteractive

# Update the package list and install necessary packages
RUN apt update && \
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
    bash
RUN apt-get install -y pkg-config

# install Rust
RUN mkdir -m777 /opt/rust /opt/.cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/.cargo PATH=/opt/.cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    bash "/opt/.cargo/env"

# install tidyVCF tool
RUN cargo install tidyvcf

# Install fastqc-rs
RUN apt-get install -y libssl-dev && \
    cargo install fastqc-rs

# Install HTSLib
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar -vxjf htslib-1.17.tar.bz2 && \
    rm htslib-1.17.tar.bz2 && \
    cd htslib-1.17 && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Install SAMTools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 &&\
    tar -vxjf samtools-1.17.tar.bz2 && \
    rm samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Install BCFTools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 &&\
    tar -vxjf bcftools-1.17.tar.bz2 && \
    rm bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz && \
    tar -xzf bedtools-2.31.0.tar.gz && \
    cd bedtools2 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bedtools2 bedtools-2.31.0.tar.gz

# Add the above to PATH
ENV PATH="$PATH:/htslib-1.17:/bcftools-1.17:/samtools-1.17"

# Install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz \
    && tar -xvzf seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin/ \
    && rm seqkit_linux_amd64.tar.gz

# Install csvtk
RUN wget https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz \
    && tar -xvzf csvtk_linux_amd64.tar.gz \
    && mv csvtk /usr/local/bin/ \
    && rm csvtk_linux_amd64.tar.gz

# Install vsearch
RUN wget https://github.com/torognes/vsearch/archive/v2.26.1.tar.gz && \
    tar xzf v2.26.1.tar.gz && \
    cd vsearch-2.26.1 && \
    ./autogen.sh && \
    ./configure CFLAGS="-O3" CXXFLAGS="-O3" && \
    make && \
    make install

# Download and extract BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz

# Add BBTools to PATH
ENV PATH="/bbmap:${PATH}"

# Install snpEff
RUN cd /opt && \
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    rm snpEff_latest_core.zip && \
    chmod +x snpEff/exec/*
ENV PATH=$PATH:/opt/snpEff/exec

# make sure the image uses the correct version of Python
RUN cd /usr/bin && \
    unlink python && \
    ln -s /usr/bin/python3.12 python && \
    unlink python3 && \
    ln -s /usr/bin/python3.12 python3

# Install MultiQC and a few other things
RUN apt install -y python3.12-distutils
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3
# RUN python3 -m pip install multiqc
COPY requirements.txt /opt/
RUN python3 -m pip install -r /opt/requirements.txt

# install R
RUN apt install -y r-base && \
    Rscript -e "install.packages('tidyverse', clean = TRUE)"

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
