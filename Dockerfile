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

# Install SKESA (and the boost cpp libraries)
RUN cd opt && \
    wget https://sourceforge.net/projects/boost/files/boost/1.83.0/boost_1_83_0.tar.gz && \
    tar -xvzf boost_1_83_0.tar.gz && \
    rm boost_1_83_0.tar.gz && \
    cd boost_1_83_0 && \
    ./bootstrap.sh && \
    apt-get install -y libbz2-dev && \
    ./b2 && \
    ./b2 install
ENV BOOST_PATH=/opt/boost_1_83_0
RUN apt-get install -y build-essential libtool autoconf cmake && \
    cd /opt && \
    git clone https://github.com/ncbi/SKESA && \
    cd SKESA/ && \
    make
ENV PATH=$PATH:/opt/SKESA:/opt/boost_1_83_0/boost

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

# Install MultiQC
RUN pip install multiqc

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

# Download and install SPAdes
# RUN wget http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz && \
#     tar -xzf SPAdes-3.15.2-Linux.tar.gz && \
#     rm SPAdes-3.15.2-Linux.tar.gz

# Add SPAdes to PATH
# ENV PATH="/SPAdes-3.15.2-Linux/bin:${PATH}"
RUN echo "alias python=python3" >> ~/.bashrc

# Install amplicon_sorter
RUN pip install --no-cache-dir \
    icecream poetry biopython matplotlib && \
    pip wheel --no-cache-dir --use-pep517 "edlib (==1.3.9)" && \
    cd /opt && \
    git clone https://github.com/nrminor/amplicon_sorter.git && \
    cd amplicon_sorter && \
    git checkout dev && \
    chmod +x amplicon_sorter.py
ENV PATH="$PATH:/opt/amplicon_sorter"

# Install read_zap_report
RUN cd /opt && \
    git clone https://github.com/nrminor/read-zap-report.git && \
    cd read-zap-report && \
    pip install -r requirements.txt && \
    chmod +x read_zap_report.py
ENV PATH="$PATH:/opt/read-zap-report/"

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

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
