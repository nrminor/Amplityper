# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Avoid prompts from APT during build
ARG DEBIAN_FRONTEND=noninteractive

# Update the package list and install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    make \
    wget \
    curl \
    gzip \
    pigz \
    zstd \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    autotools-dev \
    autoconf \
    default-jre \
    perl \
    python-is-python3 \
    python3-pip \
    cython \
    pypy-dev \
    bash

# Download and extract BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz

# Add BBTools to PATH
ENV PATH="/bbmap:${PATH}"

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

# install Rust
RUN mkdir -m777 /opt/rust /opt/.cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/.cargo PATH=/opt/.cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    bash "/opt/.cargo/env"

# install tidyVCF tool
RUN cargo install tidyvcf

# Install amplicon_sorter
RUN pip install icecream poetry biopython matplotlib edlib && \
    cd /opt && \
    git clone https://github.com/nrminor/amplicon_sorter.git && \
    cd amplicon_sorter && \
    git checkout dev && \
    poetry install && \
    chmod +x amplicon_sorter.py
ENV PATH="$PATH:/opt/amplicon_sorter"

# Install read_zap_report
RUN cd /opt && \
    git clone https://github.com/nrminor/read-zap-report.git && \
    cd read-zap-report && \
    pip install -r requirements.txt && \
    chmod +x read_zap_report.py

# Add SPAdes to PATH
# ENV PATH="/SPAdes-3.15.2-Linux/bin:${PATH}"
RUN echo "alias python=python3" >> ~/.bashrc

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
