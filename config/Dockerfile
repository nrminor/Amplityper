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

# Add the above to PATH
ENV PATH="$PATH:/htslib-1.17:/bcftools-1.17:/samtools-1.17"

# Install PyPy
RUN apt-get install -y pypy3
# RUN pypy3 -m ensurepip
# RUN pypy3 -m pip install biopython argparse

# Install the required Python modules
# RUN pip3 install biopython argparse cython pysam
RUN export HTSLIB_CONFIGURE_OPTIONS=--enable-plugins && \
    pip3 install biopython argparse cython pysam

# Download and install SPAdes
RUN wget http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz && \
    tar -xzf SPAdes-3.15.2-Linux.tar.gz && \
    rm SPAdes-3.15.2-Linux.tar.gz

# Add SPAdes to PATH
ENV PATH="/SPAdes-3.15.2-Linux/bin:${PATH}"
RUN echo "alias python=python3" >> ~/.bashrc

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
