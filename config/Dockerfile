# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Avoid prompts from APT during build
ARG DEBIAN_FRONTEND=noninteractive

# Update the package list and install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    wget \
    gzip \
    zstd \
    default-jre \
    perl \
    python3 \
    python3-pip \
    pypy-dev \
    bash \
    samtools \
    bcftools

# Download and extract BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz

# Add BBTools to PATH
ENV PATH="/bbmap:${PATH}"

# Install PyPy
RUN apt-get install -y pypy3
# RUN pypy3 -m ensurepip
# RUN pypy3 -m pip install biopython argparse

# Install the required Python modules
RUN pip3 install biopython argparse pysam

# Download and install SPAdes
RUN wget http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz && \
    tar -xzf SPAdes-3.15.2-Linux.tar.gz && \
    rm SPAdes-3.15.2-Linux.tar.gz

# Add SPAdes to PATH
ENV PATH="/SPAdes-3.15.2-Linux/bin:${PATH}"

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
