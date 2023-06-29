# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Avoid prompts from APT during build
ARG DEBIAN_FRONTEND=noninteractive

# Update the package list and install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    default-jre \
    perl \
    python3 \
    bash

# Create a directory for BBTools
WORKDIR /bbmap

# Download and extract BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz && \
    ./build.sh

# Add BBTools to PATH
ENV PATH="/bbmap:${PATH}"

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
