FROM ubuntu:20.04

LABEL version=""
LABEL description="bam2bakR environment installation"
LABEL maintainer="isaac.vock@gmail.com"

# Setup ubuntu basic softwares
RUN apt-get update \
 && apt-get install -y wget git nano gcc g++ libz-dev bedtools \
 && rm -rf /var/lib/apt/lists/*

# Install miniconda