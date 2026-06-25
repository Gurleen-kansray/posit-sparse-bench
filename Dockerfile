FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    g++ \
    cmake \
    git \
    python3 \
    python3-pip \
    wget \
    curl \
    libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --branch v3.80 --depth 1 https://github.com/stillwater-sc/universal.git /universal

ENV UNIVERSAL_INCLUDE=/universal/include

WORKDIR /work
COPY posit-sparse-bench/src/ ./src/
COPY posit-sparse-bench/Makefile ./Makefile
COPY posit-sparse-bench/data/ ./data/

RUN mkdir -p results

RUN make INCLUDES=-I/universal/include
