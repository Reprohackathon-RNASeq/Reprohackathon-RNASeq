FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    wget \
    python3 \
    python3-pip \
    unzip \
    build-essential \
    python3-dev \
    gcc \
    cutadapt

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC/fastqc /usr/local/bin/

RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.zip && \
    unzip 0.6.10.zip && \
    mv TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    chmod +x /usr/local/bin/trim_galore

ENTRYPOINT ["trim_galore"]

# Pour build l'image :
    #docker build -t trim_galore:latest .