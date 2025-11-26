FROM ubuntu:latest

# Install the necessary dependencies
RUN apt-get update && \
    apt-get install -y openjdk-11-jre-headless wget unzip perl libio-compress-perl && \
    apt-get clean

# Download and install FastQC
RUN mkdir -p /usr/local/fastqc && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    mv FastQC/ /usr/local/fastqc/FastQC_0.11.9 && \
    chmod 755 /usr/local/fastqc/FastQC_0.11.9/fastqc && \
    rm -rf fastqc_v0.11.9.zip

# Add the binary path to the PATH
ENV PATH=/usr/local/fastqc/FastQC_0.11.9:${PATH}

# Verification
RUN fastqc --version