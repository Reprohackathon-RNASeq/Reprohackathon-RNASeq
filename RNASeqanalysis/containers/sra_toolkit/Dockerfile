FROM ubuntu:latest

# Install the necessary dependencies
RUN apt-get update && \
    apt-get install --yes --quiet curl tar gzip uuid-runtime perl wget && \
    apt-get clean 

# Download and install the latest Ubuntu version of the SRA Toolkit in a fixed folder
RUN mkdir -p /usr/local/sratoolkit && \
    curl -L https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz | \
    tar xz --strip-components=1 -C /usr/local/sratoolkit

# Add the binary path to the PATH
ENV PATH=/usr/local/sratoolkit/bin:${PATH}

# Verification
RUN fasterq-dump --version
 