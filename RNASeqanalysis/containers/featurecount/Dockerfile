FROM ubuntu:22.04

# Install the necessary dependencies
RUN apt-get update && \
    apt-get install -y wget tar && \
    rm -rf /var/lib/apt/lists/*

# Download and install Subread in a fixed folder
RUN mkdir -p /usr/local/subread && \
    wget https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download -O subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    tar -xzf subread-1.4.6-p3-Linux-x86_64.tar.gz -C /usr/local/subread --strip-components 1 && \
    rm subread-1.4.6-p3-Linux-x86_64.tar.gz

# Add the binary path to the PATH
ENV PATH="/usr/local/subread/bin:${PATH}"

# Verification
RUN featureCounts -v
