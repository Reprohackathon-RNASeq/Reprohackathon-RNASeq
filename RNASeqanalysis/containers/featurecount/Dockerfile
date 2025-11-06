FROM ubuntu:22.04

# Install the necessary dependencies
RUN apt-get update && apt-get install -y wget

# Download and install Subread 1.4.6-p3 in a fixed folder
RUN mkdir /usr/local/subread && \
    wget -O /tmp/subread.tar.gz "https://versaweb.dl.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz" && \
    tar -zxvf /tmp/subread.tar.gz -C /usr/local/subread --strip-components=1 && \
    rm /tmp/subread.tar.gz

# Add the binary path to the PATH
ENV PATH=/usr/local/subread/bin:${PATH}

# Verification
RUN featureCounts -v
