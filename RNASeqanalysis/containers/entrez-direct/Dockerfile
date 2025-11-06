FROM ubuntu:latest

# Install the necessary dependencies
RUN apt-get update && \
    apt-get install -y wget perl && \
    apt-get clean

# Install Entrez Direct in /usr/local/bin
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xvzf edirect.tar.gz && \
    mkdir -p /usr/local/entrez-direct/bin && \
    mv edirect/* /usr/local/entrez-direct/bin/ && \
    rm -rf edirect edirect.tar.gz

ENV PATH="/usr/local/entrez-direct/bin:${PATH}"

# Verification
RUN esearch -version
