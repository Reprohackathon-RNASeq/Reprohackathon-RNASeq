#Testé, approuvé et mis sur docker hub!

FROM ubuntu:latest

# Installer wget et perl
RUN apt-get update && \
    apt-get install -y wget perl && \
    apt-get clean

# Installer Entrez Direct dans /usr/local/bin
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xvzf edirect.tar.gz && \
    mkdir -p /usr/local/entrez-direct/bin && \
    mv edirect/* /usr/local/entrez-direct/bin/ && \
    rm -rf edirect edirect.tar.gz

ENV PATH="/usr/local/entrez-direct/bin:${PATH}"

# Vérification
RUN esearch -version
