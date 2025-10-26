FROM ubuntu:latest

# Installer wget et perl
RUN apt-get update && \
    apt-get install -y wget perl && \
    apt-get clean

# Installer Entrez Direct dans /usr/local/bin
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xvzf edirect.tar.gz && \
    mv edirect/* /usr/local/bin/ && \
    rm -rf edirect edirect.tar.gz

# Vérification
RUN esearch -version

