FROM ubuntu:latest

# Installation de wget et perl
RUN apt-get update && \
    apt-get install -y wget perl && \
    apt-get clean

# Installation d'Entrez Direct (esearch / efetch)
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xvzf edirect.tar.gz && \
    mv edirect /opt/edirect && \
    rm edirect.tar.gz

# Ajout des binaires au PATH
ENV PATH="/opt/edirect:${PATH}"

#Verification
RUN entrez-direct --version
