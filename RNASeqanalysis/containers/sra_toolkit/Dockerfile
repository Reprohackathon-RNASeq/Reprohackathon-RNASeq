FROM ubuntu:latest

# Installer les dependances necessaires
RUN apt-get update && \
    apt-get install --yes --quiet curl tar gzip uuid-runtime && \
    apt-get clean 

# Telecharger et installer la derniere version Ubuntu du SRA Toolkit
RUN curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz | tar xz -C /usr/local

# Ajouter le chemin des binaires au PATH
ENV PATH=/usr/local/sratoolkit.3.2.1-ubuntu64/bin:${PATH}

#Commandes bash
#docker build -t sra-toolkit .
#docker run -it sra-toolkit bash

# Vérification
#prefetch --version
#fasterq-dump --version
