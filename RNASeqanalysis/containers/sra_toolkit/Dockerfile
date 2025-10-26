#Testé, approuvé et mis sur docker hub!

FROM ubuntu:latest

# Installer les dépendances nécessaires
RUN apt-get update && \
    apt-get install --yes --quiet curl tar gzip uuid-runtime perl wget && \
    apt-get clean 

# Télécharger et installer la dernière version Ubuntu du SRA Toolkit dans un dossier fixe
RUN mkdir -p /usr/local/sratoolkit && \
    curl -L https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz | \
    tar xz --strip-components=1 -C /usr/local/sratoolkit

# Ajouter le chemin des binaires au PATH
ENV PATH=/usr/local/sratoolkit/bin:${PATH}

#Vérification
RUN fasterq-dump --version
