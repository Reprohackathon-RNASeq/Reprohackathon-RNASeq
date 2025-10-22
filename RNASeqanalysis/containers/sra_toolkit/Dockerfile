FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y pigz wget tar && \
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu22.04-x64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu22.04-x64.tar.gz && \
    mv sratoolkit.* /opt/sratoolkit && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm sratoolkit.current-ubuntu22.04-x64.tar.gz

WORKDIR /data

# Définis l’ID d’étude en variable ENV 
ENV SRA_PROJECT=SRP227811

# Ajoute un script qui télécharge tous les Runs de l’étude au lancement
COPY download_sra.sh /usr/local/bin/download_sra.sh

ENTRYPOINT ["/usr/local/bin/download_sra.sh"]

## Pour run le build:
# docker build --platform linux/amd64 -t sratoolkit_dl .

