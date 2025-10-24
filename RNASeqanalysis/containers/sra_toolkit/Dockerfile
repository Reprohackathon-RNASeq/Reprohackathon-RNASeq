FROM ubuntu:22.04

# Installe les outils nécessaires
RUN apt-get update && apt-get install -y \
    wget curl tar gzip python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Télécharge et installe le SRA Toolkit (version 3.1.1)
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.3.1.1-ubuntu64.tar.gz && \
    mv sratoolkit.3.1.1-ubuntu64 /opt/sratoolkit && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm sratoolkit.3.1.1-ubuntu64.tar.gz

# Crée le dossier de travail
WORKDIR /data

# Copie le script de téléchargement
COPY download_sra.sh /usr/local/bin/download_sra.sh
RUN chmod +x /usr/local/bin/download_sra.sh

ENTRYPOINT ["/usr/local/bin/download_sra.sh"]
