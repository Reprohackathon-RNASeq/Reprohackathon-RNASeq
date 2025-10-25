#Testé et approuvé !

FROM ubuntu:latest

# Installer les dependances necessaires
RUN apt-get update && \
    apt-get install --yes --quiet curl tar gzip uuid-runtime perl wget && \
    apt-get clean 

# Telecharger et installer la derniere version Ubuntu du SRA Toolkit
RUN curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz | tar xz -C /usr/local

# Ajouter le chemin des binaires au PATH
ENV PATH=/usr/local/sratoolkit.3.2.1-ubuntu64/bin:${PATH}

# Installer Entrez Direct (esearch / efetch)
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xvzf edirect.tar.gz && \
    mv edirect /opt/edirect && \
    rm edirect.tar.gz

# Ajouter le chemin des binaires au PATH
ENV PATH="/opt/edirect:${PATH}"

#Commandes bash
#docker build -t sra-toolkit .
#docker run -it sra-toolkit bash

# Vérification
#prefetch --version
#fasterq-dump --version
#esearch -version

#Pour un nettoyage complet de Docker 
#docker system prune -a --volumes -f

#Pour envoyer l'image vers Docker Hub
#docker login
#docker tag sra-toolkit your_dockerhub_username/sra-toolkit:latest
#docker push your_dockerhub_username/sra-toolkit:latest

#Vérifier l'image sur Docker Hub
#docker pull your_dockerhub_username/sra-toolkit:latest
#docker run -it your_dockerhub_username/sra-toolkit:latest bash
