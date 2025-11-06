FROM ubuntu:latest

# Install the necessary dependencies
RUN apt-get update && \
    apt-get install --yes --quiet curl tar gzip uuid-runtime perl wget unzip && \
    apt-get clean 

# Download and install the version of Bowtie required
RUN mkdir /usr/local/bowtie && \
    wget -O /tmp/bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip && \
    unzip /tmp/bowtie.zip -d /usr/local/bowtie && \
    rm /tmp/bowtie.zip

# Add the binary path to the PATH
ENV PATH="/usr/local/bowtie:/usr/local/bowtie/bowtie-0.12.7:${PATH}"

# Vérification du contenu
RUN ls -R /usr/local/bowtie

# Vérification de l'installation
RUN /usr/local/bowtie/bowtie-0.12.7/bowtie --version
