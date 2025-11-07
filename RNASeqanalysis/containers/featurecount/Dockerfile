FROM ubuntu:22.04

RUN apt-get update && apt-get install -y wget

# Télécharger et installer Subread 1.4.6-p3
RUN wget https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download -O subread-1.4.6-p3-Linux-x86_64.tar.gz && \
#  paramètre utilisés -t gene -g ID -s 1

ENTRYPOINT ["featureCounts"]

