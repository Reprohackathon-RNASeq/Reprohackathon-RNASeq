FROM ubuntu:22.04

RUN apt-get update && apt-get install -y wget

# Télécharger et installer Subread 1.4.6-p3
RUN wget https://downloads.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz -O subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    tar -xzf subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    cp -r subread-1.4.6-p3-Linux-x86_64/bin/* /usr/local/bin/ && \
    chmod +x /usr/local/bin/* && \
    rm -rf subread-1.4.6-p3-Linux-x86_64* subread-1.4.6-p3-Linux-x86_64.tar.gz
#  paramètre utilisés -t gene -g ID -s 1

ENTRYPOINT ["featureCounts"]
