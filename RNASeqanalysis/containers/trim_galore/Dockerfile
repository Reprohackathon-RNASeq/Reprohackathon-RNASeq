FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    wget \
    python3.6 \
    python3.6-dev \
    python3-pip \
    unzip \
    build-essential \
    gcc

RUN ln -sf python3.6 /usr/bin/python3 && \
    ln -sf python3.6 /usr/bin/python

RUN python3.6 -m pip install --upgrade pip
RUN python3.6 -m pip install cutadapt==1.11

RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.4.zip && \
    unzip 0.4.4.zip && \
    mv TrimGalore-0.4.4/trim_galore /usr/local/bin/ && \
    chmod +x /usr/local/bin/trim_galore

ENTRYPOINT ["trim_galore"]

# Pour build l'image :
    # docker build -t trim_galore:latest .