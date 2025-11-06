# Older version of TrimGalore for compatibility reasons
FROM ubuntu:18.04

# Install the necessary dependencies
RUN apt-get update && apt-get install -y \
    wget \
    python3.6 \
    python3.6-dev \
    python3-pip \
    unzip \
    build-essential \
    gcc

# Configure the Python Interpreter
RUN ln -sf python3.6 /usr/bin/python3 && \
    ln -sf python3.6 /usr/bin/python

# Install Cutadapt
RUN python3.6 -m pip install --upgrade pip
RUN python3.6 -m pip install cutadapt==1.11

# Download and install the version of TrimGalore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.4.zip && \
    unzip 0.4.4.zip && \
    mv TrimGalore-0.4.4/trim_galore /usr/local/bin/ && \
    chmod +x /usr/local/bin/trim_galore

# Verification
RUN trim_galore --version