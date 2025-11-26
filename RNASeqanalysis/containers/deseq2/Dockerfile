# Base image
FROM ubuntu:20.04 

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update 
RUN apt-get install -y \
    curl \
    make \
    zlib1g-dev \
    zlib1g \
    gfortran \
    g++ \
    gcc \
    aptitude \
    fort77 \
    libbz2-dev \
    liblzma-dev \
    libpcre3 \
    libpcre3-dev \
    libcurl4 \
    libcurl4-openssl-dev \
    default-jdk \
    texlive-latex-base \
    libxml2-dev \
    libpng-dev \
    graphviz \
    libgraphviz-dev \
    libssl-dev \
    r-base-dev \
    libcairo-dev

# Additional system dependency
RUN aptitude install -y libreadline-dev

# Download specific R version: 3.4.1
RUN curl -o R_3-4-1.tar.gz https://cran.rstudio.com/src/base/R-3/R-3.4.1.tar.gz \
    && tar -xzf R_3-4-1.tar.gz \
    && rm R_3-4-1.tar.gz 

# Compile and install R
RUN cd R-3.4.1/ && ./configure && make && make install

# Install CRAN packages
RUN R -e "install.packages(c('tidyverse','ggrepel','pheatmap','Cairo'), repos='https://packagemanager.posit.co/cran/__linux__/focal/2018-01-11')"

# Install BiocInstaller for installing Bioconductor packages
RUN R -e "install.packages('BiocInstaller', repos='https://bioconductor.org/packages/3.5/bioc/')"

# Install dependencies for Rgraphviz
RUN R -e "BiocInstaller::biocLite('graph')"

# Download Rgraphviz source (version 2.54)
RUN curl -O https://www.bioconductor.org/packages/release/bioc/src/contrib/Rgraphviz_2.54.0.tar.gz \
    && tar -xzf Rgraphviz_2.54.0.tar.gz \
    && rm Rgraphviz_2.54.0.tar.gz

# Install Rgraphviz from source
RUN R -e "install.packages('Rgraphviz/', repos = NULL, type='source')" \
    && rm -rf Rgraphviz

# Install DESeq2 (version 1.16.1)
RUN R -e "options(repos = c(CRAN='https://packagemanager.posit.co/cran/__linux__/focal/2018-01-11')); BiocInstaller::biocLite('DESeq2')"

# Install EnrichmentBrowser
RUN R -e "BiocInstaller::biocLite('EnrichmentBrowser')"

# Test that DESeq2 loads correctly
RUN Rscript -e 'library("DESeq2")'
