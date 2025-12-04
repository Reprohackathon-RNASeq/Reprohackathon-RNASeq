# ReproHackathon - RNA-Seq Analysis of Staphylococcus aureus persisters

## Project Description

This project is part of a **ReproHackathon** aiming to reproduce the RNA-Seq analysis from the article:

> **"Intracellular Staphylococcus aureus persisters upon antibiotic exposure"**  
> *Nature Communications* (2020) 11:2200  
> DOI: [10.1038/s41467-020-15966-7](https://doi.org/10.1038/s41467-020-15966-7)

Analysing RNA-Seq data to identify differentially expressed genes in intracellular persister cells of *Staphylococcus aureus* after antibiotic treatment.

This is the result of a course from the AMI2B master's program at Paris-Saclay University.

## Biological Context

Bacterial **persisters** are phenotypic variants that:
- Enter a **transient non-dividing state**
- Develop **antibiotic tolerance**
- Contribute to **chronic infections** and **therapeutic failures**
- Constitute a **reservoir for relapses**

## RNA-Seq Analysis Pipeline Overview

The workflow is implemented using Nextflow and encapsulated in Docker containers to ensure reproducibility. It includes:
* Retrieval of SRA sequencing data
* Quality trimming of reads
* Reference genome preparation and indexing
* Read alignment with Bowtie
* Gene-level quantification with featureCounts
* Differential expression analysis with DESeq2
* Annotation of differentially expressed genes using AureoWiki and KEGG resources

![Pipeline Overview](Pipeline.png)

## Installation and Execution
    
### Prerequisites

The pipeline is designed to run on a virtual machine with **Docker, Nextflow and Java** running.

### Execution

Once on the VM, please run these commands to start the process.

```
cd /home/ubuntu/data/mydatalocal
git clone https://github.com/Reprohackathon-RNASeq/Reprohackathon-RNASeq.git
cd Reprohackathon-RNASeq/RNASeqanalysis
chmod +x run.sh
./run.sh
```

## Contributors

Tom BELLIVIER
Tom GORTANA
Marie MEIER
Donatien WALLAERT
