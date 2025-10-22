#!/bin/bash
set -e

# 1. Télécharge le tableau SraRunInfo.csv depuis NCBI automatiquement
curl -o SraRunInfo.csv "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${SRA_PROJECT}&display=runinfo"

# 2. Liste tous les run IDs associés à l'étude
cut -d',' -f1 SraRunInfo.csv | grep SRR > SRR_list.txt

echo "Liste des RUNs à télécharger :"
cat SRR_list.txt

# 3. Pour chaque run, télécharge le FASTQ
while read run; do
    echo "Téléchargement: $run"
    fasterq-dump $run --threads 4 --progress
done < SRR_list.txt

echo "Téléchargement terminé !"