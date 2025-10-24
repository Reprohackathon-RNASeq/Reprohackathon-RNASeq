#!/bin/bash
set -e

# Lire le SRR depuis la variable d'environnement
if [ -n "$SRA_RUN" ]; then
    echo "Téléchargement du Run unique : $SRA_RUN"
    fasterq-dump "$SRA_RUN" --threads 4 --progress
else
    # Sinon télécharger tous les runs d'un projet
    echo "Téléchargement de tous les runs de l’étude : $SRA_PROJECT"
    curl -o SraRunInfo.csv "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${SRA_PROJECT}&display=runinfo"
    cut -d',' -f1 SraRunInfo.csv | grep SRR > SRR_list.txt

    echo "Liste des RUNs à télécharger :"
    cat SRR_list.txt

    while read run; do
        echo "Téléchargement: $run"
        fasterq-dump "$run" --threads 4 --progress
    done < SRR_list.txt
fi

echo "Téléchargement terminé !"