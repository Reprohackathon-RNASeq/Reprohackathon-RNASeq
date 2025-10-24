#!/bin/bash
set -e

if [ -z "$1" ]; then
  echo "❌ Erreur : merci de fournir un identifiant SRA à télécharger."
  echo "Usage: docker run ... sra-downloader <SRA_ID>"
  exit 1
fi

SRA_ID="$1"
OUTDIR=/data

echo "🚀 Téléchargement direct de $SRA_ID en FASTQ..."
fasterq-dump "$SRA_ID" -O "$OUTDIR" --split-files --progress | head -n 4000 > /fastq/${SRA_ID}_subset.fastq


echo "✅ Téléchargement terminé. Fichiers enregistrés dans $OUTDIR"
