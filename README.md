# 🧫 ReproHackaton - Analyse RNA-Seq de Staphylococcus aureus persisters

## 📋 Description du Projet

Ce projet s'inscrit dans le cadre d'un **ReproHackaton** visant à reproduire l'analyse RNA-Seq de l'article :

> **"Intracellular Staphylococcus aureus persisters upon antibiotic exposure"**  
> *Nature Communications* (2020) 11:2200  
> DOI: [10.1038/s41467-020-15966-7](https://doi.org/10.1038/s41467-020-15966-7)

### 🎯 Objectif Scientifique
Analyser les gènes différentiellement exprimés chez les **persisters intracellulaires** de *Staphylococcus aureus* après exposition aux antibiotiques, comparés aux bactéries contrôles.

## 🧬 Contexte Biologique

Les **persisters bactériens** sont des variants phénotypiques qui :
- Entrent dans un état **non-divisant transitoire**
- Développent une **tolérance aux antibiotiques**
- Contribuent aux **infections chroniques** et aux **échecs thérapeutiques**
- Constituent un **réservoir pour les rechutes**

## 🛠️ Workflow Implementé

### 📊 Pipeline d'Analyse RNA-Seq

Téléchargement FASTQ (SRA) → Trimming → Mapping → Comptage → Analyse Différentielle

### 🔧 Outils Utilisés

| Étape | Outil | Version | Conteneur |
|-------|-------|---------|-----------|
| Téléchargement | `fasterq-dump` | 3.0.7 | ✅ |
| Qualité | `FastQC` | 0.12.1 | ✅ |
| Trimming | `TrimGalore` | 0.6.10 | ✅ |
| Mapping | `Bowtie2` | 2.5.1 | ✅ |
| Comptage | `featureCounts` | 2.0.3 | ✅ |
| Analyse DESeq2 | `R` + `DESeq2` | 4.3.1 | ✅ |

## 🗂️ Structure du Projet

projet_reprohackaton/
├── 📁 workflows/
│ ├── snakefile.smk # Workflow Snakemake principal
│ └── config.yaml # Configuration
├── 📁 containers/
│ ├── Dockerfile.fastq # Conteneur téléchargement
│ ├── Dockerfile.rnaseq # Conteneur analyse RNA-Seq
│ └── Dockerfile.deseq2 # Conteneur analyse statistique
├── 📁 scripts/
│ ├── download_sra.py # Script téléchargement SRA
│ └── deseq2_analysis.R # Analyse différentielle
├── 📁 data/
│ ├── raw/ # Données brutes FASTQ
│ ├── processed/ # Données traitées
│ └── reference/ # Génome de référence
├── 📁 results/
│ ├── counts/ # Matrices de comptage
│ ├── differential/ # Résultats DESeq2
│ └── figures/ # Graphiques et plots
├── 📄 README.md # Ce fichier
├── 📄 run.sh # Script d'exécution
└── 📄 environment.yml # Environnement Conda


## 🚀 Installation et Exécution

### Prérequis
```bash
# 1. Installer Conda/Mamba
conda install -c conda-forge mamba

# 2. Installer Snakemake
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake

# 3. Cloner le repository
git clone https://github.com/Reprohackaton-RNASeq/Reprohackaton-RNASeq.git
cd projet_reprohackaton
```
### Exécution rapide 
```bash
# Lancer tout le workflow
./run.sh

# Ou exécuter manuellement
snakemake --use-conda --cores 8 --resources mem_mb=16000
```
### Avec docker
```bash
# Build des conteneurs
docker build -t reprohackaton-rnaseq -f containers/Dockerfile.rnaseq .

# Exécution
docker run -v $(pwd)/data:/data reprohackaton-rnaseq
```

## 📊 Résultats Attendus

### 🔍 Sorties Principales

Matrice de comptage des reads par gène
Liste des gènes différentiellement exprimés (padj < 0.05)
MA-plots et volcano plots de visualisation
Heatmaps des profils d'expression
Enrichissement fonctionnel (voies KEGG)

### 🧪 Validation Reproductible

✅ Tous les outils conteneurisés
✅ Environnements reproductibles (Conda)
✅ Code versionné (Git)
✅ Documentation complète
👥 Équipe

## Étudiants M2 AMI2B - ReproHackaton 2025

[Donatien WALLAERT]
[Tom GORTANA]
[Marie MEIER]
[Tom BELLIVIER]

## 📚 Références

Article principal : Peyrusson et al. (2020) Nature Communications
Workflow : Snakemake best practices
Analyse RNA-Seq : DESeq2 vignette
Génome référence : S. aureus NCTC 8325 (CP000253.1)


## 📄 License

MIT License - Voir le fichier LICENSE pour plus de détails.