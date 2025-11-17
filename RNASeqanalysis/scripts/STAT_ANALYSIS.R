###############################################################################
# STAT ANALYSIS
###############################################################################
# --- Contenu du fichier scripts/STAT_ANALYSIS.R (Compatible R 3.4.1) ---

# Charger UNIQUEMENT les librairies essentielles (qui doivent être la bonne version)
library(DESeq2)
library(GEOquery) # Si votre R 3.4.1 le supporte
# Note : Éviter dplyr, stringr, readr

# 1. LECTURE DES ARGUMENTS DE LA LIGNE DE COMMANDE (Non modifiée)
args <- commandArgs(trailingOnly = TRUE)
coldata_path <- args[1]
count_matrix_path <- args[2]
output_file <- args[3]

# --- 2. PRÉPARATION DES DONNÉES ---

# A. Charger et Harmoniser le colData (Utiliser read.table pour l'importation)
coldata <- read.table(coldata_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE)
# S'assurer que la colonne Condition est bien nommée (si votre colData final ne contient que la colonne Condition)

# Nettoyage des niveaux de facteur (Utilisation de gsub de base R)
coldata$Condition <- factor(gsub("Intracellular persister", "IP", coldata$Condition))
coldata$Condition <- relevel(coldata$Condition, ref = "control")

# B. Charger la Matrice de Comptage
# Utiliser le chemin et spécifier que les 7 premières colonnes sont des métadonnées.
counts_matrix_raw <- read.table(count_matrix_path, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# Retirer les colonnes de description (COLONNES 1 à 6 après l'ID de gène)
# On garde les colonnes après la 6e (les comptages), et on utilise as.matrix
counts_matrix <- counts_matrix_raw[, 7:ncol(counts_matrix_raw)]
counts_matrix <- as.matrix(counts_matrix)

# C. Harmonisation Finale des Noms de Colonnes (Utilisation de gsub de base R)
# Simplification du nom de colonne: SRRxxxxxx.bam -> SRRxxxxxx
colnames(counts_matrix) <- gsub("\\.bam$", "", colnames(counts_matrix))
rownames(coldata) <- gsub("_trimmed\\.fq\\.map$", "", rownames(coldata)) # Si besoin

# Réordonner la matrice selon l'ordre strict de colData
counts_matrix <- counts_matrix[, rownames(coldata)]


# --- 3. EXÉCUTION DE L'ANALYSE DESEQ2 ---

# Créer le DESeqDataSet (Round est nécessaire)
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = coldata,
    design = ~ Condition 
)

# Filtrage des gènes à faible expression (base R)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Exécuter l'analyse complète
dds <- DESeq(dds)

# Extraire les résultats
res <- results(dds, contrast = c("Condition", "IP", "control"), alpha = 0.05)


# --- 4. GÉNÉRATION DU MA-PLOT (Base R) ---

# Ouvrir le périphérique graphique PNG (Le nom du fichier de sortie)
png(output_file, width = 800, height = 800, res = 150)

# Créer le MA Plot
DESeq2::plotMA(res, 
               main = "Supp. Fig. 3.: MA-plot of complete RNA-seq dataset", 
               ylim = c(-4, 4))

# Fermer le périphérique graphique
dev.off()