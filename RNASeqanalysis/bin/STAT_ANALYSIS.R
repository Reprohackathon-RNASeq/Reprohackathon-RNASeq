#!/usr/bin/env Rscript
#chmod +x bin/STAT_ANALYSIS.R

###############################################################################
# STAT ANALYSIS
###############################################################################

# Load required libraries
library(DESeq2) 
library(KEGGREST)

# Expect 7 arguments: coldata file, count matrix file, gene map file, output matrix stat file, output MA plot file, output PCA plot file, output translation plot file
args <- commandArgs(trailingOnly = TRUE)
coldata_path <- args[1]
count_matrix_path <- args[2]
gene_map_path <- args[3]
output_matrix_stat <- args[4]
output_ma_plot <- args[5]
output_pca_plot <- args[6]
output_translation_plot <- args[7]

###############################################################################
# 1. LOAD AND HARMONIZE COL DATA
###############################################################################

# Load coldata: First column (SRR_ID) must become row names
coldata <- read.table(coldata_path, header = TRUE, sep = "\t", 
                      row.names = 1, stringsAsFactors = TRUE)

# Convert IP and control to factors.
coldata$Condition <- factor(coldata$Condition)
# Set 'control' as the reference level for log2FC calculation
coldata$Condition <- relevel(coldata$Condition, ref = "control") 

###############################################################################
# 2. LOAD, CLEAN, AND HARMONIZE COUNT MATRIX
###############################################################################

# Load the count matrix
# Skip the header comments
counts_matrix_raw <- read.table(count_matrix_path, header = TRUE, row.names = 1, 
                                sep = "\t", stringsAsFactors = FALSE)
counts_matrix_raw <- as.matrix(counts_matrix_raw)

# Filter out descriptive columns and keep only count columns
# We take columns 7 to 12 (the columns with SRR IDs).
counts_matrix <- counts_matrix_raw[, 6:11]

# Remove the suffix ".bam" from column names to match SRR IDs in coldata
colnames(counts_matrix) <- gsub("\\.bam$", "", colnames(counts_matrix))

# Now that names match (SRRxxxxxx == SRRxxxxxx), we reorder the matrix columns 
counts_matrix <- counts_matrix[, rownames(coldata)] 

# Ensure counts are numeric
counts_matrix <- apply(counts_matrix, 2, as.numeric) 
# Get back the rownames and colnames lost during apply
rownames(counts_matrix) <- rownames(counts_matrix_raw)
colnames(counts_matrix) <- rownames(coldata)

###############################################################################
# 3. DESEQ2 ANALYSIS
###############################################################################

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix), # Needs integer counts
  colData = coldata,
  design = ~ Condition 
)

# Run the DESeq2 analysis (Normalization, Dispersion estimation, Wald test)
# The article mentioned using default parameters (which this call does)
dds <- DESeq(dds)

# Extract results (IP vs control)
res <- results(dds, contrast = c("Condition", "IP", "control"), alpha = 0.05)


###############################################################################
# 4. STATISTICS MATRIX EXPORT
###############################################################################

# Convert results to a clean data frame
res_df_clean <- as.data.frame(res)

# Filter out rows where Log2FC or padj are NA
res_df_clean <- res_df_clean[complete.cases(res_df_clean[, c("log2FoldChange", "padj", "baseMean")]), ]
# Add Locus_Tag column from row names
res_df_clean$Locus_Tag <- rownames(res_df_clean)
# Remove "gene-" prefix from Locus_Tag
res_df_clean$Locus_Tag <- gsub("^gene-", "", res_df_clean$Locus_Tag)
# Reorder columns to have Locus_Tag first
res_df_clean <- res_df_clean[, c("Locus_Tag", setdiff(names(res_df_clean), "Locus_Tag"))]

# Save the statistics matrix to a TSV file
write.table(
    x = res_df_clean, 
    file = output_matrix_stat,
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE,
    col.names = TRUE
)

###############################################################################
# 5. MA-PLOT GENERATION
###############################################################################
# Prepare data for the plot (Define M and A)
M <- res_df_clean$log2FoldChange
A <- log2(res_df_clean$baseMean) # Convert average expression to log2 (MA-plot format)

# Define thresholds and colors
padj_cutoff <- 0.05 

# Define colors with transparency
red_transp <- rgb(1, 0, 0, 0.5)    
black_transp <- rgb(0, 0, 0, 0.5)
point_colors_transp <- ifelse(res_df_clean$padj < padj_cutoff, red_transp, black_transp)

# Open the PNG graphics device
png(output_ma_plot, width = 850, height = 800, res = 150, type = "cairo")

# MA-plot with custom grid and styles
par(
  mar = c(5, 5, 4, 2) + 0.1, # Margins
  font.lab = 1, # Bold axis titles
  col.axis = "gray30", # Color of axis labels
  col.lab = "black", # Color of axis titles
  fg = "gray30", # Color of axis borders
  cex.axis = 0.6, # Size of axis labels
  cex.lab = 0.8, # Size of axis titles
  mgp = c(0.9, 0.3, 0) # Distance of axis titles
)

plot(
  x = A, 
  y = M, 
  type = "n", # Do not draw points (to draw the grid first)
  ylim = c(-4, 4), 
  xlab = "Mean of normalized counts", 
  ylab = expression(log[2]~"fold change"),
  xaxt = 'n', # Do not draw the default X axis
  yaxt = 'n', # Do not draw the default Y axis
  bty = "n", # Remove the box around the plot
  
  # Aff fine grid in GGPLOT2 style
  panel.first = {
    # Set the background color of the inner panel (light gray)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
         col = "gray90", border = NA)
    
    # Add white grid lines
    abline(h = seq(-4, 4, by = 1), col = "white", lty = 1) # Horizontal white lines
    abline(v = log2(c(1, 10, 100, 1000, 10000, 100000)), col = "white", lty = 1) # Vertical white lines
  }
)

# Axes customization
breaks_y <- seq(-4, 4, by = 2)

axis(
  side = 2, # Axe Y
  at = breaks_y, # Labels positions
  labels = breaks_y,
  las = 1, # Horizontal labels
  lwd = 0, # Remove axis line 
  lwd.ticks = 0.5, # Thickness of small ticks
  tcl = -0.2, # Reduce length of ticks
  col.axis = "black", # Color of label text
  #mgp = c(2.2, 0.5, 0)
)

breaks_log2 <- log2(c(1, 100, 10000)) 

axis(
  side = 1, # Axe X
  at = breaks_log2, 
  labels = c(expression(10^0), expression(10^2), expression(10^4)),
  lwd = 0,            
    lwd.ticks = 1,
  tcl = -0.2,
  #cex.axis = 0.8,
  col.axis = "black",
  mgp = c(0.8, 0.1, 0) # Distance of axis titles
)

# Plot the points with colors and shapes
points(
  x = A[M < 4 & M > -4], 
  y = M[M < 4 & M > -4], 
  col = point_colors_transp[M < 4 & M > -4], 
  pch = 20, # Solid circle
  cex = 0.5 # Point size
)

points(
  x = A[M >= 4], 
  y = rep(4, length(M[M >= 4])), 
  col = point_colors_transp[M >= 4], 
  pch = 24, # Upward triangle
  bg = point_colors_transp[M >= 4], # Fill color
  cex = 0.4  # Point size
)

points(
  x = A[M <= -4], 
  y = rep(-4, length(M[M <= -4])), 
  col = point_colors_transp[M <= -4], 
  pch = 25, # Downward triangle
  bg = point_colors_transp[M <= -4], # Fill color
  cex = 0.4  # Point size
)

# Add the reference line
abline(h = 0, col = "gray50", lty = 2)

###############################################################################
# 6. PCA-PLOT GENERATION
###############################################################################

# VST (Variance Stabilizing Transformation) 
# 'blind = FALSE' to uses design (~ Condition) to improve transformation by taking the model into account
vsd <- vst(dds, blind = FALSE)

# PCA data extraction with plotPCA() (from deseq2)
data_pca <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)

# % of variance explained
percentVar <- round(100 * attr(data_pca, "percentVar"))

# PCA plot
png(output_pca_plot, width = 850, height = 800, res = 150, type = "cairo")

par(mgp = c(0.9, 0.3, 0))

plot(
  data_pca$PC1, data_pca$PC2,
  col = as.numeric(data_pca$Condition) + 1,
  pch = 19,
  main = "PCA of Samples (VST-transformed Data)",
  xlab = paste0("PC1: ", percentVar[1], "% variance"),
  ylab = paste0("PC2: ", percentVar[2], "% variance"),
  xaxt = "n",   # removing axes to redraw more cleanly
  yaxt = "n",
  bty = "n",
  
  # grey grid
  panel.first = {
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "gray90", border = NA)
    
    # white lines in the grid
    abline(h = pretty(data_pca$PC2, n = 10), col = "white", lty = 1)
    abline(v = pretty(data_pca$PC1, n = 10), col = "white", lty = 1)
  }
)

# Position the text to the left for the control group and to the right for the IP
pos_vec <- ifelse(data_pca$Condition == "control", 2, 4)

text(
  data_pca$PC1,
  data_pca$PC2,
  labels = rownames(data_pca),
  pos = pos_vec,
  cex = 0.5
)

# adding the legend
legend("bottomleft", 
       legend = levels(data_pca$Condition), 
       col = 1:length(levels(data_pca$Condition)) + 1, 
       pch = 19)

###############################################################################
# 7. KEGG GENE EXTRACTION
###############################################################################

# Traduction (sao03010 - Ribosome, facteurs)
translation_kegg_links <- KEGGREST::keggLink("sao", "sao03010")
translation2_kegg_links <- KEGGREST::keggLink("sao", "br:sao03012")

# Métabolisme des Acides Aminés (sao00970)
aa_metab_kegg_links <- KEGGREST::keggLink("sao", "sao00970")
# Aminoacyl-tRNA synthetases (sao01613)

ko_synthetases <- c(
    # Phe, Leu, Ile, Met, Val
    "K01889", "K01890", "K01869", "K01870", "K01874", "K01873",
    # Ser, Pro, Thr, Ala, Tyr, Gln
    "K01875", "K01881", "K01868", "K01872", "K01866", "K01886",
    # Glu (Glu/Pro), Lys, Pyl
    "K01885", "K09698", "K04566", "K04567", "K11627",
    # His, Asn, Asp
    "K01892", "K01893", "K01876", "K09759", "K01884",
    # Cys, Trp, Arg
    "K01883", "K01867", "K01887", "K01880",
    # Gly, Sep, (Inconnu)
    "K01878", "K01879", "K14164", "K07587", "K01894" 
)

aa_trna_ssynthetases_links <- KEGGREST::keggLink("sao", ko_synthetases)

# Traduction
translation_tags <- toupper(gsub("sao:", "", translation_kegg_links))
translation2_tags <- toupper(gsub("sao:", "", translation2_kegg_links))
# Métabolisme des AA
aa_metab_tags <- toupper(gsub("sao:", "", aa_metab_kegg_links))
# Aminoacyl-tRNA synthetases
aa_trna_ssynthetases_tags <- toupper(gsub("sao:", "", aa_trna_ssynthetases_links))
# Combine les gènes
genes_of_interest <- unique(c(translation_tags, aa_metab_tags, translation2_tags, aa_trna_ssynthetases_tags))

# Garder uniquement les gènes de traduction
res_translation <- res_df_clean[res_df_clean$Locus_Tag %in% genes_of_interest, ]


cat("Nombre de gènes dans sao03010 :", length(translation2_tags), "\n")
cat("Nombre de gènes dans sao03012 :", length(translation2_tags), "\n")
cat("Nombre de gènes dans sao00970 :", length(genes_of_interest), "\n")
cat("Nombre de gènes dans sao01613 :", length(aa_trna_ssynthetases_tags), "\n")
cat("Nombre de gènes dans res_translation :", nrow(res_translation), "\n")

###############################################################################
# 8. GENES TO PLOT
###############################################################################
gene_map <- read.csv(gene_map_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
cat(colnames(gene_map))

# Fusionner res_translation et gene_map par Locus_Tag
res_translation <- merge(
  res_translation, 
  gene_map, 
  by.x = "Locus_Tag",
  by.y = "Locus_Tag",
  all.x = TRUE  # garde tous les gènes de res_translation
)

# res_labeled$gene_name contient les noms pour les gènes présents dans gene_map
# Les autres sont NA

# Vérifier les gènes qui auront un label
cat("Gènes qui seront labellisés :\n")
cat(res_translation$gene_name[!is.na(res_translation$gene_name)], sep = "\n")

#gene_map <- read.table(gene_map_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Subset to genes of interest for labelling
#genes_interest_subset <- c("frr", "infA", "tsf", "infB", "infC", "pth")
#gene_map_filtered <- gene_map[gene_map$'pan gene symbol' %in% genes_interest_subset, ]

# Create labels vector
# locus_to_label <- res_translation$Locus_Tag %in% gene_map_filtered$'locus tag'
# labels_to_plot <- rep(NA, nrow(res_translation))
# labels_to_plot[locus_to_label] <- gene_map_filtered$`pan gene symbol`[match(
#   res_translation$Locus_Tag[locus_to_label], 
#   gene_map_filtered$'locus tag'
# )]

###############################################################################
# 9. TRANSLATION GENES MA-PLOT
###############################################################################


# Données MA
M_t <- res_translation$log2FoldChange
A_t <- log2(res_translation$baseMean)

# Couleurs (gènes significatifs/non significatifs)
point_colors_t <- ifelse(res_translation$padj < padj_cutoff, "red", "grey50")

# Ouvrir PNG
png(output_translation_plot, width = 700, height = 700, res = 150, type = "cairo")

par(
  mar = c(5,5,2,2) + 0.1,    # marges
  font.lab = 2,               # axes en gras
  font.axis = 2,              # labels axes en gras
  cex.axis = 0.9,             # taille des labels axes
  cex.lab = 1,                # taille titres axes
  col.axis = "black",
  col.lab = "black",
  mgp = c(1.5, 0.5, 0),      # distance titres axes
  lwd = 2                     # épaisseur du cadre et axes
)

# Plot vide
plot(
  x = A_t,
  y = M_t,
  type = "n",
  xlim = c(0, 20),
  ylim = c(-6, 5),
  xlab = expression(log[2]~"base Mean"),
  ylab = expression(log[2]~"fold change"),
  xaxt = "n",
  yaxt = "n",
  bty = "o"  # cadre noir
)

# Ligne pointillée horizontale
abline(h = 0, col = "black", lty = 2)


# Axes personnalisés, plus épais et rapprochés
axis(side = 1, at = seq(0, 20, by = 2), lwd = 2, lwd.ticks = 2, tcl = -0.3)
axis(side = 2, at = seq(-6, 5, by = 1), lwd = 2, lwd.ticks = 2, tcl = -0.3, las = 1)

# Points
points(
  x = A_t[M_t < 5 & M_t > -6],
  y = M_t[M_t < 5 & M_t > -6],
  col = point_colors_t[M_t < 5 & M_t > -6],
  pch = 20,
  cex = 0.6
)

gene_positions <- list(
  pth   = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="pth"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="pth"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="pth"]) - 1.5,
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="pth"],
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="pth"]) - 2,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="pth"]),

  frr   = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="frr"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="frr"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="frr"]),
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="frr"] + 1.5,
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="frr"]) - 0.5,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="frr"] + 1.8),

  infA  = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="infA"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="infA"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="infA"]),
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="infA"] + 1.5,
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="infA"]) - 0.5,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="infA"] + 1.8),
               
  tsf   = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="tsf"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="tsf"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="tsf"]) + 1.5,
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="tsf"] + 1.5,
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="tsf"]) + 2,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="tsf"] + 1.8),

  infC  = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="infC"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="infC"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="infC"]) + 1.5,
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="infC"] - 1.5,
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="infC"]) + 2,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="infC"] - 1.8),

  infB  = list(x0 = log2(res_translation$baseMean[res_translation$gene_name=="infB"]),
               y0 = res_translation$log2FoldChange[res_translation$gene_name=="infB"],
               x1 = log2(res_translation$baseMean[res_translation$gene_name=="infB"]) + 1.5,
               y1 = res_translation$log2FoldChange[res_translation$gene_name=="infB"] - 1.5,
               text_x = log2(res_translation$baseMean[res_translation$gene_name=="infB"]) + 2,
               text_y = res_translation$log2FoldChange[res_translation$gene_name=="infB"] - 1.8)
)

# Dessiner les traits et labels
for(g in names(gene_positions)) {
  pos <- gene_positions[[g]]
  arrows(
    x0 = pos$x0, y0 = pos$y0, 
    x1 = pos$x1, y1 = pos$y1,
    length = 0,       # pas de flèche
    col = "black", 
    lwd = 2           # trait plus gras
  )
  text(
    x = pos$text_x, y = pos$text_y, 
    labels = g, 
    cex = 1.1, 
    font = 2, 
    col = "black"
  )
}

legend("bottomleft", 
       legend = c("Significant", "Non Significant"),
       col = c("red", "black"),
       pch = 20,
       pt.cex = 1,
       cex = 0.8,      # taille de la légende
       bty = "n")       # pas de bordure

legend("bottomright", 
       legend = c("AA-tRNA synthetases"), 
       col = c("black"), 
       pch = 1,         
       pt.cex = 0.9,
       bty = "n",
       cex = 0.8)


res_translation$is_aa_trna_ssynthetases <- res_translation$Locus_Tag %in% aa_trna_ssynthetases_tags
aa_trna_ssynthetases_points <- res_translation[res_translation$is_aa_trna_ssynthetases & 
                                   res_translation$log2FoldChange < 5 & 
                                   res_translation$log2FoldChange > -6, ]

points(
  x = log2(aa_trna_ssynthetases_points$baseMean),
  y = aa_trna_ssynthetases_points$log2FoldChange,
  pch = 1,       
  col = "black", 
  lwd = 2,      
  cex = 0.9
)
# Close the graphics device
dev.off()
