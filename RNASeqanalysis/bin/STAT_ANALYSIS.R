#!/usr/bin/env Rscript
#chmod +x bin/STAT_ANALYSIS.R

###############################################################################
# STAT ANALYSIS
###############################################################################

# Load required libraries
library(DESeq2) 
library(KEGGREST)

# Expect 6 arguments: coldata file, count matrix file, gene map file, output matrix stat file, output MA plot file, output PCA plot file
args <- commandArgs(trailingOnly = TRUE)
coldata_path <- args[1]
count_matrix_path <- args[2]
gene_map_path <- args[3]
output_matrix_stat <- args[4]
output_ma_plot <- args[5]
output_pca_plot <- args[6]

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
# 4. GENE MAPPING AND STATISTICS MATRIX EXPORT
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
# 7. TRANSLATION PCA-PLOT GENERATION
###############################################################################

# Close the graphics device
dev.off()
