#!/usr/bin/env Rscript
#chmod +x bin/RESULTS_COMPARISON.R

###############################################################################
# RESULTS COMPARISON - Base R Version
###############################################################################

# Load necessary libraries
library(dplyr)

# Expect 3 arguments: matrix stat file, results article file, output path
args <- commandArgs(trailingOnly = TRUE)
matrix_stat_path <- args[1]
results_article_path <- args[2]
output_path <- args[3]

###############################################################################
# 1. DATA PREPARATION AND MERGING
###############################################################################

# Load DESeq2 results
matrix_stat <- read.table(matrix_stat_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load reference article results
matrix_article <- read.table(results_article_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Join the two matrices on Locus_Tag
comparison_data <- matrix_stat %>%
  inner_join(matrix_article, c("Locus_Tag" = "Name"), suffix = c("", ".ref"))

# Convert Log2FC and Padj columns to numeric
comparison_data$log2FoldChange <- as.numeric(comparison_data$log2FoldChange)
comparison_data$log2FoldChange.ref <- as.numeric(comparison_data$log2FoldChange.ref)
comparison_data$padj <- as.numeric(comparison_data$padj)
comparison_data$padj.ref <- as.numeric(comparison_data$padj.ref)

# Remove rows with NA
comparison_data <- comparison_data[complete.cases(comparison_data[, c("log2FoldChange","log2FoldChange.ref")]), ]

###############################################################################
# 2. CORRELATION AND SIGNIFICANCE CATEGORIZATION
###############################################################################

# Pearson correlation
rho <- cor(comparison_data$log2FoldChange, comparison_data$log2FoldChange.ref, method = "pearson")

# Categorize significance
comparison_data$category <- "no"
comparison_data$category[comparison_data$padj < 0.05 & comparison_data$padj.ref < 0.05] <- "in both"
comparison_data$category[comparison_data$padj < 0.05 & comparison_data$padj.ref >= 0.05] <- "only in our data"
comparison_data$category[comparison_data$padj >= 0.05 & comparison_data$padj.ref < 0.05] <- "only in the article"

# Colors
col_vec <- c("in both" = "darkred", "only in our data" = "orange", 
             "only in the article" = "darkblue", "no" = "gray50")
point_colors <- col_vec[comparison_data$category]

###############################################################################
# 3. SCATTER PLOT (Base R)
###############################################################################

png(output_path, width = 7, height = 6, units = "in", res = 300, type = "cairo")

plot(
  comparison_data$log2FoldChange.ref, 
  comparison_data$log2FoldChange, 
  col = point_colors, 
  pch = 20,
  xlab = "Reference Log2FC",
  ylab = "Calculated Log2FC",
  main = paste0("Log2 Fold Changes Correlation\nR = ", round(rho, 3)),
  cex = 0.7
)

# Add grid
grid(col = "gray85", lty = 1)

# Dashed diagonal
abline(a = 0, b = 1, lty = 2, col = "black")

# ---- LEGENDS ----

# Title above categories
legend(
  "topright",
  legend = "padj < 0.05",
  bty = "n",
  cex = 0.9,
  text.font = 2
)

# Actual category legend
legend(
  "topright",
  legend = names(col_vec),
  col = col_vec,
  pch = 20,
  cex = 0.8,
  inset = c(0, 0.10)
)

dev.off()
