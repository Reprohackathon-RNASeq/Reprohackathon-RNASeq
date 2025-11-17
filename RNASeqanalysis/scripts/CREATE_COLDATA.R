###############################################################################
# GET COLDATA
###############################################################################

# Load required libraries
library(GEOquery)
library(dplyr)
library(stringr)
library(readr) 
library(tibble)

# Expect 3 arguments: geo_id file, sra file, output file
args <- commandArgs(trailingOnly = TRUE)
geo_file_path <- args[1]
sra_file_path <- args[2]
output_file <- args[3]

# Read the simplified GEO table
colData_raw <- read_tsv(geo_file_path, col_names = c("GSM_ID", "Condition"))

# Read the SRA data table
sra_data_table <- read_csv(sra_file_path)

# Table mapping SRR_ID AND GSM_ID
mapping_table <- sra_data_table %>%
  dplyr::select(
    SRR_ID = Run, 
    GSM_ID = SampleName
  ) %>%
  dplyr::filter(!is.na(SRR_ID))


# Mapping
colData_mapped <- colData_raw %>% left_join(mapping_table, by = "GSM_ID")

# Select and reorder columns
colData <- colData_mapped %>% dplyr::select(SRR_ID, Condition)

# Final colData
write.table(colData, file=output_file, quote=FALSE, sep='\t', row.names=FALSE)
###############################################################################