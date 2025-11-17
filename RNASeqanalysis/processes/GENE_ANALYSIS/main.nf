process GENE_ANALYSIS {
    publishDir "results/gene_analysis", mode: 'copy', overwrite: true

    input:
    path coldata_file
    path counts_file

    output:
    file 'coldata.tsv', emit: coldata_file

    script:
    """
    cat <<-EOF > CREATE_COLDATA.R

    ###############################################################################
    #GET COLDATA
    ###############################################################################

    library(GEOquery)
    library(dplyr)
    library(stringr)
    library(readr) 
    library(tibble)

    geo_accession <- "${geo_id}"
    sra_file_path <- "${sra_data_file}"

    # Table mapping SRR_ID AND GSM_ID
    mapping_table <- sra_run_table %>%
        dplyr::select(
            SRR_ID = Run, 
            GSM_ID = SampleName
        ) %>%
        dplyr::filter(!is.na(SRR_ID))

    # Download the GEO Series data
    gse <- getGEO(geo_accession, GSEMatrix = TRUE, AnnotGPL = FALSE)
    pData <- pData(gse[[1]])

    colData_raw <- data.frame(
        GSM_ID = rownames(pData),
        Condition = pData$title,
        stringsAsFactors = FALSE
        )

    # Removing the suffix from the SampleTitle column
    colData_raw$Condition <- str_replace(colData_raw$Condition, " replicate [0-9]+$", "")

    #Mapping
    colData_mapped <- colData_raw %>% left_join(mapping_table, by = "GSM_ID")

    colData <- colData_mapped %>%
        dplyr::select(
            SRR_ID,    
            Condition) %>%
        dplyr::relocate(SRR_ID, .before = Condition)
    
    # Final colData
    write.table(colData, file='coldata.tsv', quote=FALSE, sep='\t', row.names=FALSE)
    ###############################################################################

    EOF
    
    Rscript CREATE_COLDATA.R
    """
}