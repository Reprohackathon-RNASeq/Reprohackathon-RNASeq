nextflow.enable.dsl=2

include { GET_SRR } from "./processes/GET_SRR/"
include { DOWNLOAD_FASTQ } from "./processes/DOWNLOAD_FASTQ/"
include { QUALITY_FASTQ } from "./processes/QUALITY_FASTQ/"
include { TRIM_SEQUENCE } from "./processes/TRIM_SEQUENCE/"
include { GET_DATA_GENOME } from "./processes/GET_DATA_GENOME/" 
include { INDEX_REF_GENOME } from "./processes/INDEX_REF_GENOME/"
include { MAPPING_BOWTIE } from "./processes/MAPPING_BOWTIE/"
include { FEATURECOUNTS } from "./processes/FEATURECOUNTS/"
include { GET_GEO_TABLE } from "./processes/GET_GEO_TABLE/"
include { GET_SRA_DATA } from "./processes/GET_SRA_DATA/"
include { CREATE_COLDATA } from "./processes/CREATE_COLDATA/"
include { STAT_ANALYSIS } from "./processes/STAT_ANALYSIS/"

workflow {

    // Define the channel for SRA IDs based on the parameters provided (sra_run ou sra_project)
    if (params.sra_run) {
        ch_sra_ids = Channel.value(params.sra_run)

    } else if (params.sra_project) {
        ch_sra_ids = GET_SRR(params.sra_project).srr_list
            .map { file ->
                file.text.readLines().collect { it.trim() }
            }
            .flatten()
    } else {
        error "You must provide either --sra_run or --sra_project."
    }

    // Download the FASTQ files for all SRA IDs
    ch_fastq_files = DOWNLOAD_FASTQ(ch_sra_ids).fastq_files 

    // Perform quality check on the downloaded FASTQ files
    //QUALITY_FASTQ(ch_fastq_files)

    // Trim the downloaded FASTQ files
    ch_trimmed_sequences = TRIM_SEQUENCE(ch_fastq_files)

    // Get the reference genome
    ch_data_genome = GET_DATA_GENOME(params.ref_genome)
    ch_ref_genome = ch_data_genome.ref_genome_file
    ch_annotations = ch_data_genome.gff_gile

    // Index genome with bowtie
    ch_indexed_genome = INDEX_REF_GENOME(ch_ref_genome).indexed_genome

    // Map the trimmed sequences to the reference genome
    ch_mapping = MAPPING_BOWTIE(ch_trimmed_sequences.combine(ch_indexed_genome))

    // Generate count matrix using featureCounts
    ch_count_matrix = FEATURECOUNTS(ch_mapping.mapped_reads.collect(), ch_annotations).count_matrix

    // Get GEO ID from SRA Project
    if (params.sra_project) {
        ch_geo_table = GET_GEO_TABLE(params.geo_id).geo_table
        ch_sra_data = GET_SRA_DATA(params.sra_project).sra_data_file
    }

    // Create colData file
    ch_script_R_coldata = Channel.value(file("${projectDir}/scripts/CREATE_COLDATA.R"))
    ch_coldata_file = CREATE_COLDATA(ch_geo_table.combine(ch_sra_data).combine(ch_script_R_coldata)).coldata_file

    ch_script_R_analysis = Channel.value(file("${projectDir}/scripts/STAT_ANALYSIS.R"))
    ch_ma_plot = STAT_ANALYSIS(ch_coldata_file.combine(ch_count_matrix).combine(ch_script_R_analysis)).ma_plot_file
}


