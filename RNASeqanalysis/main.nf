nextflow.enable.dsl=2

include { GET_SRR } from "./processes/GET_SRR/"
include { DOWNLOAD_FASTQ } from "./processes/DOWNLOAD_FASTQ/"
include { QUALITY_FASTQ } from "./processes/QUALITY_FASTQ/"
include { TRIM_SEQUENCE } from "./processes/TRIM_SEQUENCE/"
include { GET_REF_GENOME } from "./processes/GET_REF_GENOME/" 
include { INDEX_REF_GENOME } from "./processes/INDEX_REF_GENOME/"
include { MAPPING_BOWTIE } from "./processes/MAPPING_BOWTIE/"
include { FEATURECOUNTS } from "./processes/FEATURECOUNTS/

params.sra_run = null
params.sra_project = null
params.ref_genome = null

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
    ch_ref_genome = GET_REF_GENOME(params.ref_genome).ref_genome_file

    // Index genome with bowtie
    ch_indexed_genome = INDEX_REF_GENOME(ch_ref_genome).indexed_genome

    // Map the trimmed sequences to the reference genome
    ch_mapping = MAPPING_BOWTIE(ch_trimmed_sequences.combine(ch_indexed_genome))

    // COOUNT FEATURES 
    counted = FEATURECOUNTS(ch_mapping)
}

