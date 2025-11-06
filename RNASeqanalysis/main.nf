nextflow.enable.dsl=2

include { GET_SRR } from "./processes/GET_SRR/"
include { DOWNLOAD_FASTQ } from "./processes/DOWNLOAD_FASTQ/"
include { TRIM_SEQUENCE } from "./processes/TRIM_SEQUENCE/"
include { GET_REF_GENOME } from "./processes/GET_REF_GENOME/main" 
include { FEATURECOUNTS } from "./processes/FEATURECOUNTS/main" 

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

    // Trim the downloaded FASTQ files
    ch_trimmed_seq = TRIM_SEQUENCE(ch_fastq_files)

    FEATURECOUNTS(ch_trimmed_seq)

    // Get the reference genome
    ch_ref_genome = GET_REF_GENOME(params.ref_genome).ref_genome_file

}
