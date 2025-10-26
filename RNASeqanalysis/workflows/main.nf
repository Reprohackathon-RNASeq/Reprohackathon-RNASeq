nextflow.enable.dsl=2

include { GET_SRR } from "./processes/GET_SRR/"
include { DOWNLOAD_FASTQ } from "./processes/DOWNLOAD_FASTQ/"

params.sra_run = null
params.sra_project = null

workflow {

    def ch_sra_ids

    if (params.sra_run) {
        ch_sra_ids = Channel.value(params.sra_run)

    } else if (params.sra_project) {
        ch_sra_ids = GET_SRR(params.sra_project)
            .map { file ->
                file.text.readLines().collect { it.trim() }
            }
            .flatten()
    } else {
        error "Vous devez fournir soit --sra_run, soit --sra_project"
    }

    DOWNLOAD_FASTQ(ch_sra_ids)


}
