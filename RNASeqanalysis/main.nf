nextflow.enable.dsl=2

include { GET_SRR } from "./processes/GET_SRR/"
include { DOWNLOAD_FASTQ } from "./processes/DOWNLOAD_FASTQ/"
include { TRIM_GALORE } from './processes/TRIM_GALORE/main.nf'

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

    // 1. Téléchargement des fastq à partir des SRR diffusés
    def ch_fastq = DOWNLOAD_FASTQ(ch_sra_ids)

    // 2. Trim Galore sur tous les fastq téléchargés
    TRIM_GALORE(ch_fastq)


}


