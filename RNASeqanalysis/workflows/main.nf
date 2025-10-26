#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process GET_SRR {
    container 'mariemeier/sra-toolkit:latest'

    input:
    val sra_project

    output:
    path 'SRR_list.txt', emit: srr_list

    script:
    """
    esearch -db sra -query ${sra_project} | efetch -format runinfo | cut -d',' -f1 | grep SRR > SRR_list.txt
    """
}

// =======================
// PROCESSUS : TÉLÉCHARGEMENT SRA
// =======================
process DOWNLOAD_FASTQ {
    container 'mariemeier/sra-toolkit:latest'

    input:
    val sra_id

    output:
    path "${sra_id}*.fastq"  // tous les FASTQ générés
    publishDir "results/fastq", mode: 'copy'

    script:
    """
    echo "Téléchargement FASTQ pour ${sra_id}"
    fasterq-dump ${sra_id} -O ./ --split-files --progress
    echo "Téléchargement terminé : ${sra_id}"
    """
}


process DOWNLOAD_SRA {
    container 'mariemeier/sra-toolkit:latest'

    input:
    val sra_id

    output:
    path "${sra_id}"
    publishDir "results/raw_sra", mode: 'copy'

    script:
    """
    echo "Démarrage du téléchargement pour : ${sra_id}"
    prefetch ${sra_id}
    echo "Fichier SRA téléchargé dans le dossier ${sra_id}/"
    """
}

// =======================
// WORKFLOW
// =======================
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

    // Télécharger directement les FASTQ
    DOWNLOAD_FASTQ(ch_sra_ids)
}
