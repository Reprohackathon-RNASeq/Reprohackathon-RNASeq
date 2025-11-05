process DOWNLOAD_FASTQ {
    tag "$sra_id"

    input:
    val sra_id
    
    output:
    path "${sra_id}.fastq.gz", emit: fastq_files

    publishDir "data/fastq", mode: 'copy', overwrite: true

    script:
    """
    echo "Téléchargement FASTQ pour ${sra_id}"
    prefetch ${sra_id}
    fasterq-dump --threads ${task.cpus} --progress ${sra_id} 
    gzip ${sra_id}.fastq 
    echo "Téléchargement terminé : ${sra_id}"
    """
}