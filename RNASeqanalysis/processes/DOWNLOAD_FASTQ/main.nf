process DOWNLOAD_FASTQ {
    publishDir "data/fastq", mode: 'copy', overwrite: true
    
    tag "$sra_id"

    input:
    val sra_id
    
    output:
    path "${sra_id}.fastq.gz", emit: fastq_files

    script:
    """
    echo "Téléchargement FASTQ pour ${sra_id}"
    prefetch ${sra_id}
    fasterq-dump --threads ${task.cpus} --progress ${sra_id} 
    gzip ${sra_id}.fastq 

    # Si le mode test est activé, tronquer le fichier à 10 000 lectures
    if [ "${params.test}" = "true" ]; then
        echo "Running in TEST MODE — keeping only first 10,000 reads"
        zcat ${sra_id}.fastq.gz | head -n 40000 | gzip > ${sra_id}_test.fastq.gz
        mv ${sra_id}_test.fastq.gz ${sra_id}.fastq.gz
    fi
    
    echo "Téléchargement terminé : ${sra_id}"
    """
}