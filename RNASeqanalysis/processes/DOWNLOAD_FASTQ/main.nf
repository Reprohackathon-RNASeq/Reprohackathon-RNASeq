process DOWNLOAD_FASTQ {
    cpus = 3 
    
    publishDir "data/fastq", mode: 'copy', overwrite: true
    
    tag "$sra_id"

    input:
    val sra_id
    
    output:
    path "${sra_id}.fastq.gz", emit: fastq_files

    script:
    """
    echo "Downloading FASTQ for ${sra_id}"
    prefetch ${sra_id}
    
    if [ "${params.test}" == "true" ]; then
        # Mode Test : Tronque et compresse.
        echo "Running in TEST MODE — keeping only first 10,000 reads"
        
        fasterq-dump --threads ${task.cpus} -t . --progress ${sra_id} 
        gzip ${sra_id}.fastq
        
        zcat ${sra_id}.fastq.gz | head -n 10000 | gzip > ${sra_id}_test.fastq.gz
        mv ${sra_id}_test.fastq.gz ${sra_id}.fastq.gz
        
    else
        fasterq-dump --threads ${task.cpus} -t . --progress ${sra_id} 
        gzip ${sra_id}.fastq
    fi

    echo "DOWNLOADING completed: ${sra_id}"
    """
}
