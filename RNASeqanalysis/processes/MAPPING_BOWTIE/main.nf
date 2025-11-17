process MAPPING_BOWTIE {
  cpus = 3 
  
  publishDir "results/mapping", mode: 'copy', overwrite: true

  input:
    tuple path(trimmed_fastq), path(indexed_genome)

  output:
    path "*.bam", emit: mapped_reads

  script:
    def srr_id = trimmed_fastq.baseName.replaceAll('_trimmed.fq.gz', '').replaceAll('_trimmed.fq', '')
    """
    # Unzip directly and map with Bowtie
    gunzip -c ${trimmed_fastq} | \
    bowtie -S ./index_ref_genome/indexed_ref_genome - | \
    samtools view -@ ${task.cpus} -bS -o ${srr_id}.bam
    """
}