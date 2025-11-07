process MAPPING_BOWTIE {
  publishDir "results/mapping", mode: 'copy', overwrite: true

  input:
    tuple path(trimmed_fastq), path(indexed_genome)

  output:
    path "*.map", emit: mapped_reads

  script:
    """
    gunzip -c ${trimmed_fastq} > reads.fastq
    bowtie ./index_ref_genome/indexed_ref_genome reads.fastq > ${trimmed_fastq.baseName.replace('.fq.gz','')}.map
    """
}