process MAPPING_BOWTIE {
  publishDir "results/mapping", mode: 'copy', overwrite: true

  input:
    tuple path(trimmed_fastq), path(indexed_genome)

  output:
    path "*.sam", emit: mapped_reads

  script:
    """
    # Unzip directly and map with Bowtie
    gunzip -c ${trimmed_fastq} | bowtie -S ./index_ref_genome/indexed_ref_genome - > ${trimmed_fastq.baseName.replace('.fq.gz','')}.sam
    """
}