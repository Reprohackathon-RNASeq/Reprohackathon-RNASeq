process MAPPING_BOWTIE {
  publishDir "results/mapping", mode: 'copy', overwrite: true

  input:
    path indexed_genome
    path trimmed_fastq

  output:
    path "*.sam", emit: mapped_reads

  script:
    """
    bowtie -x ${indexed_genome}/indexed_ref_genome -U ${trimmed_fastq} -S \$(basename ${trimmed_fastq} .fq.gz).sam
    """
}