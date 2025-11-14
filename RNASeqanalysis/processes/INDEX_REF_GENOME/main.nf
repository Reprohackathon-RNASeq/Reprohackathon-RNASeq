process INDEX_REF_GENOME {
  cpus = 2
  
  publishDir "results/index", mode: 'copy', overwrite: true

  input:
    path ref_genome

  output:
    path "index_ref_genome", emit: indexed_genome

  script:
    """
    # Create output directory for the index
    mkdir -p index_ref_genome
    # apply bowtie-build to index the reference genome
    bowtie-build ${ref_genome} index_ref_genome/indexed_ref_genome
    """
}