process INDEX_REF_GENOME {
  publishDir "results/index", mode: 'copy', overwrite: true

  input:
    path ref_genome

  output:
    path "index_ref_genome/*"

  script:
    """
    mkdir -p index_ref_genome
    bowtie-build ${ref_genome} index_ref_genome/indexed_ref_genome
    """
}