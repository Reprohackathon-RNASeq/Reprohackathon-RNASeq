process TRIM_SEQUENCE {
  publishDir "results/counted", mode: 'copy', overwrite: true

  input:
    path ch_trimmed_seq

  output:
    path "*__counted.fq.gz"

    // paramètres utilisés -t gene -g ID -s 1
  script:
    """
    featureCounts -o __counted.fq.gz -t gene -g ID -s 1 

    """
}