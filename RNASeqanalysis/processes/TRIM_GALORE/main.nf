process TRIM_GALORE {
  publishDir "data/trimmed", mode: 'copy'
  input:
    tuple val(sample_id), path(fastq)

  output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz")

  script:
    """
    trim_galore -q 20 --phred33 --length 25 $fastq
    """
}