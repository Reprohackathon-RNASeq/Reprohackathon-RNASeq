process QUALITY_FASTQ {
  publishDir "results/quality/fastq", mode: 'copy', overwrite: true

  input:
    path fastq_files

  output:
    path "fastqc_report/*"

  script:
    """
    mkdir fastqc_report
    fastqc -o fastqc_report -f fastq ${fastq_files}
    """
}
