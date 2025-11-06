process TRIM_SEQUENCE {
  publishDir "results/trimmed", mode: 'copy', overwrite: true

  input:
    path ch_fastq_files

  output:
    path "*_trimmed.fq.gz"

  script:
    """
    # Remove .gz extension
    BASENAME=\$(basename ${ch_fastq_files} .gz)

    # Remove .fastq extension to get SRA ID
    SRA_ID=\${BASENAME%.fastq}

    trim_galore -q 20 --phred33 --length 25 \${SRA_ID}.fastq.gz 

    # Rename output file
    # mv \${SRA_ID}.fastq.gz_trimmed.fq.gz \${SRA_ID}_trimmed.fq.gz
    """
}