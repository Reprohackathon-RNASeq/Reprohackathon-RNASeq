process FEATURECOUNTS {
  cpus = 3 
  
  publishDir "results/count_matrix", mode: 'copy', overwrite: true

  input:
    path mapped_reads
    path gff_file

  output:
    path "counts.txt", emit: count_matrix // The final count matrix
    path "*.summary", emit: counts_summary // Summary file for each sample

  script:
    """
    featureCounts -t gene -g ID -s 1 -a ${gff_file} -o counts.txt ${mapped_reads}
    """
}
