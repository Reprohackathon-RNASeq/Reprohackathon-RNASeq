process COUNT_READS {
  publishDir "results/count_matrix", mode: 'copy', overwrite: true

  input:
    path mapped_reads

  output:
    path "*_counts.txt"

  script:
    """
    featureCounts -o _counts.fq.gz -t gene -g ID -s 1 

    """
}

featureCounts --extraAttributes Name -t gene -g ID -F GTF -T <#CPUS> -a <GFF> -o counts.txt <BAM FILES>