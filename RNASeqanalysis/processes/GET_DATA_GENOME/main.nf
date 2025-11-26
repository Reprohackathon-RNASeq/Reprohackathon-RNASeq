process GET_DATA_GENOME {
    publishDir "data/ref_genome", mode: 'copy', overwrite: true

    input:
    val ref_genome

    output:
    path 'reference.fasta', emit: ref_genome_file
    path 'annotation.gff', emit: gff_file

    script:
    """
    esearch -db nucleotide -query $ref_genome | efetch -format fasta > reference.fasta
    esearch -db nuccore -query $ref_genome | efetch -format gff3 > annotation.gff
    """
}
