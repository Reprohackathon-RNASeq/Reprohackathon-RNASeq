process GET_REF_GENOME {
    publishDir "data/ref_genome", mode: 'copy', overwrite: true

    input:
    val ref_genome

    output:
    path '', emit: ref_genome_file

    script:
    """
    esearch -db nucleotide -query $ref_genome | efetch -format fasta > reference.fasta
    """
}