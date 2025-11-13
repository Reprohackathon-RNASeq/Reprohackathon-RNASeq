process GET_GEO_ID {
  input:
    val sra_project

    output:
    path geo_id, emit: geo_id_channel

    script:
    """
    GEO_ACCESSION=\$(esearch -db bioproject -query "$sra_project" | \
        elink -target gds | \
        esummary | \
        xtract -pattern DocumentSummary -element Accession)
        
    echo \$GEO_ACCESSION
    """
}