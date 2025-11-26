process CREATE_COLDATA {
    publishDir "results/coldata", mode: 'copy', overwrite: true

    input:
    tuple path(geo_id), path(sra_data_file)
    
    output:
    path 'coldata.tsv', emit: coldata_file

    script:
    """
    CREATE_COLDATA.R ${geo_id} ${sra_data_file} coldata.tsv
    """

}