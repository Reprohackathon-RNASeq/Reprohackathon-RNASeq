process CREATE_COLDATA {
    publishDir "results/coldata", mode: 'copy', overwrite: true

    input:
    tuple path(geo_id), path(sra_data_file), path(script_R)
    
    output:
    path 'coldata.tsv', emit: coldata_file

    script:
    """
    Rscript ${script_R} ${geo_id} ${sra_data_file} coldata.tsv
    """

}