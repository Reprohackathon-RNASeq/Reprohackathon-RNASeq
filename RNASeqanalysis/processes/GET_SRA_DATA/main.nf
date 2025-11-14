process GET_SRA_DATA {
    publishDir "data/sra_data_raw", mode: 'copy', overwrite: true
    
    input:
    val sra_project 

    output:
    path 'sra_data_complete.csv', emit: sra_data_file 

    script:
    """ 
    esearch -db sra -query "${sra_project}" | efetch -format runinfo > "sra_data_complete.csv"
    """
}
