process GET_SRR {
    tag "$sra_project"

    input:
    val sra_project

    output:
    path 'SRR_list.txt', emit: srr_list, overwrite: true

    script:
    """
    esearch -db sra -query ${sra_project} | efetch -format runinfo | cut -d',' -f1 | grep SRR > SRR_list.txt
    """
}
