process DESeq2 {
    publishDir "results/traitement", mode: 'copy', overwrite: true

    input:
        path counted

    output:
        path "traitement", emit: traitement

    script:
        """
        mkdir -p traitement
        Rscript mon_script.R ${counted} traitement/deseq2_results.txt
        """
}