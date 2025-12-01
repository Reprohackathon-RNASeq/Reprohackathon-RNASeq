process STAT_ANALYSIS {
  publishDir "results/stat_analysis", mode: 'copy', overwrite: true

  input:
    tuple path(coldata), path(count_matrix)
    path gene_names

  output:
    path "matrix_stat.tsv", emit: matrix_stat_file
    path "ma_plot.png", emit: ma_plot_file
    path "pca_plot.png", emit: pca_plot_file
    path "translation_plot.png", emit: translation_plot_file

  script:
    """
    STAT_ANALYSIS.R ${coldata} ${count_matrix} ${gene_names} matrix_stat.tsv ma_plot.png pca_plot.png translation_plot.png
    """
}
