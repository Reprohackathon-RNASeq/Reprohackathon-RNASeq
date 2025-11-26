process STAT_ANALYSIS {
  publishDir "results/plots", mode: 'copy', overwrite: true

  input:
    tuple path(coldata), path(count_matrix)

  output:
    path "ma_plot.png", emit: ma_plot_file
    path "pca_plot.png", emit: pca_plot_file

  script:
    """
    STAT_ANALYSIS.R ${coldata} ${count_matrix} ma_plot.png pca_plot.png
    """
}
