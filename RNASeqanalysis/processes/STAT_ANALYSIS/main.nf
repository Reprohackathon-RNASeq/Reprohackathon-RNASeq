process STAT_ANALYSIS {
  publishDir "results/ma_plot", mode: 'copy', overwrite: true

  input:
    tuple path(coldata), path(count_matrix), path(script_R)

  output:
    path "ma_plot.png", emit: ma_plot_file

  script:
    """
    Rscript ${script_R} ${coldata} ${count_matrix} ma_plot.png
    """
}
