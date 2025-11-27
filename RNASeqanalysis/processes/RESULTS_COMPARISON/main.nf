process RESULTS_COMPARISON {
  publishDir "results/results_comparison", mode: 'copy', overwrite: true

  input:
    path(matrix_stat)
    path(results_article)

  output:
    path "comparison_plot.png", emit: comparison_plot_file

  script:
    """
    RESULTS_COMPARISON.R ${matrix_stat} ${results_article} comparison_plot.png
    """
}
