utils::globalVariables(c("Var1", "Var2", "Freq"))
#' Calculate the correlation matrix of the input dataframe
#'
#' @param x (vector or tibble). A numeric vector, matrix or data frame.
#' @param y (vector, tibble or NULL). A numeric vector, matrix or data frame with compatible dimensions with `x`. Default is NULL.
#' @param use (character). A character string. The method to use for computing correlations. Default is "pairwise.complete.obs".
#' @param method (character). A character string. The correlation method to use. Default is "pearson".
#'
#' @return cor_matrix (matrix). A matrix of protein-protein correlations.
#' @keywords internal
correlate <- function(x, y = NULL, use = "pairwise.complete.obs", method = "pearson") {

  cor_matrix <- round(
    stats::cor(x, y, use = "pairwise.complete.obs", method = "pearson"),
    2
  )

  return(cor_matrix)
}


#' Create a correlation heatmap and report the protein-protein correlations above a threshold
#'
#' The function calculates the correlation matrix of the input dataframe.
#' It filters out the protein pairs with correlation values above or below the threshold.
#' It then creates a heatmap of the correlation matrix.
#'
#' @param x (vector or tibble). A numeric vector, matrix or data frame.
#' @param y (vector, tibble or NULL). A numeric vector, matrix or data frame with compatible dimensions with `x`. Default is NULL.
#' @param use (character). A character string. The method to use for computing correlations. Default is "pairwise.complete.obs".
#' @param method (character). A character string. The correlation method to use. Default is "pearson".
#' @param threshold (numeric). The reporting protein-protein correlation threshold. Default is 0.8.
#' @param cluster_rows (logical). Whether to cluster the rows. Default is TRUE.
#' @param cluster_cols (logical). Whether to cluster the columns. Default is TRUE.
#'
#' @return A list containing the following elements:
#'   - cor_matrix (matrix). A matrix of protein-protein correlations.
#'   - cor_results (tibble). A tibble with the filtered protein pairs and their correlation values.
#'   - cor_plot (plot). A heatmap of protein-protein correlations.
#' @export
#'
#' @examples
#' df <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::select(-DAid)
#' cor_results <- create_corr_heatmap(df, threshold = 0.7)
create_corr_heatmap <- function(x,
                                y = NULL,
                                use = "pairwise.complete.obs",
                                method = "pearson",
                                threshold = 0.8,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE) {

  cor_matrix <- correlate(x,
                          y = NULL,
                          use = "pairwise.complete.obs",
                          method = "pearson")

  cor_long <- as.data.frame(as.table(cor_matrix), .name_repair = "minimal", stringsAsFactors = FALSE)

  cor_results <- cor_long |>
    dplyr::filter(Var1 != Var2) |>
    dplyr::filter(Freq > threshold | Freq < -threshold) |>
    dplyr::arrange(dplyr::desc(Freq)) |>
    dplyr::rename(Protein1 = Var1, Protein2 = Var2, Correlation = Freq)

  cor_plot <- tidyheatmaps::tidyheatmap(cor_long,
                                 rows = Var1,
                                 columns = Var2,
                                 values = Freq,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 show_selected_row_labels = c(""),
                                 show_selected_col_labels = c(""),
                                 color_legend_min = -1,
                                 color_legend_max = 1,
                                 treeheight_row = 20,
                                 treeheight_col = 20
  )

  return(list("cor_matrix" = cor_matrix, "cor_results" = cor_results, "cor_plot" = cor_plot))
}
