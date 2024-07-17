utils::globalVariables(c("v1", "v2", "val", "x", "y"))
#' Cluster data
#'
#' This function clusters data and returns a data frame with the clustered data.
#' This data can be used to plot a heatmap with ggplot2.
#'
#' @param df (tibble). A data frame with the data to cluster.
#' @param distance_method (character). The distance method to use. Default is "euclidean".
#' @param clustering_method (character). The clustering method to use. Default is "ward.D2".
#' @param cluster_rows (logical). Whether to cluster rows. Default is TRUE.
#' @param cluster_cols (logical). Whether to cluster columns. Default is TRUE.
#' @param wide (logical). Whether the data is wide or long. Default is TRUE.
#'
#' @return (list). A list with the following elements:
#'  - clustered_data (tibble). A data frame with the clustered data.
#'  - hc_rows (hclust). The hierarchical clustering object for rows.
#'  - hc_cols (hclust). The hierarchical clustering object for columns.
#' @export
#'
#' @examples
#' clean_df <- example_data |> dplyr::select(DAid, Assay, NPX)
#' clustered_df <- cluster_data(clean_df, wide = FALSE)
cluster_data <- function(df,
                         distance_method = "euclidean",
                         clustering_method = "ward.D2",
                         cluster_rows = T,
                         cluster_cols = T,
                         wide = T) {
  if (isFALSE(wide)) {
    wide_data <- widen_data(df) |>
        tibble::column_to_rownames(var = names(df)[1])
  } else {
    wide_data <- df |>
      tibble::column_to_rownames(var = names(df)[1])
  }

    order_row <- rownames(wide_data)
    order_col <- colnames(wide_data)

    if(isTRUE(cluster_rows)) {
      hc_rows <- wide_data |>
        stats::dist(method = distance_method) |>
        stats::hclust(method = clustering_method)
      order1 <- hc_rows$labels[hc_rows$order]
    } else {
      hc_rows <- NULL
      order1 <- order_row
    }

    if(isTRUE(cluster_cols)) {
      hc_cols <- wide_data |>
        t() |>
        stats::dist(method = distance_method) |>
        stats::hclust(method = clustering_method)
      order2 <- hc_cols$labels[hc_cols$order]
    } else {
      hc_cols <- NULL
      order2 <- order_col
    }

    clustering_results <- wide_data |>
      tibble::rownames_to_column() |>
      dplyr::rename(v1 = 1) |>
      tidyr::gather(v2, val, -1) |>
      dplyr::rename(x = v1, y = v2, value = val) |>
      dplyr::mutate(x = factor(x, levels = order1),
                    y = factor(y, levels = order2)) |>
      dplyr::arrange(x, y)

    clustered_data <- tibble::as_tibble(clustering_results)
    return(list("clustered_data" = clustered_data, "hc_rows" = hc_rows, "hc_cols" = hc_cols))
}


#' Create dendrograms (under development - it should be inside a ggplot heatmap function)
#'
#' This function creates dendrograms for the rows and columns of the clustered data.
#'
#' @param hc_rows (hclust). The hierarchical clustering object for rows.
#' @param hc_cols (hclust). The hierarchical clustering object for columns.
#'
#' @return (list). A list with the following elements:
#'  - dend_rows (ggplot). A ggplot object with the dendrogram for rows.
#'  - dend_cols (ggplot). A ggplot object with the dendrogram for columns.
#' @keywords internal
create_dendrograms <- function(hc_rows, hc_cols) {
  if (!is.null(hc_rows)) {
    dend_rows <- ggdendro::ggdendrogram(hc_rows, labels = F, leaf_labels = F, rotate = T) +
      ggplot2::scale_y_reverse() +
      ggdendro::theme_dendro()
  }

  if (!is.null(hc_cols)) {
    dend_cols <- ggdendro::ggdendrogram(hc_cols, labels = F, leaf_labels = F) +
      ggdendro::theme_dendro()
  }
    return(list("dend_rows" = dend_rows, "dend_cols" = dend_cols))
}
