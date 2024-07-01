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
#' @return clustering_results (data.frame). A data frame with the clustered data.
#' @export
#'
#' @examples
#' clean_df <- example_data |> dplyr::select(DAid, Assay, NPX)
#' clustered_df <- cluster_data(clean_df, wide = FALSE)
cluster_data <- function(df, distance_method = "euclidean", clustering_method = "ward.D2",
                         cluster_rows = T, cluster_cols = T, wide = T) {

    wide_data <- widen_data(df, wide) |>
        tibble::column_to_rownames(var = names(df)[1])

    order_row <- rownames(wide_data)
    order_col <- colnames(wide_data)

    if(isTRUE(cluster_rows)) {
      order1 <- wide_data |>
        stats::dist(method = distance_method) |>
        stats::hclust(method = clustering_method) |>
        with(labels[order])
    } else {
      order1 <- order_row
    }

    if(isTRUE(cluster_cols)) {
      order2 <- wide_data |>
        t() |>
        stats::dist(method = distance_method) |>
        stats::hclust(method = clustering_method) |>
        with(labels[order])
    } else {
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

    return(clustering_results)
}
