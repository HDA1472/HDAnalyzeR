utils::globalVariables(c("DAid", "Assay", "NPX", "Var1", "Var2", "Freq"))
#' Calculate the percentage of NAs in each column
#'
#' The function calculates the percentage of NAs in each column of the input dataframe.
#' It filters out the columns with 0% missing data and returns the rest in descending order.
#'
#' @param df (tibble). The input dataframe
#'
#' @return na_percentage (tibble). A tibble with the column names and the percentage of NAs in each column
#' @export
#'
#' @examples
#' na_percentages <- calc_na_percentage_col(example_metadata)
#' print(na_percentages)
calc_na_percentage_col <- function(df) {

  na_percentage <- df |>
    dplyr::summarise_all(~ round(sum(is.na(.)) / dplyr::n() * 100, 1)) |>
    tidyr::gather(key = "column", value = "na_percentage") |>
    dplyr::filter(na_percentage > 0) |>  # Filter out columns with no NAs
    dplyr::arrange(dplyr::desc(na_percentage))

  return(na_percentage)
}


#' Calculate the percentage of NAs in each row
#'
#' The function calculates the percentage of NAs in each row of the input dataframe.
#' It filters out the rows with 0% missing data and returns the rest in descending order.
#'
#' @param df (tibble). The input dataframe
#'
#' @return na_percentage (tibble). A tibble with the DAids and the percentage of NAs in each row
#' @export
#'
#' @examples
#' na_percentages <- calc_na_percentage_row(example_metadata)
#' print(na_percentages)
calc_na_percentage_row <- function(df) {

  na_percentage <- df |>
    dplyr::rowwise() |>
    dplyr::mutate(na_percentage = round(sum(is.na(dplyr::across(dplyr::everything())))/ncol(df) * 100, 1)) |>
    dplyr::ungroup() |>
    dplyr::filter(na_percentage > 0) |>
    dplyr::arrange(dplyr::desc(na_percentage)) |>
    dplyr::select(DAid, na_percentage)

  return(na_percentage)
}


#' Check normality of the data
#'
#' @param df (tibble). The input dataframe
#'
#' @return normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status
#' @export
#'
#' @examples
#' df <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX")
#' normality_results <- check_normality(df)
check_normality <- function(df) {

  future::plan(future::multicore)

  # Perform Shapiro-Wilk test in parallel
  p_values <- future.apply::future_lapply(df |> dplyr::select(-dplyr::any_of(c("DAid"))), function(column) {
    stats::shapiro.test(column)$p.value
  })

  adjusted_p_values <- stats::p.adjust(unlist(p_values), method = "BH")

  # Determine normality based on adjusted p-values
  normality <- adjusted_p_values > 0.05

  normality_results <- tibble::tibble(
    Protein = names(p_values),
    p_value = unlist(p_values),
    adjusted_p_value = adjusted_p_values,
    is_normal = normality
  )

  return(normality_results)
}


#' Create a correlation heatmap and report the protein-protein correlations above a threshold
#'
#' The function calculates the correlation matrix of the input dataframe.
#' It filters out the protein pairs with correlation values above or below the threshold.
#' It then creates a heatmap of the correlation matrix.
#'
#' @param df (tibble). The input dataframe
#' @param threshold (numeric). The reporting protein-protein correlation threshold. Default is 0.8
#'
#' @return A list containing the following elements:
#'   - cor_matrix (matrix). A matrix of protein-protein correlations
#'   - cor_results (tibble). A tibble with the filtered protein pairs and their correlation values
#'   - p (ggplot). A heatmap of protein-protein correlations
#' @export
#'
#' @examples
#' df <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX")
#' cor_results <- create_corr_heatmap(df, threshold = 0.7)
create_corr_heatmap <- function(df, threshold) {

  cor_matrix <- round(
    stats::cor(df |> dplyr::select(-dplyr::any_of(c("DAid"))),
               use = "pairwise.complete.obs",
               method = "pearson"),
    2
  )

  cor_long <- as.data.frame(as.table(cor_matrix), .name_repair = "minimal", stringsAsFactors = FALSE)

  cor_results <- cor_long |>
    dplyr::filter(Var1 != Var2) |>
    dplyr::filter(Freq > threshold | Freq < -threshold) |>
    dplyr::arrange(dplyr::desc(Freq)) |>
    dplyr::rename(Protein1 = Var1, Protein2 = Var2, Correlation = Freq)

  p <- tidyheatmaps::tidyheatmap(cor_long,
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

  return(list(cor_matrix = cor_matrix, cor_results = cor_results, p = p))
}


#' Summarize the quality control results of Olink data
#'
#' @param df (tibble). The input dataframe
#' @param wide (logical). Whether the input dataframe is in wide format. Default is TRUE
#' @param threshold (numeric). The reporting protein-protein correlation threshold. Default is 0.8
#'
#' @return A list containing the following elements:
#'   - na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column
#'   - na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row
#'   - normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status
#'   - cor_matrix (matrix). A matrix of protein-protein correlations
#'   - cor_results (tibble). A tibble with the filtered protein pairs and their correlation values
#'   - heatmap (ggplot). A heatmap of protein-protein correlations
#' @export
#'
#' @examples
#' qc_summary_data(example_data, wide = FALSE, threshold = 0.7)
qc_summary_data <- function(df, wide = T, threshold = 0.8) {

  if (isFALSE(wide)) {
    df <- df |>
      dplyr::select(DAid, Assay, NPX) |>
      tidyr::pivot_wider(names_from = "Assay", values_from = "NPX")
  }

  sample_n <- nrow(df)
  protein_n <- ncol(df) - 1
  na_percentage_col <- calc_na_percentage_col(df)
  na_percentage_row <- calc_na_percentage_row(df)
  normality_results <- check_normality(df)
  cor <- create_corr_heatmap(df, threshold = threshold)
  cor_matrix <- cor$cor_matrix
  cor_results <- cor$cor_results
  p <- cor$p

  return(
    list(
      na_percentage_col = na_percentage_col,
      na_percentage_row = na_percentage_row,
      normality_results = normality_results,
      cor_matrix = cor_matrix,
      cor_results = cor_results,
      heatmap = p
      )
    )
}
