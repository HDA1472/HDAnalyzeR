utils::globalVariables(c("DAid", "Assay", "NPX", "adj.P.Val"))
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

  adj.P.Val <- stats::p.adjust(unlist(p_values), method = "BH")

  # Determine normality based on adjusted p-values
  normality <- adj.P.Val > 0.05

  normality_results <- tibble::tibble(
    Protein = names(p_values),
    p_value = unlist(p_values),
    adj.P.Val = adj.P.Val,
    is_normal = normality
  ) |>
    dplyr::arrange((adj.P.Val))

  return(normality_results)
}


#' Print the summary of the quality control results
#'
#' The function prints the summary of the quality control results of the input dataframe.
#'
#' @param sample_n (numeric). The number of samples
#' @param var_n (numeric). The number of variables
#' @param na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column
#' @param na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row
#' @param normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status
#' @param cor_results (tibble). A tibble with the filtered protein pairs and their correlation values
#' @param heatmap (plot). A heatmap of protein-protein correlations
#' @param threshold (numeric). The reporting protein-protein correlation threshold
#'
#' @return NULL
#' @export
#'
#' @examples
#' summary_results <- qc_summary_data(example_data, wide = FALSE, threshold = 0.7, report = FALSE)
#' print_summary(summary_results$sample_n, summary_results$var_n, summary_results$na_percentage_col,
#'               summary_results$na_percentage_row, summary_results$normality_results,
#'               summary_results$cor_results, summary_results$heatmap, 0.7)
print_summary <- function(sample_n, var_n, na_percentage_col, na_percentage_row,
                          normality_results = F, cor_results = F, heatmap = F,  threshold = F) {

  print("Summary:")
  print("Note: In case of long output, only the first 10 rows are shown. To see the rest display the object with view()")
  print(paste0("Number of samples: ", sample_n))
  print(paste0("Number of variables: ", var_n))
  print("--------------------------------------")
  print("NA percentage in each column:")
  print(na_percentage_col)
  print("--------------------------------------")
  print("NA percentage in each row:")
  print(na_percentage_row)
  print("--------------------------------------")
  if (!isFALSE(normality_results)) {
    print("Normality test results:")
    print(normality_results)
    print("--------------------------------------")
  }
  if (!isFALSE(cor_results)) {
    print(paste0("Protein-protein correlations above ", threshold, ":"))
    print(cor_results)
    print("--------------------------------------")
  }
  if (!isFALSE(heatmap)) {
    print("Correlation heatmap:")
    print(heatmap)
    print("--------------------------------------")
  }

  invisible(NULL)
}


#' Summarize the quality control results of Olink data
#'
#' The function summarizes the quality control results of the input dataframe.
#' It can handles both long and wide dataframes
#'
#' @param df (tibble). The input dataframe
#' @param wide (logical). Whether the input dataframe is in wide format. Default is TRUE
#' @param threshold (numeric). The reporting protein-protein correlation threshold. Default is 0.8
#' @param report (logical). Whether to print the summary. Default is TRUE
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
qc_summary_data <- function(df, wide = T, threshold = 0.8, report = T) {

  df <- widen_data(df, wide)

  sample_n <- nrow(df)
  protein_n <- ncol(df) - 1
  na_percentage_col <- calc_na_percentage_col(df)
  na_percentage_row <- calc_na_percentage_row(df)
  normality_results <- check_normality(df)
  cor <- create_corr_heatmap(df |> dplyr::select(-dplyr::any_of(c("DAid"))),
                             threshold = threshold)
  cor_matrix <- cor$cor_matrix
  cor_results <- cor$cor_results
  p <- cor$p

  if (isTRUE(report)) {
    print_summary(sample_n, protein_n, na_percentage_col, na_percentage_row,
                  normality_results, cor_results, p, threshold)
  }

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


#' Summarize the quality control results of metadata
#'
#' The function summarizes the quality control results of the input dataframe.
#'
#' @param df (tibble). The input dataframe
#' @param report (logical). Whether to print the summary. Default is TRUE
#'
#' @return A list containing the following elements:
#'   - na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column
#'   - na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row
#' @export
#'
#' @examples
#' qc_summary_metadata(example_metadata)
qc_summary_metadata <- function(df, report = T) {

  sample_n <- nrow(df)
  var_n <- ncol(df) - 1
  na_percentage_col <- calc_na_percentage_col(df)
  na_percentage_row <- calc_na_percentage_row(df)

  if (isTRUE(report)) {
    print_summary(sample_n, var_n, na_percentage_col, na_percentage_row)
  }

  return(
    list(
      na_percentage_col = na_percentage_col,
      na_percentage_row = na_percentage_row
    )
  )
}
