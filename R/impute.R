#' Impute via k-nearest neighbors
#'
#' This function imputes missing values in a dataset using the k-nearest neighbors algorithm.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation.
#'
#' @param wide_data (tibble). The input dataframe
#' @param k (integer). The number of neighbors to consider for imputation
#' @param exclude_cols (string or vector of strings). The columns to exclude from imputation
#' @param show_na_percentage (TRUE or NULL). If TRUE, the percentage of missing values in each column is displayed
#'
#' @return imputed_data (tibble). The dataframe with missing values imputed
#' @export
#'
#' @examples
#' test_data <- tibble::tibble(
#'   A = c(80, 44, NA, 50, 29),
#'   B = c(30, NA, 85, 70, 54),
#'   C = c(7, 10, 25, 74, 49)
#' )
#' imputed_data <- impute_knn(test_data, k = 3)
impute_knn <- function(wide_data, k = 5, exclude_cols = c("DAid", "Disease"),
                       show_na_percentage = TRUE) {

  data_in <- wide_data |>
    dplyr::select(-dplyr::any_of(exclude_cols))

  if (show_na_percentage) {
    na_percentages <- calc_na_percentage(data_in)
    print(na_percentages)
  }

  recipe <- recipes::recipe(~ ., data = data_in) |>
    recipes::step_impute_knn(recipes::all_predictors(), neighbors = k)

  imputed_data <- recipe |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  cols <- wide_data |>
    dplyr::select(dplyr::any_of(exclude_cols))
  imputed_data <- dplyr::bind_cols(cols, imputed_data)

  return(imputed_data)
}
