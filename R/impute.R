#' Impute via Median
#'
#' This function imputes missing values in a dataset using the median of each column.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation.
#'
#' @param wide_data (tibble). The input dataframe
#' @param exclude_cols (string or vector of strings). The columns to exclude from imputation
#' @param show_na_percentage (logical). If T, the percentage of missing values in each column is displayed
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
#' imputed_data <- impute_median(test_data)
impute_median <- function(wide_data, exclude_cols = c("DAid", "Disease"),
                          show_na_percentage = T) {

  data_in <- wide_data |>
    dplyr::select(-dplyr::any_of(exclude_cols))

  if (isTRUE(show_na_percentage)) {
    na_percentages <- calc_na_percentage_col(data_in)
    print(na_percentages)
  }

  recipe <- recipes::recipe(~ ., data = data_in) |>
    recipes::step_impute_median(recipes::all_predictors())

  imputed_data <- recipe |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

  cols <- wide_data |>
    dplyr::select(dplyr::any_of(exclude_cols))
  imputed_data <- dplyr::bind_cols(cols, imputed_data)

  return(imputed_data)
}


#' Impute via k-nearest neighbors
#'
#' This function imputes missing values in a dataset using the k-nearest neighbors algorithm.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation.
#'
#' @param wide_data (tibble). The input dataframe
#' @param k (integer). The number of neighbors to consider for imputation
#' @param exclude_cols (string or vector of strings). The columns to exclude from imputation
#' @param show_na_percentage (logical). If T, the percentage of missing values in each column is displayed
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

  if (isTRUE(show_na_percentage)) {
    na_percentages <- calc_na_percentage_col(data_in)
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


#' Impute via missForest
#'
#' This function imputes missing values in a dataset using the `missForest` method.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation.
#'
#' @param wide_data (tibble). The input dataframe
#' @param maxiter (integer). The maximum number of iterations
#' @param ntree (integer). The number of trees to grow
#' @param parallelize (string). The type of parallelization to use. Options are "no", "variables", or "forests"
#' @param ncores (integer). The number of cores to use for parallelization
#' @param exclude_cols (string or vector of strings). The columns to exclude from imputation
#' @param show_na_percentage (logical). If T, the percentage of missing values in each column is displayed
#'
#' @return A data frame with imputed values.
#' @export
#'
#' @examples
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::slice_head(n = 100)
#' imputed_data <- impute_missForest(test_data, maxiter = 1, ntree = 50, parallelize = "no")
impute_missForest <- function(wide_data, maxiter = 10, ntree = 100, parallelize = "variables",
                              ncores = 4, exclude_cols = c("DAid", "Disease"),
                              show_na_percentage = T) {

  data_in <- wide_data |>
    dplyr::select(-dplyr::any_of(exclude_cols))

  if (isTRUE(show_na_percentage)) {
    na_percentages <- calc_na_percentage_col(data_in)
    print(na_percentages)
  }
  if (parallelize == "no") {
    ncores <- 1
  } else {
    doParallel::registerDoParallel(cores = ncores)
  }
  set.seed(123)
  data_in <- as.data.frame(data_in)  # Convert to data frame for missForest
  imputed_data <- missForest::missForest(data_in, maxiter = maxiter, ntree = ntree,
                                         verbose = T, parallelize = parallelize)$ximp
  imputed_data <- tibble::as_tibble(imputed_data)

  cols <- wide_data |>
    dplyr::select(dplyr::any_of(exclude_cols))
  imputed_data <- dplyr::bind_cols(cols, imputed_data)

  return(imputed_data)
}
