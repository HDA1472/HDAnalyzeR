#' Impute via Median
#'
#' `impute_median()` imputes missing values in a dataset using the median of each column.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation.
#'
#' @param olink_data The input dataset.
#' @param wide If TRUE, the data is in wide format.
#' @param exclude_cols The columns to exclude from imputation.
#' @param show_na_percentage If TRUE, the percentage of missing values in each column is displayed.
#'
#' @return The imputed dataset.
#' @export
#'
#' @details This is the fastest but usually least accurate imputation method.
#' @examples
#' # Data before imputation
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::slice_head(n = 100)
#' test_data
#'
#' # Data after imputation
#' impute_median(test_data)
impute_median <- function(olink_data,
                          wide = TRUE,
                          exclude_cols = c("DAid", "Disease"),
                          show_na_percentage = TRUE) {

  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
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
#' `impute_knn()` imputes missing values in a dataset using the k-nearest neighbors method.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation. The user can also specify
#' the number of neighbors to consider for imputation.
#'
#' @param olink_data The input dataset.
#' @param wide If TRUE, the data is in wide format.
#' @param k The number of neighbors to consider for imputation.
#' @param exclude_cols The columns to exclude from imputation.
#' @param show_na_percentage If TRUE, the percentage of missing values in each column is displayed.
#'
#' @return The imputed dataset.
#' @export
#'
#' @examples
#' # Data before imputation
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::slice_head(n = 100)
#' test_data
#'
#' # Data after imputation
#' impute_knn(test_data, k = 3)
impute_knn <- function(olink_data,
                       wide = TRUE,
                       k = 5,
                       exclude_cols = c("DAid", "Disease"),
                       show_na_percentage = TRUE) {

  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
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
#' `impute_missForest()` imputes missing values in a dataset using the `missForest` method.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation. The user can also specify
#' the maximum number of iterations, the number of trees to grow, the type of parallelization
#' ("no", "variables", or "forests"), as well as the number of cores to use for parallelization.
#'
#' @param olink_data The input dataset.
#' @param wide If TRUE, the data is in wide format.
#' @param maxiter The maximum number of iterations.
#' @param ntree  The number of trees to grow.
#' @param parallelize The type of parallelization to use. Options are "no", "variables", or "forests".
#' @param ncores The number of cores to use for parallelization.
#' @param exclude_cols The columns to exclude from imputation.
#' @param show_na_percentage If TRUE, the percentage of missing values in each column is displayed.
#'
#' @return The imputed dataset.
#' @export
#'
#' @examples
#' # Data before imputation
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::slice_head(n = 100)
#' test_data
#'
#' # Data after imputation
#' impute_missForest(test_data, maxiter = 1, ntree = 50, parallelize = "no")
impute_missForest <- function(olink_data,
                              wide = TRUE,
                              maxiter = 10,
                              ntree = 100,
                              parallelize = "variables",
                              ncores = 4,
                              exclude_cols = c("DAid", "Disease"),
                              show_na_percentage = TRUE) {

  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
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


#' Impute via MICE
#'
#' `impute_mice()` imputes missing values in a dataset using the MICE algorithm.
#' It allows the user to exclude certain columns from imputation and can also display the
#' percentage of missing values in each column before imputation. The user can also specify
#' the number of imputed datasets to create, the maximum number of iterations, and the imputation
#' method to use.
#'
#' @param olink_data The input dataset.
#' @param wide If TRUE, the data is in wide format.
#' @param m The number of imputed datasets to create. Default is 5.
#' @param maxit The maximum number of iterations. Default is 5.
#' @param method The imputation method to use. Type `methods(mice::mice)` after `library(mice)` for all options. Default is "pmm".
#' @param exclude_cols The columns to exclude from imputation.
#' @param show_na_percentage If TRUE, the percentage of missing values in each column is displayed.
#'
#' @return The imputed dataset.
#' @export
#'
#' @examples
#' # Data before imputation
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = "Assay", values_from = "NPX") |>
#'   dplyr::slice_head(n = 100)
#' test_data
#'
#' # Data after imputation
#' impute_mice(test_data)
impute_mice <- function(olink_data,
                        wide = TRUE,
                        m = 5,
                        maxit = 5,
                        method = "pmm",
                        exclude_cols = c("DAid", "Disease"),
                        show_na_percentage = TRUE) {

  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
  data_in <- wide_data |>
    dplyr::select(-dplyr::any_of(exclude_cols))

  if (isTRUE(show_na_percentage)) {
    na_percentages <- calc_na_percentage_col(data_in)
    print(na_percentages)
  }

  data_in <- mice::mice(data_in, m = m, maxit = maxit, method = method, seed = 123, printFlag = F)
  imputed_data <- mice::complete(data_in)

  cols <- wide_data |>
    dplyr::select(dplyr::any_of(exclude_cols))
  imputed_data <- dplyr::bind_cols(cols, imputed_data)

  return(imputed_data)
}
