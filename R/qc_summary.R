utils::globalVariables(c("DAid"))
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


qc_summary <- function() {

}
