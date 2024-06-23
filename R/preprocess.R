utils::globalVariables(c("PlateID", "Cohort", "Assay_Warning", "QC_Warning", "Exclude Sample"))
#' Clean data
#'
#' The function cleans the data by filtering out rows based on the specified criteria.
#' It keeps only the specified columns. It removes rows with NAs in the DAid and NPX columns.
#'
#' @param df_in (tibble). The input dataframe
#' @param keep_cols (string or vector of strings). The columns to keep in the output dataframe
#' @param cohort (string or vector of strings). The cohort to keep
#' @param exclude_plates (string or vector of strings). The plates to exclude
#' @param filter_assay (TRUE or NULL). If TRUE only rows with Assay_Warning == "PASS" are kept, else NULL
#' @param filter_qc (TRUE or NULL). If TRUE only rows with QC_Warning == "PASS" are kept, else NULL
#'
#' @return df_out (tibble). The cleaned dataframe
#' @export
#'
#' @examples
#' df <- clean_data(example_data, exclude_plates = c("Plate1", "Plate2"), filter_assay = TRUE)
clean_data <- function(df_in, keep_cols = c("DAid", "Assay", "NPX"), cohort = NULL,
                       exclude_plates = NULL, filter_assay = NULL, filter_qc = NULL) {


  df_out <- df_in |>
    dplyr::filter(if ("PlateID" %in% colnames(df_in)) {
                    !(PlateID %in% exclude_plates)
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if ("Cohort" %in% colnames(df_in)) {
                    is.null(cohort) | Cohort %in% cohort
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if ("Assay_Warning" %in% colnames(df_in)) {
                    is.null(filter_assay) | Assay_Warning %in% filter_assay
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if ("QC_Warning" %in% colnames(df_in)) {

                    is.null(filter_qc) | QC_Warning %in% filter_qc
                  } else {
                    TRUE
                  }) |>
    dplyr::select(dplyr::any_of(keep_cols))

    df_out <- remove_na(df_out, c("DAid", "NPX"))

    return(df_out)
}


#' Clean metadata
#'
#' The function cleans the metadata by filtering out rows based on the specified criteria.
#' It keeps only the specified columns. It removes rows with NAs in the DAid and Disease columns.
#'
#' @param df_in (tibble). The input metadata
#' @param keep_cols (string or vector of strings). The columns to keep in the output metadata
#' @param exclude_sample (string). The value in the "Exclude Sample" column to exclude
#'
#' @return df_out (tibble). The cleaned metadata
#' @export
#'
#' @examples
#' df <- clean_metadata(example_metadata, exclude_sample = "yes")
clean_metadata <- function(df_in, keep_cols = c("DAid", "Disease", "Sex", "Age", "BMI"), exclude_sample = NULL) {

  df_out <- df_in |>
    dplyr::filter(if ("Exclude Sample" %in% colnames(df_in)) {
      is.null(exclude_sample) | !(`Exclude Sample` %in% exclude_sample)
    } else {
      TRUE
    }) |>
    dplyr::select(dplyr::any_of(keep_cols))

  df_out <- remove_na(df_out, c("DAid", "Disease"))

  return(df_out)
}


#' Generate and store wide and join dataframes
#'
#' The function generates wide and join dataframes from the long data and metadata.
#' It stores them in the data/processed/data_metadata directory in RDA format.
#'
#' @param long_data (tibble). The long data
#' @param metadata (tibble). The metadata
#'
#' @return list(wide_data, join_data) (list of tibbles). The wide and join dataframes
#' @export
#'
#' @examples
#' clean_data <- clean_data(example_data, keep_cols = c("DAid", "Assay", "NPX"))
#' clean_metadata <- clean_metadata(example_metadata, keep_cols = c("DAid", "GROUP", "Age"))
#' result_df <- generate_df(clean_data, clean_metadata)
#' wide_data <- result_df[[1]]
#' join_data <- result_df[[2]]
#' # Clean up the created directory
#' unlink("data", recursive = TRUE)
generate_df <- function(long_data, metadata) {

  wide_data <- tidyr::pivot_wider(long_data, names_from = "Assay", values_from = "NPX")

  join_data <- wide_data |>
    dplyr::left_join(metadata, by = "DAid")

  create_dir("data/processed/data_metadata", date = TRUE)
  save_df(long_data, "data/processed/data_metadata", "long_data", "rda")
  save_df(wide_data, "data/processed/data_metadata", "wide_data", "rda")
  save_df(metadata, "data/processed/data_metadata", "metadata", "rda")
  save_df(join_data, "data/processed/data_metadata", "join_data", "rda")

  return(list(wide_data, join_data))
}
