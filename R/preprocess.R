utils::globalVariables(c("PlateID", "Cohort", "Assay_Warning", "QC_Warning"))
#' Clean data
#'
#' The function cleans the data by filtering out rows based on the specified criteria.
#' It keeps only the specified columns.
#'
#' @param df_in The input dataframe
#' @param keep_cols The columns to keep in the output dataframe
#' @param cohort The cohort to keep
#' @param exclude_plates The plates to exclude
#' @param filter_assay If TRUE only rows with Assay_Warning == "PASS" are kept, else NULL
#' @param filter_qc If TRUE only rows with QC_Warning == "PASS" are kept, else NULL
#'
#' @return df_out The cleaned dataframe
#' @export
#'
#' @examples
#' df <- clean_data(example_data, exclude_plates = c("Plate1", "Plate2"), filter_assay = TRUE)
clean_data <- function(df_in, keep_cols = c("DAid", "Assay", "NPX"), cohort = NULL,
                       exclude_plates = NULL, filter_assay = NULL, filter_qc= NULL) {

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
                    is.null(filter_assay) | Assay_Warning == "PASS"
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if ("QC_Warning" %in% colnames(df_in)) {
                    is.null(filter_qc) | QC_Warning == "PASS"
                  } else {
                    TRUE
                  }) |>
    dplyr::select(dplyr::any_of(keep_cols))

    return(df_out)
}


clean_metadata <- function(df_in, keep_cols = c("DAid", "Disease", "Sex", "Age", "BMI"), exclude_sample = NULL) {

  df_out <- df_in |>
    dplyr::filter(if ("Exclude Sample" %in% colnames(df_in)) {
      is.null(exclude_sample) | `Exclude Sample` == "PASS"
    } else {
      TRUE
    }) |>
    dplyr::select(dplyr::any_of(keep_cols))

  return(df_out)
}
