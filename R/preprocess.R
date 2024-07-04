utils::globalVariables(c("PlateID", "Cohort", "Assay_Warning", "Exclude Sample"))
#' Replace values with NA
#'
#' The function replaces the specified values with NA in the input vector.
#'
#' @param df_in (tibble). The input dataframe.
#' @param replace_w_na (vector). The values to replace with NA. Default is c(0, "0", "", "Unknown", "unknown", "none", NA, "na").
#'
#' @return df_out (tibble). The dataframe with specified values replaced with NA.
#' @keywords internal
replace_with_na <- function(df_in, replace_w_na = c(0, "0", "", "Unknown", "unknown", "none", NA, "na")) {
  df_out <- df_in |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(. %in% replace_w_na, NA, .)))

  return(df_out)
}


#' Clean data
#'
#' The function cleans the data by filtering out rows based on the specified criteria.
#' It keeps only the specified columns. It removes rows with NAs in the DAid and NPX columns.
#'
#' @param df_in (tibble). The input dataframe.
#' @param keep_cols (vector). The columns to keep in the output dataframe.
#' @param cohort (vector). The cohort to keep.
#' @param filter_plates (vector). The plates to exclude.
#' @param filter_assays (vector). The assays to filter out.
#' @param filter_assay_warning (logical). If TRUE, only the rows with Assay_Warning == "PASS" are kept.
#' @param remove_na_cols (vector). The columns to check for NAs and remove respective rows.
#' @param replace_w_na (vector). The values to replace with NA. Default is c(0, "0", "", "Unknown", "unknown", "none", NA, "na").
#'
#' @return df_out (tibble). The cleaned dataframe.
#' @export
#'
#' @examples
#' df <- clean_data(example_data, filter_plates = c("Plate1", "Plate2"), filter_assay_warning = TRUE)
clean_data <- function(df_in, keep_cols = c("DAid", "Assay", "NPX"), cohort = NULL, filter_plates = NULL,
                       filter_assays = NULL, filter_assay_warning = F, remove_na_cols = c("DAid", "NPX"),
                       replace_w_na = c(0, "0", "", "Unknown", "unknown", "none", NA, "na")) {

  df_out <- df_in |>
    dplyr::filter(if ("PlateID" %in% colnames(df_in)) {
                    !(PlateID %in% filter_plates)
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if ("Cohort" %in% colnames(df_in)) {
                    isFALSE(cohort) | Cohort %in% cohort
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(if (filter_assay_warning == T && "Assay_Warning" %in% colnames(df_in)) {
                    (Assay_Warning == "PASS")
                  } else {
                    TRUE
                  }) |>
    dplyr::filter(!(Assay %in% filter_assays)) |>
    dplyr::select(dplyr::any_of(keep_cols))

    if (!is.null(replace_w_na)) {
      df_out <- replace_with_na(df_out, replace_w_na)
    }

    if (!is.null(remove_na_cols)) {
      rows_before <- nrow(df_out)
      df_out <- stats::na.omit(df_out, target.colnames = remove_na_cols)
      rows_after <- nrow(df_out)
      if (rows_before != rows_after) {
        message("Removed ", rows_before - rows_after, " rows with NAs based on ", remove_na_cols)
      }
    }

    return(df_out)
}


#' Clean metadata
#'
#' The function cleans the metadata by filtering out rows based on the specified criteria.
#' It keeps only the specified columns. It removes rows with NAs in the DAid and Disease columns.
#'
#' @param df_in (tibble). The input metadata.
#' @param keep_cols (vector). The columns to keep in the output metadata.
#' @param filter_samples (vector). The samples to filter out.
#' @param remove_na_cols (vector). The columns to check for NAs and remove respective rows.
#' @param replace_w_na (vector). The values to replace with NA. Default is c("Unknown", "unknown", "none", NA, "na").
#'
#' @return df_out (tibble). The cleaned metadata.
#' @export
#'
#' @examples
#' df <- clean_metadata(example_metadata)
clean_metadata <- function(df_in, keep_cols = c("DAid", "Disease", "Sex", "Age", "BMI"),
                           filter_samples = NULL, remove_na_cols = c("DAid", "Disease"),
                           replace_w_na = c("Unknown", "unknown", "none", NA, "na")) {

  df_out <- df_in |>
    dplyr::filter(!(DAid %in% filter_samples)) |>
    dplyr::select(dplyr::any_of(keep_cols))

  if (!is.null(replace_w_na)) {
    df_out <- replace_with_na(df_out, replace_w_na)
  }

  if (!is.null(remove_na_cols)) {
    rows_before <- nrow(df_out)
    df_out <- stats::na.omit(df_out, target.colnames = remove_na_cols)
    rows_after <- nrow(df_out)
    if (rows_before != rows_after) {
      message("Removed ", rows_before - rows_after, " rows with NAs based on ", remove_na_cols)
    }
  }

  return(df_out)
}


#' Generate and store wide and join dataframes
#'
#' The function generates wide and join dataframes from the long data and metadata.
#' It stores them in the data/processed/data_metadata directory in RDA format.
#'
#' @param long_data (tibble). The long data.
#' @param metadata (tibble). The metadata.
#' @param join (logical). If TRUE, the dataframes are joined with metadata.
#' @param metadata_cols (vector). The metadata columns to join with data.
#' @param save (logical). If TRUE, the dataframes are saved.
#'
#' @return A list containing the following elements:
#'  - wide_data (tibble). The wide data.
#'  - join_data (tibble). The joined data with metadata. If join is FALSE, this is not returned.
#' @export
#'
#' @examples
#' clean_data <- clean_data(example_data, keep_cols = c("DAid", "Assay", "NPX"))
#' clean_metadata <- clean_metadata(example_metadata, keep_cols = c("DAid", "GROUP", "Age"))
#' result_df <- generate_df(clean_data, clean_metadata)
#' wide_data <- result_df$wide_data
#' join_data <- result_df$join_data
#' # Clean up the created directory
#' unlink("data", recursive = TRUE)
generate_df <- function(long_data, metadata = NULL, join = T, metadata_cols = c("DAid", "Disease", "Sex", "Age", "BMI"), save = T) {

  wide_data <- widen_data(long_data, wide = F)

  if (isTRUE(join)) {
    join_data <- wide_data |>
      dplyr::left_join(metadata |> dplyr::select(dplyr::any_of(metadata_cols)), by = "DAid")
    if (isTRUE(save)) {
      dir_name <- create_dir("data/processed")
      save_df(wide_data, "wide_data", dir_name, date = T, "rda")
      save_df(join_data, "join_data", dir_name, date = T, "rda")
    }

    return(list("wide_data" = wide_data, "join_data" = join_data))

  } else {
    join_data <- NULL
    if (isTRUE(save)) {
      dir_name <- create_dir("data/processed")
      save_df(wide_data, "wide_data", dir_name, date = T, "rda")
    }
    return(wide_data)
  }
}
