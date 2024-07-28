utils::globalVariables(c("PlateID", "Cohort", "Assay_Warning", "Exclude Sample"))
#' Replace specific values with NA
#'
#' `replace_with_na()` replaces the specified values in the input vector with NAs.
#'
#' @param df_in The input dataframe.
#' @param replace_w_na The values to replace with NA. Default is c(0, "0", "", "Unknown", "unknown", "none", NA, "na").
#'
#' @return The tibble with the specified values replaced with NA.
#' @keywords internal
replace_with_na <- function(df_in, replace_w_na = c(0, "0", "", "Unknown", "unknown", "none", NA, "na")) {
  df_out <- df_in |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(. %in% replace_w_na, NA, .)))

  return(df_out)
}


#' Preprocess data
#'
#' `clean_data()` preprocesses the data by filtering out rows based on the specified criteria.
#'   - It keeps only the data of the specified cohort, plates, and assays.
#'   - It can remove all samples with Assay_Warning != "PASS".
#'   - It keeps only the specified columns.
#'   - It removes rows with NAs in the DAid and NPX columns.
#'   - It replaces the specified values with NA.
#'
#' @param df_in The input dataframe.
#' @param keep_cols The columns to keep in the output dataframe.
#' @param cohort The cohort to keep.
#' @param filter_plates The plates to exclude.
#' @param filter_assays The assays to filter out.
#' @param filter_assay_warning If TRUE, only the rows with Assay_Warning == "PASS" are kept. Default is FALSE.
#' @param remove_na_cols The columns to check for NAs and remove respective rows. Defaults is c("DAid", "NPX").
#' @param replace_w_na The values to replace with NA. Default is c(0, "0", "", "Unknown", "unknown", "none", NA, "na").
#'
#' @return The preprocessed dataframe.
#' @export
#'
#' @examples
#' # Unprocessed data
#' example_data
#'
#' # Preprocessed data
#' clean_data(example_data, filter_plates = c("Plate1", "Plate2"), filter_assay_warning = TRUE)
clean_data <- function(df_in,
                       keep_cols = c("DAid", "Assay", "NPX"),
                       cohort = NULL,
                       filter_plates = NULL,
                       filter_assays = NULL,
                       filter_assay_warning = FALSE,
                       remove_na_cols = c("DAid", "NPX"),
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


#' Preprocess metadata
#'
#' `clean_metadata()` preprocesses the metadata by filtering out rows based on the specified criteria.
#'   - It keeps only the specified columns.
#'   - It removes rows with NAs in the DAid and Disease columns.
#'   - It replaces the specified values with NA.
#'
#' @param df_in The input metadata.
#' @param keep_cols The columns to keep in the output metadata.
#' @param remove_na_cols The columns to check for NAs and remove respective rows. Defaults is c("DAid", "Disease").
#' @param replace_w_na The values to replace with NA. Default is c("Unknown", "unknown", "none", NA, "na").
#'
#' @return The preprocessed metadata.
#' @export
#'
#' @examples
#' # Unprocessed metadata
#' example_metadata
#'
#' # Preprocessed metadata
#' clean_metadata(example_metadata)
clean_metadata <- function(df_in,
                           keep_cols = c("DAid", "Disease", "Sex", "Age", "BMI"),
                           remove_na_cols = c("DAid", "Disease"),
                           replace_w_na = c("Unknown", "unknown", "none", NA, "na")) {

  df_out <- df_in |>
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


#' Create and save wide and join dataframes
#'
#' `generate_df()` creates wide and join dataframes from the long data and metadata.
#' It saves them in the data/processed/data_metadata directory in RDA format.
#'
#' @param long_data The long data.
#' @param metadata The metadata.
#' @param join If TRUE, the dataframes are joined with metadata.
#' @param metadata_cols The metadata columns to join with data.
#' @param save If TRUE, the dataframes are saved.
#'
#' @return A list containing the following elements:
#'  - wide_data: The wide data.
#'  - join_data: The joined data with metadata. If join is FALSE, this is not returned.
#' @export
#'
#' @examples
#' # Preprocess data and metadata
#' clean_data <- clean_data(example_data, keep_cols = c("DAid", "Assay", "NPX"))
#' clean_metadata <- clean_metadata(example_metadata,
#'                                  keep_cols = c("DAid", "Disease", "Sex", "Age", "BMI"))
#'
#' # Create wide and join dataframes
#' result_df <- generate_df(clean_data, clean_metadata, save = FALSE)
#' result_df$wide_data
#' result_df$join_data
generate_df <- function(long_data,
                        metadata = NULL,
                        join = TRUE,
                        metadata_cols = c("DAid", "Disease", "Sex", "Age", "BMI"),
                        save = TRUE) {


  wide_data <- widen_data(long_data)

  if (isTRUE(join)) {
    join_data <- wide_data |>
      dplyr::left_join(metadata |> dplyr::select(dplyr::any_of(metadata_cols)), by = "DAid") |>
      dplyr::relocate(dplyr::all_of(metadata_cols), dplyr::everything())
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
