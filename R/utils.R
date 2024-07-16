utils::globalVariables(c("DAid", "Assay", "NPX"))
#' Create directory w/o system date
#'
#' The function creates a directory with the specified name. If the directory already exists, a message is printed.
#'
#' @param dir_name (string). The name of the directory to create
#' @param date (logical). If T, a directory with the current date as name will be created in the directory with `dir_name`.
#'
#' @return dir_name (string). The name of the created directory
#' @export
#'
#' @examples
#' create_dir("my_directory", date = FALSE)  # create outer directory
#' create_dir("my_directory", date = TRUE)  # create inner directory with date as name
#' # Clean up the created directory
#' unlink("my_directory", recursive = TRUE)
create_dir <- function(dir_name, date = F) {

  if (date) {
    current_date <- format(Sys.time(), "%Y_%m_%d")  # Get the current date
    dir_name <- paste0(dir_name, "/", current_date)
  }

  # Check if the directory already exists
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)

    # Check if the directory was created successfully
    if (dir.exists(dir_name)) {
    } else {
      warning(paste("Failed to create directory", dir_name))
    }

  } else {
    message(paste("Directory", dir_name, "already exists."))
  }

  return(dir_name)
}


#' Save dataframe
#'
#' The function saves a dataframe in the specified format (CSV, TSV, or RDA) in the specified directory.
#'
#' @param df (tibble). The dataframe to save
#' @param file_name (string). The name of the file to save
#' @param dir_name (string). The directory where the file will be saved
#' @param date (logical). If T, a directory with the current date as name will be created in the directory with `dir_name`.
#' @param file_type (string). The type of file to save the dataframe as. Options are "csv", "tsv", or "rda"
#'
#' @return NULL
#' @export
#'
#' @examples
#' df <- data.frame(x = 1:10, y = rnorm(10))
#' save_df(df, "sample_data", "my_data", file_type = "csv")
#' # Clean up the created directory
#' unlink("my_data", recursive = TRUE)
save_df <- function(df, file_name, dir_name, date = F, file_type = c("csv", "tsv", "rda")) {

  valid_file_types <- c("csv", "tsv", "rda")

  if (!file_type %in% valid_file_types) {
    stop("Unsupported file type: ", file_type)
  }

  # Create the directory if it doesn't exist, else store the file in the existing directory
  dir_name <- create_dir(dir_name, date = date)

  file_path <- file.path(dir_name, paste0(file_name, ".", file_type))

  if (file_type == "csv") {
    utils::write.csv(df, file_path, row.names = FALSE)
  } else if (file_type == "tsv") {
    utils::write.table(df, file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else if (file_type == "rda") {
    save(df, file = file_path)
  }

  invisible(NULL)
}


#' Import dataframe
#'
#' The function imports a dataframe from a file in CSV, TSV, RDA, RDS, XLSX, or TXT format.
#'
#' @param file_path (string). The path to the file to import
#'
#' @return df (tibble). The imported dataframe
#' @export
#'
#' @examples
#' df_out <- example_data
#' save_df(df_out, "sample_data", "my_data", file_type = "rda")
#' df_in <- load("my_data/sample_data.rda")
#' # Clean up the created directory
#' unlink("my_data", recursive = TRUE)
import_df <- function(file_path) {

  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)

  df <- switch(tolower(file_extension),
               csv = utils::read.csv(file_path, stringsAsFactors = FALSE),
               tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path, guess_max=10000000),
               stop("Unsupported file type: ", file_extension))

  df <- tibble::as_tibble(df)
  return(df)
}


#' Widen data
#'
#' The function widens the data from long to wide format.
#'
#' @param olink_data (tibble). A dataframe containing Olink data to be normalized
#' @param wide (logical). A logical value indicating whether the data is in wide format
#'
#' @return wide_data (tibble). A dataframe containing the data in wide format
#' @export
#'
#' @examples
#' wide_data <- widen_data(example_data, wide = FALSE)
widen_data <- function(olink_data, wide) {

  if (isFALSE(wide)) {
    wide_data <- olink_data |>
      dplyr::select(DAid, Assay, NPX) |>
      tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  } else {
    wide_data <- olink_data
  }

  return(wide_data)
}
