utils::globalVariables(c("DAid", "Assay", "NPX"))
#' Create directory
#'
#' `create_dir()` creates a directory with a specified name.
#' The user can choose to create another inner directory with the current date as its name.
#' If the directory already exists, a message is printed.
#'
#' @param dir_name The name of the directory to create.
#' @param date If TRUE, a directory with the current date as name will be created inside the directory with `dir_name`.
#'
#' @return The relative file path of the created directory as a string.
#' @export
#'
#' @examples
#' # Create a directory with a specified name
#' create_dir("my_directory", date = FALSE)
#' unlink("my_directory", recursive = TRUE)  # Clean up the created directory
#'
#' # Create a directory with a specified name and an inner directory with the current date as name
#' create_dir("my_directory", date = TRUE)
#' unlink("my_directory", recursive = TRUE)  # Clean up the created directory
#'
#' # Create a directory inside another directory
#' create_dir("outer_directory/inner_directory", date = FALSE)
#' unlink("outer_directory", recursive = TRUE)  # Clean up the created directory
#'
#' # Create a directory inside a pre existing one
#' create_dir("outer_directory", date = FALSE)
#' create_dir("outer_directory/inner_directory", date = FALSE)
#' unlink("outer_directory", recursive = TRUE)  # Clean up the created directory
create_dir <- function(dir_name, date = FALSE) {

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


#' Save tibble as CSV, TSV, Excel or RDA file
#'
#' `save_df()` saves a dataframe in the specified format (csv, tsv, rda, or xlsx) in a
#' specified directory. If the directory does not exist, it will be created.
#' The recommended file type is RDA.
#'
#' @param df The dataframe to save.
#' @param file_name The name of the file to save.
#' @param dir_name The directory where the file will be saved.
#' @param date If TRUE, a directory with the current date as name will be created in the directory with `dir_name`.
#' @param file_type The type of file to save the dataframe as. Options are "csv", "tsv", "rda", or "xlsx".
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Save a metadata dataframe as an RDA file
#' save_df(example_metadata, "metadata", "my_data", file_type = "rda")
#'
#' file.exists("my_data/metadata.rda")  # Check if the file exists
#' unlink("my_data", recursive = TRUE)  # Clean up the created directory
save_df <- function(df, file_name, dir_name, date = FALSE, file_type = c("csv", "tsv", "rda", "xlsx")) {

  valid_file_types <- c("csv", "tsv", "rda", "xlsx")

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
  } else if (file_type == "xlsx") {
    writexl::write_xlsx(df, file_path)
  }

  invisible(NULL)
}


#' Import dataframe from file
#'
#' `import_df()` imports a dataframe from a file. The file format can be CSV,
#' TSV, TXT, RDA, RDS, XLSX, or Parquet format. It recognizes the file format,
#' reads it and returns it as a tibble.
#'
#' @param file_path The path to the file to import.
#'
#' @return The imported dataframe as a tibble.
#' @export
#'
#' @examples
#' # Save a dataframe as an RDA file
#' save_df(example_metadata, "metadata", "my_data", file_type = "rda")
#'
#' # Import the saved RDA file again as a tibble
#' import_df("my_data/metadata.rda")
#'
#' unlink("my_data", recursive = TRUE)  # Clean up the created directory
import_df <- function(file_path) {

  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)

  df <- switch(tolower(file_extension),
               csv = readr::read_csv(file_path),
               tsv = readr::read_tsv(file_path),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path, guess_max=10000000),
               parquet = arrow::read_parquet(file_path),
               stop("Unsupported file type: ", file_extension))

  df <- tibble::as_tibble(df)
  return(df)
}


#' Widen Olink data
#'
#' `widen_data()` transforms the data from long to wide format. It should be used to
#' transform Olink data from long to wide format with Assays as columns names and NPX as values.
#' The first column contains the DAids.
#'
#' @param olink_data A tibble containing Olink data to be transformed.
#'
#' @return A tibble containing the data in wide format.
#' @export
#'
#' @examples
#' # Olink data in long format
#' example_data
#'
#' # Transform Olink data in wide format
#' widen_data(example_data)
widen_data <- function(olink_data) {

  wide_data <- olink_data |>
    dplyr::select(DAid, Assay, NPX) |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)

  return(wide_data)
}
