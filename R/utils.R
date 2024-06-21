#' Create directory w/o system date
#'
#' The function creates a directory with the specified name. If the directory already exists, a message is printed.
#'
#' @param dir_name The name of the directory to create
#' @param date Logical. If TRUE, the current date and time will be appended to the directory name
#'
#' @return NULL
#' @export
#'
#' @examples
#' create_dir("my_directory", date = FALSE)
#' create_dir("my_directory/inner_dir", date = TRUE)
#' # Clean up the created directory
#' unlink("my_directory", recursive = TRUE)
create_dir <- function(dir_name, date = FALSE) {

  if (date) {
    current_date <- format(Sys.time(), "%Y_%m_%d_%H%M%S")  # Get the current date and time
    dir_name <- paste0(dir_name, "_", current_date)
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

  invisible(NULL)

}


#' Save dataframe
#'
#' The function saves a dataframe in the specified format (CSV, TSV, or RDA) in the specified directory.
#'
#' @param df The dataframe to save
#' @param dir_name The directory where the file will be saved
#' @param file_name The name of the file to save
#' @param file_type The type of file to save the dataframe as. Options are "csv", "tsv", or "rda"
#'
#' @return NULL
#' @export
#'
#' @examples
#' df <- data.frame(x = 1:10, y = rnorm(10))
#' save_df(df, "my_data", "sample_data", "csv")
#' # Clean up the created directory
#' unlink("my_data", recursive = TRUE)
save_df <- function(df, dir_name, file_name, file_type = c("csv", "tsv", "rda")) {
  valid_file_types <- c("csv", "tsv", "rda")

  if (!file_type %in% valid_file_types) {
    stop("Unsupported file type: ", file_type)
  }

  # Create the directory if it doesn't exist, else store the file in the existing directory
  create_dir(dir_name)

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
#' @param file_path The path to the file to import
#'
#' @return df The imported dataframe
#' @export
#'
#' @examples
#' df_out <- example_data
#' save_df(df_out, "my_data", "sample_data", "rda")
#' df_in <- load("my_data/sample_data.rda")
#' # Clean up the created directory
#' unlink("my_data", recursive = TRUE)
import_df <- function(file_path) {
  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)

  data <- switch(tolower(file_extension),
                 csv = utils::read.csv(file_path, stringsAsFactors = FALSE),
                 tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
                 txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
                 rda = { load(file_path); get(ls()[1]) },
                 rds = readRDS(file_path),
                 xlsx = readxl::read_excel(file_path),
                 stop("Unsupported file type: ", file_extension))

  return(data)
}
