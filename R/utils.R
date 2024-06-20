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


#' Save dataframe in CSV, TSV, or RDA format
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
  file_type <- match.arg(file_type)

  # Create the directory if it doesn't exist, else store the file in the existing directory
  create_dir(dir_name)

  file_path <- file.path(dir_name, paste0(file_name, ".", file_type))

  if (file_type == "csv") {
    utils::write.csv(df, file_path, row.names = FALSE)
  } else if (file_type == "tsv") {
    utils::write.table(df, file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else if (file_type == "rda") {
    save(df, file = file_path)
  } else {
    stop("Unsupported file type: ", file_type)
  }

  invisible(NULL)

}
