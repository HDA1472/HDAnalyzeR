#' Create directory w/o system date
#'
#' The function creates a directory with the specified name. If the directory already exists, a message is printed.
#'
#' @param dir_name (string). The name of the directory to create
#' @param date (logical). If TRUE, the current date and time will be appended to the directory name
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
#' @param df (tibble). The dataframe to save
#' @param dir_name (string). The directory where the file will be saved
#' @param file_name (string). The name of the file to save
#' @param file_type (string). The type of file to save the dataframe as. Options are "csv", "tsv", or "rda"
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
#' @param file_path (string). The path to the file to import
#'
#' @return df (tibble). The imported dataframe
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

  df <- switch(tolower(file_extension),
               csv = utils::read.csv(file_path, stringsAsFactors = FALSE),
               tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path),
               stop("Unsupported file type: ", file_extension))

  return(df)
}


#' Check for NAs in a column and remove these rows
#'
#' The function checks for NAs in the specified column and removes the rows with NAs.
#' It returns a warning message if any rows are removed stating the number and the indexes of removed rows.
#'
#' @param df_in (tibble). The input dataframe
#' @param cols (string or vector of strings). The column to check for NAs
#'
#' @return df_out (tibble). The dataframe with NAs removed
#' @export
#'
#' @examples
#' df <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 2, 3, 4))
#' df_out <- remove_na(df, "x")
remove_na <- function(df_in, cols) {

  rows_to_omit <- integer(0)  # Keeps track of rows to omit
  warning_messages <- character(0)  # Keeps track of the warning messages

  for (col in cols) {

    if (!is.null(col) && col %in% colnames(df_in)) {
      na_rows <- which(is.na(df_in[[col]]))
      num_na_rows <- length(na_rows)
      rows_to_omit <- c(rows_to_omit, na_rows)

      if (num_na_rows > 0) {
        warning_message <- sprintf("Omitted %d rows with NAs in column '%s'. Indices of omitted rows: \n%s",
                                   num_na_rows, col, paste(na_rows, collapse = ", "))
        warning_messages <- c(warning_messages, warning_message)
      }

    }

  }

  # Omit the rows with NAs from the dataframe and show warning messages
  rows_to_omit <- unique(rows_to_omit)

  if (length(warning_messages) > 0) {
    df_out <- df_in[-rows_to_omit, ]
    warning(paste(warning_messages, collapse = "\n"))
  } else {
    df_out <- df_in
  }

  rownames(df_out) <- NULL  # Re index the rows

  return(df_out)
}


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
#' na_percentages <- calc_na_percentage(example_metadata)
#' print(na_percentages)
calc_na_percentage <- function(df) {

  na_percentage <- df |>
    dplyr::summarise_all(~ round(sum(is.na(.) / dplyr::n() * 100), 1)) |>
    tidyr::gather(key = "column", value = "na_percentage") |>
    dplyr::filter(na_percentage > 0) |>  # Filter out columns with no NAs
    dplyr::arrange(dplyr::desc(na_percentage))

  return(na_percentage)
}
