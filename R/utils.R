#' Create directory w/o system date
#'
#' The function creates a directory with the specified name. If the directory already exists, a message is printed.
#'
#' @param dir_name The name of the directory to create
#' @param date Logical. If TRUE, the current date and time will be appended to the directory name
#'
#' @return NULL
#'
#' @examples
#' create_directory("my_directory", date = FALSE)
#' create_directory("my_directory/inner_dir", date = TRUE)
create_directory <- function(dir_name, date = FALSE) {

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
