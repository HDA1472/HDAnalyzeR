test_that("Directory creation without date", {
  dir_name <- "test_directory_without_date"
  create_dir(dir_name)
  expect_true(dir.exists(dir_name), "Directory should be created")
  unlink(dir_name, recursive = TRUE)
})

test_that("Directory creation with date", {
  dir_prefix <- "test_directory_with_date"
  create_dir(dir_prefix, date = TRUE)

  # List all directories that match the prefix
  created_dirs <- list.files(pattern = paste0("^", dir_prefix, "_\\d{4}_\\d{2}_\\d{2}_\\d{6}$"))

  # Check that at least one directory with the expected prefix and date format exists
  expect_true(length(created_dirs) > 0, "Directory with date should be created")
  unlink(created_dirs, recursive = TRUE)
})

test_that("Handling existing directory", {
  dir_name <- "existing_directory"
  create_dir(dir_name)
  create_dir(dir_name)  # Attempt to create again
  expect_true(dir.exists(dir_name), "Existing directory should still exist")
  unlink(dir_name, recursive = TRUE)
})

test_that("Warning on failed directory creation", {
  dir_name <- "C:/Windows/System32/restricted_directory"  # A directory that likely requires superuser privileges
  expect_warning(create_dir(dir_name), "Failed to create directory")
  unlink(dir_name, recursive = TRUE)
})
