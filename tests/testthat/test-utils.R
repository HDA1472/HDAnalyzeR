# Test create_dir --------------------------------------------------------------
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

# Test save_df -----------------------------------------------------------------
test_that("save_df creates directory and saves CSV", {
  dir_name <- "test_csv_dir"
  file_name <- "test_csv_file"

  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
  }

  save_df(df, dir_name, file_name, "csv")
  expect_true(dir.exists(dir_name))
  expect_true(file.exists(file.path(dir_name, paste0(file_name, ".csv"))))
  unlink(dir_name, recursive = TRUE)
})


test_that("save_df saves TSV", {
  dir_name <- "test_tsv_dir"
  file_name <- "test_tsv_file"

  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
  }

  save_df(df, dir_name, file_name, "tsv")
  expect_true(dir.exists(dir_name))
  expect_true(file.exists(file.path(dir_name, paste0(file_name, ".tsv"))))
  unlink(dir_name, recursive = TRUE)
})

test_that("save_df saves RDA", {
  dir_name <- "test_rda_dir"
  file_name <- "test_rda_file"

  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
  }

  save_df(df, dir_name, file_name, "rda")
  expect_true(dir.exists(dir_name))
  expect_true(file.exists(file.path(dir_name, paste0(file_name, ".rda"))))
  unlink(dir_name, recursive = TRUE)
})

test_that("save_dataframe works with existing directory", {
  dir_name <- "existing_dir"
  file_name <- "test_existing_file"

  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }

  save_df(df, dir_name, file_name, "csv")
  expect_true(file.exists(file.path(dir_name, paste0(file_name, ".csv"))))
  unlink(dir_name, recursive = TRUE)
})

test_that("save_dataframe handles invalid file type", {
  dir_name <- "test_invalid_dir"
  file_name <- "test_invalid_file"
  invalid_file_type <- "invalid"

  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
  }

  expected_error_message <- paste("Unsupported file type:", invalid_file_type)
  expect_error(save_df(df, dir_name, file_name, invalid_file_type), expected_error_message)
  expect_false(dir.exists(dir_name))
})
