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
  df <- data.frame(x = 1:10, y = rnorm(10))
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
  df <- data.frame(x = 1:10, y = rnorm(10))
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
  df <- data.frame(x = 1:10, y = rnorm(10))
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

test_that("save_df works with existing directory", {
  df <- data.frame(x = 1:10, y = rnorm(10))
  dir_name <- "existing_dir"
  file_name <- "test_existing_file"

  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }

  save_df(df, dir_name, file_name, "csv")
  expect_true(file.exists(file.path(dir_name, paste0(file_name, ".csv"))))
  unlink(dir_name, recursive = TRUE)
})

test_that("save_df handles invalid file type", {
  df <- data.frame(x = 1:10, y = rnorm(10))
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

# Test import_df ---------------------------------------------------------------
test_that("import_df handles CSV files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.csv"
  utils::write.csv(df_out, file_name, row.names = FALSE)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})


test_that("import_df handles TSV files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.tsv"
  utils::write.table(df_out, file_name, sep='\t', row.names = FALSE, col.names = TRUE)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})


test_that("import_df handles TXT files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.txt"
  utils::write.table(df_out, file_name, row.names = FALSE, col.names = TRUE)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})


test_that("import_df handles RDS files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.rds"
  saveRDS(df_out, file = file_name)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})


test_that("import_df handles RDA files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.rda"
  save(df_out, file = file_name)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})


test_that("import_df handles XLSX files", {
  df_out <- data.frame(x = 1:10, y = rnorm(10))
  file_name <- "test_file.xlsx"
  writexl::write_xlsx(df_out, file_name)
  expect_true(file.exists(file_name))

  if (file.exists(file_name)) {
    df_in <- import_df(file_name)
    expect_true(is.data.frame(df_in))
  } else {
    skip("Test skipped: File not found.")
  }

  unlink(file_name, recursive = TRUE)
})
