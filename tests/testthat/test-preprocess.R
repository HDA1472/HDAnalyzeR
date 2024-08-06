# Test replace_with_na ---------------------------------------------------------
test_that("The specified values are replaced with NAs", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data1 <- example_data
  test_data2 <- example_data
  test_data1$DAid[random_indices] <- 0
  test_data2$DAid[random_indices] <- NA
  result <- replace_with_na(test_data1)
  expected <- test_data2
  expect_equal(result, expected)
})

# Test clean_data --------------------------------------------------------------
test_that("clean_data selects specified columns", {
  result <- clean_data(example_data)
  expected <- example_data |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data excludes specified plates", {
  result <- clean_data(example_data, filter_plates = "Run001")
  expected <- example_data |>
    dplyr::filter(!(PlateID %in% "Run001")) |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data filters by cohort", {
  # Add Cohort column to filter based on
  result <- clean_data(example_data, cohort = "UCAN")
  expected <- example_data |>
    dplyr::filter(Cohort %in% "UCAN") |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data filters by Assay_Warning", {
  result <- clean_data(example_data, filter_assay_warning = T)
  expected <- example_data |>
    dplyr::filter(Assay_Warning == "PASS") |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data handles non-existent columns gracefully", {
  # Create data only with Assay and NPX columns
  partial_data <- example_data |>
    dplyr::select(Assay, NPX)
  result <- clean_data(partial_data, keep_cols = c("DAid", "Assay", "NPX"))
  expected <- partial_data |>
    dplyr::select(any_of(c("DAid", "Assay", "NPX")))
  expect_equal(result, expected)
})


test_that("clean_data handles all parameters together", {
  # Add Cohort column
  result <- clean_data(example_data, cohort = "UCAN", filter_plates = "P2",
                       filter_assay_warning = T)
  expected <- example_data |>
    dplyr::filter(Cohort %in% "UCAN" & !(PlateID %in% "P2") & Assay_Warning == "PASS") |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("The NA columns are removed", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data <- example_data
  test_data$DAid[random_indices] <- NA

  result <- clean_data(test_data, keep_cols = c("DAid", "NPX"))

  expected <- test_data |>
    dplyr::select(DAid, NPX) |>
    dplyr::filter(!is.na(DAid))

  result <- as.character(unclass(result))
  expected <- as.character(unclass(expected))

  expect_equal(result, expected)
})


test_that("The specified values are replaced with NAs", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data1 <- example_data
  test_data2 <- example_data
  test_data1$DAid[random_indices] <- 0
  test_data2$DAid[random_indices] <- NA
  result <- clean_data(test_data1, keep_cols = c("DAid", "NPX"), remove_na_cols = NULL)
  expected <- test_data2 |>
    dplyr::select(DAid, NPX)
  expect_equal(result, expected)
})


test_that("The specified values are not replaced with NAs", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data1 <- example_data
  test_data2 <- example_data
  test_data1$DAid[random_indices] <- 0
  result <- clean_data(test_data1, keep_cols = c("DAid", "NPX"), replace_w_na = NULL, remove_na_cols = NULL)
  expected <- test_data1 |>
    dplyr::select(DAid, NPX)
  expect_equal(result, expected)
})


test_that("clean_data excludes specified assays", {
  set.seed(123)
  result <- clean_data(example_data, filter_assays = c("AARSD1", "ABL1"))
  expected <- example_data |>
    dplyr::filter(!(Assay %in% c("AARSD1", "ABL1"))) |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


# Test clean_metadata ----------------------------------------------------------
test_that("clean_metadata selects specified columns", {
  result <- clean_metadata(example_metadata, keep_cols = c("DAid", "Age"))
  expected <- example_metadata |>
    dplyr::select(DAid, Age)
  expect_equal(result, expected)
})


test_that("The NA columns are removed", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_metadata <- example_metadata
  test_metadata$DAid[random_indices] <- NA

  result <- clean_metadata(test_metadata, keep_cols = c("DAid", "Age"))

  expected <- test_metadata |>
    dplyr::select(DAid, Age) |>
    dplyr::filter(!is.na(DAid))

  result <- as.character(unclass(result))
  expected <- as.character(unclass(expected))

  expect_equal(result, expected)
})


test_that("clean_metadata handles non-existent columns gracefully", {
  result <- clean_metadata(example_metadata, keep_cols = c("DAid", "Age", "BMI", "Height"))
  expected <- example_metadata |>
    dplyr::select(any_of(c("DAid", "Age", "BMI")))
  expect_equal(result, expected)
})


test_that("The specified values are replaced with NAs", {
  result <- clean_metadata(example_metadata, keep_cols = c("DAid", "Disease"), remove_na_cols = NULL)
  expected <- example_metadata |>
    dplyr::select(DAid, Disease)
  expect_equal(result, expected)
})


# Test generate_df -------------------------------------------------------------
test_that("generate_df returns a list of a wide and a joined dataframe", {
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX)
  test_metadata <- example_metadata |>
    dplyr::select(DAid, Disease, Sex, Age, BMI)

  expected_wide <- test_data |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  expected_join <- expected_wide |>
    dplyr::left_join(test_metadata, by = "DAid") |>
    dplyr::relocate("DAid", "Disease", "Sex", "Age", "BMI", dplyr::everything())

  result_df <- generate_df(test_data, test_metadata)
  wide_data <- result_df[[1]]
  join_data <- result_df[[2]]

  expect_equal(wide_data, expected_wide)
  expect_equal(join_data, expected_join)
  unlink("data", recursive = TRUE)
})
