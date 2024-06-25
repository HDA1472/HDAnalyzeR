# Test clean_data --------------------------------------------------------------
test_that("clean_data selects specified columns", {
  result <- clean_data(example_data)
  expected <- example_data |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data excludes specified plates", {
  result <- clean_data(example_data, exclude_plates = "Run001")
  expected <- example_data |>
    dplyr::filter(!(PlateID %in% "Run001")) |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data filters by cohort", {
  # Add Cohort column to filter based on
  example_data_cohort <- example_data |>
    dplyr::mutate(Cohort = rep(c("A", "B"), length.out = dplyr::n()))
  result <- clean_data(example_data_cohort, cohort = "A")
  expected <- example_data_cohort |>
    dplyr::filter(Cohort %in% "A") |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data filters by Assay_Warning", {
  result <- clean_data(example_data, filter_assay = "PASS")
  expected <- example_data |>
    dplyr::filter(Assay_Warning == "PASS") |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("clean_data filters by QC_Warning", {
  result <- clean_data(example_data, filter_qc = c("MANUAL_WARN", "PASS"))
  expected <- example_data |>
    dplyr::filter(QC_Warning %in% c("MANUAL_WARN", "PASS")) |>
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
  example_data_cohort <- example_data |>
    dplyr::mutate(Cohort = rep(c("A", "B"), length.out = dplyr::n()))
  result <- clean_data(example_data_cohort, cohort = "A", exclude_plates = "P2",
                       filter_assay = "PASS", filter_qc = c("MANUAL_WARN", "PASS"))
  expected <- example_data_cohort |>
    dplyr::filter(Cohort %in% "A" & !(PlateID %in% "P2") & Assay_Warning == "PASS" & QC_Warning %in% c("MANUAL_WARN", "PASS")) |>
    dplyr::select(DAid, Assay, NPX)
  expect_equal(result, expected)
})


test_that("The NA columns are removed", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data <- example_data
  test_data$DAid[random_indices] <- NA
  suppressWarnings({
    result <- clean_data(test_data, keep_cols = c("DAid", "NPX"))
  })
  expected <- test_data |>
    dplyr::select(DAid, NPX) |>
    dplyr::filter(!is.na(DAid))
  expect_equal(result, expected)
})


test_that("The specified values are replaced with NAs", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_data1 <- example_data
  test_data2 <- example_data
  test_data1$DAid[random_indices] <- 0
  test_data2$DAid[random_indices] <- NA
  result <- clean_data(test_data1, keep_cols = c("DAid", "NPX"),
                       apply_replacement = TRUE, remove_na_cols = NULL)
  expected <- test_data2 |>
    dplyr::select(DAid, NPX)
  expect_equal(result, expected)
})


# Test clean_metadata ----------------------------------------------------------
test_that("clean_metadata selects specified columns", {
  result <- clean_metadata(example_metadata, keep_cols = c("DAid", "Age"))
  expected <- example_metadata |>
    dplyr::select(DAid, Age)
  expect_equal(result, expected)
})


test_that("clean_metadata excludes specified samples", {
  set.seed(123)
  test_metadata <- example_metadata |>
    dplyr::mutate(`Exclude Sample` = ifelse(runif(dplyr::n()) <= 0.05, "yes", "no"))
  result <- clean_metadata(test_metadata, keep_cols = c("DAid", "Age"), exclude_sample = "yes")
  expected <- test_metadata |>
    dplyr::filter(!(`Exclude Sample` %in% "yes")) |>
    dplyr::select(DAid, Age)
  expect_equal(result, expected)
})


test_that("The NA columns are removed", {
  random_indices <- sample(1:nrow(example_metadata), 20)
  test_metadata <- example_metadata
  test_metadata$DAid[random_indices] <- NA
  suppressWarnings({
    result <- clean_metadata(test_metadata, keep_cols = c("DAid", "Age"))
  })
  expected <- test_metadata |>
    dplyr::select(DAid, Age) |>
    dplyr::filter(!is.na(DAid))
  expect_equal(result, expected)
})


test_that("clean_metadata handles non-existent columns gracefully", {
  result <- clean_metadata(example_metadata, keep_cols = c("DAid", "Age", "BMI", "Height"))
  expected <- example_metadata |>
    dplyr::select(any_of(c("DAid", "Age", "BMI")))
  expect_equal(result, expected)
})


test_that("The specified values are replaced with NAs", {
  result <- clean_data(example_metadata, keep_cols = c("DAid", "GROUP"),
                       apply_replacement = TRUE, remove_na_cols = NULL)
  expected <- example_metadata |>
    dplyr::select(DAid, GROUP)
  expect_equal(result, expected)
})


# Test generate_df -------------------------------------------------------------
test_that("generate_df returns a list of a wide and a joined dataframe", {
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX)
  test_metadata <- example_metadata |>
    dplyr::select(DAid, Sex, Age, BMI)

  expected_wide <- test_data |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  expected_join <- expected_wide |> dplyr::left_join(test_metadata, by = "DAid")

  result_df <- generate_df(test_data, test_metadata)
  wide_data <- result_df[[1]]
  join_data <- result_df[[2]]

  expect_equal(wide_data, expected_wide)
  expect_equal(join_data, expected_join)
  unlink("data", recursive = TRUE)
})


test_that("generate_df saves dataframes properly", {
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX)
  test_metadata <- example_metadata |>
    dplyr::select(DAid, Sex, Age, BMI)

  expected_wide <- test_data |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  expected_join <- expected_wide |> dplyr::left_join(test_metadata, by = "DAid")

  result_df <- generate_df(test_data, test_metadata)
  long_data <- import_df("data/processed/data_metadata/long_data.rda")
  wide_data <- import_df("data/processed/data_metadata/wide_data.rda")
  metadata <- import_df("data/processed/data_metadata/metadata.rda")
  join_data <- import_df("data/processed/data_metadata/join_data.rda")

  expect_equal(long_data, test_data)
  expect_equal(wide_data, expected_wide)
  expect_equal(metadata, test_metadata)
  expect_equal(join_data, expected_join)
  unlink("data", recursive = TRUE)
})
