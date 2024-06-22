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
