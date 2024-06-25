# Test calc_na_percentages_col -------------------------------------------------
test_that("calc_na_percentage_col calculates NA percentages", {
  result <- calc_na_percentage_col(example_metadata)
  expected <- tibble::tibble(
    column = c("Grade"),
    na_percentage = c(91.5)
  )
  expect_equal(result, expected)
})


test_that("calc_na_percentage_col handles dataframe with no NAs", {
  result <- calc_na_percentage_col(example_data)
  expected <- tibble::tibble(
    column = character(),
    na_percentage = numeric()
  )
  expect_equal(result, expected)
})


# Test calc_na_percentages_row -------------------------------------------------
test_that("calc_na_percentage_row calculates NA percentages", {
  test_data <- tibble::tibble(
    DAid = c("1", "2", "3", "4", "5"),
    A = c(80, 44, NA, 50, 29),
    B = c(30, NA, 85, 70, 54),
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  result <- calc_na_percentage_row(test_data)
  expected <- tibble::tibble(
    DAid = c("2", "3"),
    na_percentage = c(20.0, 20.0)
  )
  expect_equal(result, expected)
})


test_that("calc_na_percentage_row handles dataframe with no NAs", {
  test_data <- tibble::tibble(
    DAid = c("1", "2", "3", "4", "5"),
    A = c(80, 44, 6, 50, 29),
    B = c(30, 5, 85, 70, 54),
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  result <- calc_na_percentage_row(test_data)
  expected <- tibble::tibble(
    DAid = character(),
    na_percentage = numeric()
  )
  expect_equal(result, expected)
})


# Test qc_summary --------------------------------------------------------------

