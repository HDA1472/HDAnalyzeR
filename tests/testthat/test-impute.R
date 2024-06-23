# Test impute_knn --------------------------------------------------------------
test_that("impute_knn imputes values in the proper way", {
  dummy_data <- tibble::tibble(
    A = c(80, 44, NA, 50, 29),  # The 2 neighbors of the 1st row are 44 and 50
    B = c(30, NA, 85, 70, 54),  # The 2 neighbors of the 2nd row are 30 and 85
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  result <- impute_knn(dummy_data, k = 2)
  expected <- tibble::tibble(
    A = c(80, 44, 47, 50, 29),
    B = c(30.0, 57.5, 85.0, 70.0, 54.0),
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  expect_equal(result, expected)
})


test_that("impute_knn handles columns excluded from imputation", {
  dummy_data <- tibble::tibble(
    ID = c(1, 2, 3, 4, 5),
    A = c(80, 44, NA, 50, 29),  # The 2 neighbors of the 1st row are 44 and 50
    B = c(30, NA, 85, 70, 54),  # The 2 neighbors of the 2nd row are 30 and 85
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  result <- impute_knn(dummy_data, k = 2, exclude_cols = "ID")
  expected <- tibble::tibble(
    ID = c(1, 2, 3, 4, 5),
    A = c(80, 44, 47, 50, 29),
    B = c(30.0, 57.5, 85.0, 70.0, 54.0),
    C = c(7, 10, 25, 74, 49),
    D = c(14, 0, 5, 9, 20)
  )
  expect_equal(result, expected)
})
