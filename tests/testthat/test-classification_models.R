# Test split_data --------------------------------------------------------------
test_that("split_data splits data properly", {
  test_data <- tibble::tibble(
    Disease = c("A", "A", "B", "A", "A", "A", "B", "B"),
    protein1 = c(1, 1, 1.6, 2, 1, 1, 1.6, 2),
    protein2 = c(0.8, 0.6, 0.5, 0.3, 0.8, 0.6, 0.5, 0.3)
  )
  result <- split_data(test_data)
  train_res <- result$train_set
  test_res <- result$test_set
  train_exp <- tibble::tibble(
    Disease = c("A", "A", "A", "B", "B"),
    protein1 = c(1, 2, 1, 1.6, 2),
    protein2 = c(0.6, 0.3, 0.6, 0.5, 0.3)
  )
  test_exp <- tibble::tibble(
    Disease = c("A", "B", "A"),
    protein1 = c(1, 1.6, 1),
    protein2 = c(0.8, 0.5, 0.8)
  )
  expect_equal(train_res, train_exp)
  expect_equal(test_res, test_exp)
})


# Test filter_sex_specific_disease ---------------------------------------------
test_that("filter_sex_specific_disease filters data properly", {
  test_data <- tibble::tibble(
    Disease = c("A", "A", "B", "A", "A", "A", "B", "B"),
    Sex = c("Female", "Male", "Female", "Male", "Female", "Male", "Female", "Female"),
    protein1 = c(1, 1, 1.6, 2, 1, 1, 1.6, 2),
    protein2 = c(0.8, 0.6, 0.5, 0.3, 0.8, 0.6, 0.5, 0.3)
  )
  result <- filter_sex_specific_disease(test_data, "B", c("A", "B"), only_female = c("B"))
  control_res <- result$control_data
  diseases_res <- result$diseases_subset
  control_exp <- tibble::tibble(
    Disease = c("A", "B", "A", "B", "B"),
    Sex = c("Female", "Female", "Female", "Female", "Female"),
    protein1 = c(1, 1.6, 1, 1.6, 2),
    protein2 = c(0.8, 0.5, 0.8, 0.5, 0.3)
  )
  diseases_exp = c("A", "B")
  expect_equal(control_res, control_res)
  expect_equal(diseases_res, diseases_exp)
})


test_that("filter_sex_specific_disease filters data and diseases vector properly", {
  test_data <- tibble::tibble(
    Disease = c("A", "A", "B", "A", "A", "A", "B", "B", "C"),
    Sex = c("Female", "Male", "Female", "Male", "Female", "Male", "Female", "Female", "Male"),
    protein1 = c(1, 1, 1.6, 2, 1, 1, 1.6, 2, 1),
    protein2 = c(0.8, 0.6, 0.5, 0.3, 0.8, 0.6, 0.5, 0.3, 1)
  )
  result <- filter_sex_specific_disease(test_data, "B", c("A", "B", "C"), only_female = c("B"), only_male = c("C"))
  control_res <- result$control_data
  diseases_res <- result$diseases_subset
  control_exp <- tibble::tibble(
    Disease = c("A", "B", "A", "B", "B"),
    Sex = c("Female", "Female", "Female", "Female", "Female"),
    protein1 = c(1, 1.6, 1, 1.6, 2),
    protein2 = c(0.8, 0.5, 0.8, 0.5, 0.3)
  )
  diseases_exp = c("A", "B")
  expect_equal(control_res, control_res)
  expect_equal(diseases_res, diseases_exp)
})


# Test make_groups -------------------------------------------------------------
test_that("make_groups creates groups properly", {
  test_data <- tibble::tibble(
    Disease = c("A", "A", "B", "A", "A", "A", "B", "B", "C", "C", "C", "C"),
    protein1 = c(1, 1, 1.6, 2, 1, 1, 1.6, 2, 1, 2, 1, 2),
    protein2 = c(0.8, 0.6, 0.5, 0.3, 0.8, 0.6, 0.5, 0.3, 0.1, 0.2, 0.1, 0.2)
  )
  result <- make_groups(test_data, c("A", "B", "C"))
  exp_a <- tibble::tibble(
    Disease = c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
    protein1 = c(1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 2),
    protein2 = c(0.8, 0.6, 0.3, 0.8, 0.6, 0.3, 0.3, 0.3, 0.2, 0.1, 0.2)
  )
  exp_b <- tibble::tibble(
    Disease = c("B", "B", "B", "A", "A", "C", "C"),
    protein1 = c(1.6, 1.6, 2, 1, 2, 1, 2),
    protein2 = c(0.5, 0.5, 0.3, 0.6, 0.3, 0.1, 0.2)
  )
  exp_c <- tibble::tibble(
    Disease = c("C", "C", "C", "C", "A", "A", "B", "B"),
    protein1 = c(1, 2, 1, 2, 1, 1, 2, 1.6),
    protein2 = c(0.1, 0.2, 0.1, 0.2, 0.8, 0.6, 0.3, 0.5)
  )
  expect_equal(result$A, exp_a)
  expect_equal(result$B, exp_b)
  expect_equal(result$C, exp_c)
})


test_that("make_groups creates groups and filters sex specific groups properly", {
  test_data <- tibble::tibble(
    Disease = c("A", "A", "B", "A", "A", "A", "B", "B", "C", "C", "C", "C"),
    Sex = c("Female", "Male", "Female", "Male", "Female", "Male", "Female", "Female", "Male", "Female", "Male", "Female"),
    protein1 = c(1, 1, 1.6, 2, 1, 1, 1.6, 2, 1, 2, 1, 2),
    protein2 = c(0.8, 0.6, 0.5, 0.3, 0.8, 0.6, 0.5, 0.3, 0.1, 0.2, 0.1, 0.2)
  )
  result <- make_groups(test_data, c("A", "B", "C"), only_female = "B")
  exp_a <- tibble::tibble(
    Disease = c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
    Sex = c("Female", "Male", "Male", "Female", "Male", "Female", "Female", "Female", "Female", "Male", "Female"),
    protein1 = c(1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 2),
    protein2 = c(0.8, 0.6, 0.3, 0.8, 0.6, 0.3, 0.3, 0.3, 0.2, 0.1, 0.2)
  )
  exp_b <- tibble::tibble(
    Disease = c("B", "B", "B", "A", "A", "C", "C"),
    Sex = c("Female", "Female", "Female", "Female", "Female", "Female", "Female"),
    protein1 = c(1.6, 1.6, 2, 1, 1, 2, 2),
    protein2 = c(0.5, 0.5, 0.3, 0.8, 0.8, 0.2, 0.2)
  )
  exp_c <- tibble::tibble(
    Disease = c("C", "C", "C", "C", "A", "A", "B", "B"),
    Sex = c("Male", "Female", "Male", "Female", "Female", "Female", "Female", "Female"),
    protein1 = c(1, 2, 1, 2, 1, 1, 1.6, 2),
    protein2 = c(0.1, 0.2, 0.1, 0.2, 0.8, 0.8, 0.5, 0.3)
  )
  expect_equal(result$A, exp_a)
  expect_equal(result$B, exp_b)
  expect_equal(result$C, exp_c)
})
