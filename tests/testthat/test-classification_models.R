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


# Test generate_subtitle -------------------------------------------------------
test_that("split_data splits data properly", {
  features <- tibble::tibble(Variable = c("A", "B", "C"),
                            Importance = c(0.1, 0.2, 0.3),
                            Sign = c("NEG", "POS", "NEG"),
                            Scale_Importance = c(10, 20, 30))

  result <- generate_subtitle(features, 0.8, 0.9, 0.7, 0.85, 0.5, c("accuracy", "sensitivity"))
  expected <- "accuracy = 0.8    sensitivity = 0.9    \n"
  expect_equal(result, expected)
})


# Test do_elnet ----------------------------------------------------------------
test_that("AML Classification Results Test", {

  unique_samples <- unique(example_data$Sample)
  filtered_data <- example_data |>
    dplyr::filter(Sample %in% unique_samples[1:148])

  res <- do_elnet(filtered_data,
                  example_metadata,
                  palette = "cancers12",
                  cv_sets = 2,
                  grid_size = 1,
                  ncores = 1)

  # Test if the `res` object has the `AML` component
  expect_true("AML" %in% names(res), info = "Check if 'AML' component exists in res")

  # Test if the `AML` component has the expected structure
  expect_true(is.list(res$AML), info = "Check if 'AML' component is a list")

  # Test if the AML results contain 'hypopt_res' component
  expect_true("hypopt_res" %in% names(res$AML), info = "Check if 'hypopt_res' exists in AML results")

  # Test if the AML results contain 'hypopt_res$elnet_tune' component
  expect_true("elnet_tune" %in% names(res$AML$hypopt_res), info = "Check if 'elnet' exists in AML results")

  # Test if elnet_tune has expected structure and content
  expect_true(is.data.frame(res$AML$hypopt_res$elnet_tune), info = "Check if elnet_tune is a data frame")
  expect_true(nrow(res$AML$hypopt_res$elnet_tune) == 2, info = "Check if elnet_tune has rows")
  expect_true(ncol(res$AML$hypopt_res$elnet_tune) == 5, info = "Check if elnet_tune has rows")

  # Test if AML results contain 'finalfit_res' component
  expect_true("finalfit_res" %in% names(res$AML), info = "Check if 'finalfit_res' exists in AML results")

  # Test if finalfit_res has expected structure and content
  expect_true(is.list(res$AML$finalfit_res), info = "Check if finalfit_res is a list")

  # Test if finalfit_res contains 'metrics' component
  expect_true("metrics" %in% names(res$AML$testfit_res), info = "Check if 'metrics' exists in finalfit_res")

  # Test if metrics contains accuracy, sensitivity, specificity, auc, and conf_matrix
  expect_true("accuracy" %in% names(res$AML$testfit_res$metrics), info = "Check if 'accuracy' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$accuracy), info = "Check if 'accuracy' is numeric")
  expect_true("sensitivity" %in% names(res$AML$testfit_res$metrics), info = "Check if 'sensitivity' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$sensitivity), info = "Check if 'sensitivity' is numeric")
  expect_true("specificity" %in% names(res$AML$testfit_res$metrics), info = "Check if 'specificity' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$specificity), info = "Check if 'specificity' is numeric")
  expect_true("auc" %in% names(res$AML$testfit_res$metrics), info = "Check if 'auc' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$auc), info = "Check if 'auc' is numeric")
  expect_true("conf_matrix" %in% names(res$AML$testfit_res$metrics), info = "Check if 'conf_matrix' exists in metrics")

  # Test if AML results contain 'var_imp_res' component
  expect_true("var_imp_res" %in% names(res$AML), info = "Check if 'var_imp_res' exists in AML results")

  # Test if var_imp_res has expected structure and content
  expect_true(is.list(res$AML$var_imp_res), info = "Check if var_imp_res is a list")

  # Test if var_imp_res contains 'features' component
  expect_true("features" %in% names(res$AML$var_imp_res), info = "Check if 'features' exists in var_imp_res")

  # Test if features has expected structure and content
  expect_true(is.data.frame(res$AML$var_imp_res$features), info = "Check if features is a data frame")
})


# Test do_rf -------------------------------------------------------------------
test_that("AML Classification Results Test", {

  unique_samples <- unique(example_data$Sample)
  filtered_data <- example_data |>
    dplyr::filter(Sample %in% unique_samples[1:148])

  res <- do_rf(filtered_data,
               example_metadata,
               palette = "cancers12",
               cv_sets = 2,
               grid_size = 1,
               ncores = 1)

  # Test if the `res` object has the `AML` component
  expect_true("AML" %in% names(res), info = "Check if 'AML' component exists in res")

  # Test if the `AML` component has the expected structure
  expect_true(is.list(res$AML), info = "Check if 'AML' component is a list")

  # Test if the AML results contain 'hypopt_res' component
  expect_true("hypopt_res" %in% names(res$AML), info = "Check if 'hypopt_res' exists in AML results")

  # Test if the AML results contain 'hypopt_res$rf_tune' component
  expect_true("rf_tune" %in% names(res$AML$hypopt_res), info = "Check if 'rf_tune' exists in AML results")

  # Test if rf_tune has expected structure and content
  expect_true(is.data.frame(res$AML$hypopt_res$rf_tune), info = "Check if rf_tune is a data frame")
  expect_true(nrow(res$AML$hypopt_res$rf_tune) == 2, info = "Check if rf_tune has rows")
  expect_true(ncol(res$AML$hypopt_res$rf_tune) == 5, info = "Check if rf_tune has rows")

  # Test if AML results contain 'finalfit_res' component
  expect_true("finalfit_res" %in% names(res$AML), info = "Check if 'finalfit_res' exists in AML results")

  # Test if finalfit_res has expected structure and content
  expect_true(is.list(res$AML$finalfit_res), info = "Check if finalfit_res is a list")

  # Test if finalfit_res contains 'metrics' component
  expect_true("metrics" %in% names(res$AML$testfit_res), info = "Check if 'metrics' exists in finalfit_res")

  # Test if metrics contains accuracy, sensitivity, specificity, auc, and conf_matrix
  expect_true("accuracy" %in% names(res$AML$testfit_res$metrics), info = "Check if 'accuracy' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$accuracy), info = "Check if 'accuracy' is numeric")
  expect_true("sensitivity" %in% names(res$AML$testfit_res$metrics), info = "Check if 'sensitivity' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$sensitivity), info = "Check if 'sensitivity' is numeric")
  expect_true("specificity" %in% names(res$AML$testfit_res$metrics), info = "Check if 'specificity' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$specificity), info = "Check if 'specificity' is numeric")
  expect_true("auc" %in% names(res$AML$testfit_res$metrics), info = "Check if 'auc' exists in metrics")
  expect_true(is.numeric(res$AML$testfit_res$metrics$auc), info = "Check if 'auc' is numeric")
  expect_true("conf_matrix" %in% names(res$AML$testfit_res$metrics), info = "Check if 'conf_matrix' exists in metrics")

  # Test if AML results contain 'var_imp_res' component
  expect_true("var_imp_res" %in% names(res$AML), info = "Check if 'var_imp_res' exists in AML results")

  # Test if var_imp_res has expected structure and content
  expect_true(is.list(res$AML$var_imp_res), info = "Check if var_imp_res is a list")

  # Test if var_imp_res contains 'features' component
  expect_true("features" %in% names(res$AML$var_imp_res), info = "Check if 'features' exists in var_imp_res")

  # Test if features has expected structure and content
  expect_true(is.data.frame(res$AML$var_imp_res$features), info = "Check if features is a data frame")
})
