# Test do_limma_de -------------------------------------------------------------
test_that("do_limma_de performs DE properly", {
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX) |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
    dplyr::left_join(example_metadata |>
                       dplyr::select(dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI"))),
                     by = "DAid") |>
    dplyr::select(DAid, Disease, Sex, Age, BMI, 2:11)

  result <- do_limma_de(test_data, "AML") |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 2)))

  expected <- tibble::tibble(
    Assay = c("ABL1", "ACTA2", "ACAN", "ACE2", "ACP6", "ACAA1", "AARSD1", "ACOX1", "ACTN4", "ACP5"),
    logFC = c(0.75, 0.5, -0.38, 0.47, -0.31, -0.16, 0.11, -0.01, 0, 0.02),
    CI.L = c(0.32, 0.22, -0.61, 0.14, -0.59, -0.57, -0.22, -0.3, -0.21, -0.21),
    CI.R = c(1.17, 0.79, -0.15, 0.79, -0.02, 0.25, 0.43, 0.28, 0.21, 0.25),
    AveExpr = c(1.81, 1.61, 0.57, 0.93, 1.13, 1.01, 3.13, 0.5, 0.39, 0.92),
    t = c(3.47, 3.44, -3.28, 2.82, -2.12, -0.78, 0.63, -0.07, -0.03, 0.17),
    P.Value = c(0, 0, 0, 0, 0.03, 0.44, 0.53, 0.95, 0.98, 0.87),
    adj.P.Val = c(0, 0, 0, 0.01, 0.07, 0.73, 0.75, 0.98, 0.98, 0.98),
    B = c(-0.48, -0.61, -1.12, -2.44, -4.1, -5.92, -6.02, -6.26, -6.27, -6.27),
    Disease = rep("AML", 10),
    sig = c("significant up", "significant up", "significant down", "significant up",
            "not significant", "not significant", "not significant", "not significant",
            "not significant", "not significant")
  )

  expect_equal(result, expected)
})


test_that("do_limma_de corrects for Sex", {
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX) |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
    dplyr::left_join(example_metadata |>
                       dplyr::select(dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI"))),
                     by = "DAid") |>
    dplyr::select(DAid, Disease, Sex, Age, BMI, 2:11)

  result <- do_limma_de(test_data, "AML", correct_sex = TRUE) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 2)))

  expected <- tibble::tibble(
    Assay = c("ACTA2", "ACAN", "ABL1", "ACP6", "ACE2", "ACAA1", "AARSD1", "ACOX1", "ACTN4", "ACP5"),
    logFC = c(0.5, -0.39, 0.66, -0.35, 0.33, -0.19, 0.06, -0.04, -0.01, -0.01),
    CI.L = c(0.21, -0.62, 0.25, -0.64, 0.01, -0.59, -0.27, -0.33, -0.22, -0.24),
    CI.R = c(0.79, -0.16, 1.08, -0.06, 0.66, 0.22, 0.39, 0.26, 0.19, 0.22),
    AveExpr = c(1.61, 0.57, 1.81, 1.13, 0.93, 1.01, 3.13, 0.5, 0.39, 0.92),
    t = c(3.41, -3.36, 3.13, -2.39, 2.05, -0.9, 0.37, -0.24, -0.13, -0.08),
    P.Value = c(0, 0, 0, 0.02, 0.04, 0.37, 0.71, 0.81, 0.9, 0.94),
    adj.P.Val = c(0, 0, 0.01, 0.04, 0.08, 0.61, 0.94, 0.94, 0.94, 0.94),
    B = c(-0.68, -0.82, -1.52, -3.47, -4.17, -5.76, -6.08, -6.17, -6.19, -6.21),
    Disease = rep("AML", 10),
    sig = c("significant up", "significant down", "significant up", "significant down",
            "not significant", "not significant", "not significant", "not significant",
            "not significant", "not significant")
  )

  expect_equal(result, expected)
})


# Test do_ttest_de -------------------------------------------------------------
test_that("do_ttest_de performs DE properly", {
  de_res <- matrix(nrow=0, ncol=4)
  colnames(de_res) <- c("Assay", "P.Value", "logFC", "Disease")
  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX) |>
    dplyr::left_join(example_metadata |>
                       dplyr::select(dplyr::any_of(c("DAid", "Disease"))),
                     by = "DAid")

  result <- tibble::as_tibble(do_ttest_de(de_res, test_data, "AML", "ABL1"), .name_repair = "minimal")
  expected <- tibble::tibble(
    Assay = "ABL1",
    P.Value = "0.00137578481974897",
    logFC = "0.902077248474924",
    Disease = "AML"
  )

  expect_equal(result, expected)
})


# Test run_de ------------------------------------------------------------------
test_that("run_de performs DE properly", {
  first_10_unique_assays <- example_data |>
    dplyr::distinct(Assay) |>
    dplyr::slice(1:10) |>
    dplyr::pull(Assay)

  test_data <- example_data |>
    dplyr::select(DAid, Assay, NPX) |>
    dplyr::filter(Assay %in% first_10_unique_assays)

  result <- run_de(test_data, example_metadata, wide = F, volcano = F)
  result <- result$AML |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 2)))

  expected <- tibble::tibble(
    Assay = c("ABL1", "ACTA2", "ACAN", "ACE2", "ACP6", "ACAA1", "AARSD1", "ACOX1", "ACTN4", "ACP5"),
    logFC = c(0.75, 0.5, -0.38, 0.47, -0.31, -0.16, 0.11, -0.01, 0, 0.02),
    CI.L = c(0.32, 0.22, -0.61, 0.14, -0.59, -0.57, -0.22, -0.3, -0.21, -0.21),
    CI.R = c(1.17, 0.79, -0.15, 0.79, -0.02, 0.25, 0.43, 0.28, 0.21, 0.25),
    AveExpr = c(1.81, 1.61, 0.57, 0.93, 1.13, 1.01, 3.13, 0.5, 0.39, 0.92),
    t = c(3.47, 3.44, -3.28, 2.82, -2.12, -0.78, 0.63, -0.07, -0.03, 0.17),
    P.Value = c(0, 0, 0, 0, 0.03, 0.44, 0.53, 0.95, 0.98, 0.87),
    adj.P.Val = c(0, 0, 0, 0.01, 0.07, 0.73, 0.75, 0.98, 0.98, 0.98),
    B = c(-0.48, -0.61, -1.12, -2.44, -4.1, -5.92, -6.02, -6.26, -6.27, -6.27),
    Disease = rep("AML", 10),
    sig = c("significant up", "significant up", "significant down", "significant up",
            "not significant", "not significant", "not significant", "not significant",
            "not significant", "not significant")
  )

  expect_equal(result, expected)
})
