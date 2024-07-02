utils::globalVariables(c("adj.P.Val", "P.Value", "logFC", "sig", "sig.label", "Sex", "Disease", "case", "control"))
#' Differential expression analysis with limma
#'
#' This function performs differential expression analysis using limma.
#' It can correct the results with sex, age, and BMI.
#' The output dataframe includes the logFC, the p-values, as well as the adjusted p-values with FDR.
#'
#' @param join_data (tibble). A tibble with the Olink data in wide format joined with metadata.
#' @param disease (character). The disease of interest
#' @param correct_sex (logical). Correct the results with Sex. Default is FALSE
#' @param correct_age (logical). Correct the results with Age. Default is FALSE
#' @param correct_bmi (logical). Correct the results with BMI. Default is FALSE
#' @param only_female (character or vector). The female specific diseases. Default is NULL
#' @param only_male (character or vector). The male specific diseases. Default is NULL
#'
#' @return de_res (tibble). A tibble with the differential expression results
#' @export
#'
#' @examples
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
#'   dplyr::left_join(example_metadata |>
#'                      dplyr::select(dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI"))),
#'                    by = "DAid")
#'
#' do_limma_de(test_data, "AML")
do_limma_de <- function(join_data,
                        disease,
                        correct_sex = F,
                        correct_age = F,
                        correct_bmi = F,
                        only_female = NULL,
                        only_male = NULL) {

  # Filter for Sex if disease is Sex specific
  if(!is.null(only_female) & disease %in% only_female) {
    join_data <- join_data |>
      dplyr::filter(Sex == "Female")
  } else {
    join_data <- join_data
  }

  if(!is.null(only_male) & disease %in% only_male) {
    join_data <- join_data |>
      dplyr::filter(Sex == "Male")
  } else {
    join_data <- join_data
  }

  join_data <- join_data |>
    dplyr::mutate(Disease = ifelse(Disease == disease, "1_Case", "0_Control"))

  # Design a model - add Disease, and Sex, Age, BMI
  formula <- "~0 + as.factor(Disease)"

  if (correct_sex) {
    formula <- paste(formula, "+ as.factor(Sex)")
  }
  if (correct_age) {
    formula <- paste(formula, "+ Age")
  }
  if (correct_bmi) {
    formula <- paste(formula, "+ BMI")
  }

  design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
  colnames(design) <- c("control", "case", if (correct_sex) "Sex", if (correct_age) "Age", if (correct_bmi) "BMI")

  contrast <- limma::makeContrasts(Diff = case - control, levels = design)

  # Fit linear model to each protein assay
  dat_fit <- join_data |>
    dplyr::select(-dplyr::any_of(c("Disease", "Sex", "Age", "BMI"))) |>
    tibble::column_to_rownames("DAid") |>
    t()

  fit <- limma::lmFit(dat_fit, design = design, method = "robust", maxit = 10000)

  contrast_fit <- limma::contrasts.fit(fit, contrast)

  # Apply empirical Bayes smoothing to the SE
  ebays_fit <- limma::eBayes(contrast_fit)

  # Extract DE results
  de_results <- limma::topTable(ebays_fit,
                                n = nrow(ebays_fit$p.value),
                                adjust.method = "fdr",
                                confint = TRUE)

  de_res <- de_results |>
    tibble::as_tibble(rownames = "Assay") |>
    dplyr::mutate(Disease = disease) |>
    dplyr::mutate(sig = dplyr::case_when(adj.P.Val < 0.05 & logFC < 0 ~ "significant down",
                                         adj.P.Val < 0.05 & logFC > 0 ~ "significant up",
                                         T ~ "not significant"))

  return(de_res)
}


#' Differential expression analysis with t-test
#'
#' This function performs differential expression analysis using t-test.
#'
#' @param de_res (matrix). A matrix to store the results
#' @param long_data (tibble). A tibble with the Olink data in long format.
#' @param disease (character). The disease of interest
#' @param assay (character). The assay of interest
#' @param only_female (character or vector). The female specific diseases. Default is NULL
#' @param only_male (character or vector). The male specific diseases. Default is NULL
#'
#' @return de_res(matrix). A matrix with the differential expression results
#' @export
#'
#' @examples
#' de_res <- matrix(nrow=0, ncol=4)
#' colnames(de_res) <- c("Assay", "P.Value", "logFC", "Disease")
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   dplyr::left_join(example_metadata |>
#'                      dplyr::select(dplyr::any_of(c("DAid", "Disease"))),
#'                    by = "DAid")
#'
#' do_ttest_de(de_res, test_data, "AML", "ABL1")
do_ttest_de <- function(de_res,
                        long_data,
                        disease,
                        assay,
                        only_female = NULL,
                        only_male = NULL) {

  # Filter for Sex if disease is Sex specific
  if(!is.null(only_female) & disease %in% only_female) {
    long_data <- long_data |>
      dplyr::filter(Sex == "Female")
  } else if(!is.null(only_male) & disease %in% only_male) {
    long_data <- long_data |>
      dplyr::filter(Sex == "Male")
  } else {
    long_data <- long_data
  }

  case_group <- long_data |>
    dplyr::filter(Disease == disease, Assay == assay) |>
    dplyr::pull(NPX)

  control_group <- long_data |>
    dplyr::filter(Disease != disease, Assay == assay) |>
    dplyr::pull(NPX)

  test_res <- stats::t.test(case_group, control_group)
  p.val <- test_res$p.value
  difference <- mean(case_group, na.rm = T) - mean(control_group, na.rm = T)

  de_res <- rbind(de_res, c(assay, p.val, difference, disease))

  return(de_res)
}


#' Create volcano plots
#'
#' This function creates volcano plots for the differential expression results.
#'
#' @param disease (character). The disease of interest
#' @param de_result (tibble). The differential expression results
#'
#' @return p (plot). A ggplot object with the volcano plot
#' @export
#'
#' @examples
#' test_data <- tibble::tibble(
#'   Assay = c("AARSD1", "ABL1", "ACAA1", "ACAN", "ACE2"),
#'   logFC = c(0.5, -1.2, 0.3, -0.8, 0.6),
#'   adj.P.Val = c(0.01, 0.03, 0.2, 0.05, 0.001),
#'   sig = c("significant up",
#'           "significant down",
#'           "not significant",
#'           "significant down",
#'           "significant up")
#' )
#' create_volcano("AML", test_data)
create_volcano <- function(disease, de_result) {
  top.sig.down <- de_result |>
    dplyr::filter(adj.P.Val < 0.05 & logFC < 0) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::pull(Assay)

  top.sig.up <- de_result |>
    dplyr::filter(adj.P.Val < 0.05 & logFC > 0) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::pull(Assay)

  top.sig.prot <- c(top.sig.up[1:40], top.sig.down[1:10])

  tab <- de_result |>
    dplyr::mutate(sig.label = ifelse(Assay %in% top.sig.prot, "top sig", 0))

  num.sig.up <- length(top.sig.up)
  num.sig.down <- length(top.sig.down)

  p <- de_result |>
    ggplot2::ggplot(ggplot2::aes(x = logFC,
                                 y = -log10(adj.P.Val),
                                 color = sig,
                                 label = Assay
    )) +
    ggplot2::geom_point(size = 1, alpha = 0.4)+
    ggrepel::geom_text_repel(data = subset(tab, sig.label == "top sig")) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
    ggplot2::scale_color_manual(values = c("grey", "blue", "red")) +
    ggplot2::ggtitle(label = paste0(disease, ""),
                     subtitle = paste0("Num significant up = ", num.sig.up,
                                       "\nNum significant down = ", num.sig.down)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.subtitle = ggplot2::element_text(size = 10, face = "italic"))

  return(p)
}


#' Run differential expression analysis
#'
#' This function runs differential expression analysis using limma or t-test.
#' It can generate and save volcano plots.
#'
#' @param olink_data (tibble). A tibble with the Olink data in wide format.
#' @param metadata (tibble). A tibble with the metadata
#' @param correct_sex (logical). Correct the results with Sex. Default is FALSE
#' @param correct_age (logical). Correct the results with Age. Default is FALSE
#' @param correct_bmi (logical). Correct the results with BMI. Default is FALSE
#' @param wide (logical). If the data is in wide format. Default is TRUE
#' @param only_female (character or vector). The female specific diseases. Default is NULL
#' @param only_male (character or vector). The male specific diseases. Default is NULL
#' @param stat (character). The statistical test to use (limma or ttest). Default is "limma"
#' @param volcano (logical). Generate volcano plots. Default is TRUE
#' @param save (logical). Save the volcano plots. Default is FALSE
#'
#' @return de_results (list). A list with the differential expression results and volcano plots
#'   - de_results: A list with the differential expression results
#'   - volcano_plots: A list with the volcano plots
#' @export
#'
#' @examples
#' test_data <- example_data |> dplyr::select(DAid, Assay, NPX)
#' run_de(test_data, example_metadata, wide = FALSE)
run_de <- function(olink_data,
                   metadata,
                   correct_sex = F,
                   correct_age = F,
                   correct_bmi = F,
                   wide = T,
                   only_female = NULL,
                   only_male = NULL,
                   stat = "limma",
                   volcano = T,
                   save = F) {

  # Prepare Olink data and merge them with metadata
  if (isFALSE(wide)) {
    join_data <- olink_data |>
      dplyr::select(1:3) |>
      tidyr::pivot_wider(names_from = 2, values_from = 3) |>
      dplyr::left_join(
        metadata |> dplyr::select(dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI"))),
        by = "DAid")
  } else {
    join_data <- olink_data |> dplyr::left_join(
      metadata |> dplyr::select(dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI"))),
      by = "DAid")
  }

  # Run differential expression analysis
  levels <- unique(join_data$Disease)
  if (stat == "limma") {
    de_results <- lapply(levels,
                         function(disease) do_limma_de(join_data,
                                                       disease,
                                                       correct_sex,
                                                       correct_age,
                                                       correct_bmi,
                                                       only_female,
                                                       only_male))
  } else if (stat == "ttest") {
    long_data <- join_data |>
      dplyr::select(-dplyr::any_of(c("Age", "BMI"))) |>
      tidyr::pivot_longer(!c("DAid", "Disease", "Sex"), names_to = "Assay", values_to = "NPX")

    assays <- unique(long_data$Assay)

    de_results <- lapply(levels, function(disease) {
      de_res <- matrix(nrow=0, ncol=4)
      colnames(de_res) <- c("Assay", "P.Value", "logFC", "Disease")
      de_res_list <- lapply(assays, function(assay) {
        de_res <- do_ttest_de(de_res,
                              long_data,
                              disease,
                              assay,
                              only_female,
                              only_male)
      })

      combined_de_res <- do.call(rbind, de_res_list)
      combined_de_res <- as.data.frame(combined_de_res) |>
        dplyr::mutate(P.Value = as.numeric(P.Value),
                      logFC = as.numeric(logFC))
      combined_de_res$adj.P.Val <- stats::p.adjust(combined_de_res$P.Value, method = "fdr")
      combined_de_res <- combined_de_res |>
        dplyr::mutate(sig = dplyr::case_when(
          adj.P.Val < 0.05 & logFC < 0 ~ "significant down",
          adj.P.Val < 0.05 & logFC > 0 ~ "significant up",
          TRUE ~ "not significant"
        ))
    })
  }
  names(de_results) <- levels

  # Generate (and save) volcano plots
  if (volcano) {

    volcano_plots <- lapply(levels,
                            function(disease) create_volcano(disease, de_results[[disease]]))
    names(volcano_plots) <- levels

    if (isTRUE(save)) {
      dir_name <- create_dir("results/volcano_plots", date = T)
      for (i in 1:length(levels)) {
        ggplot2::ggsave(volcano_plots[[i]], filename = paste0(dir_name, "/", levels[i], "_volcano.png"), width = 10, height = 8)
      }
    }
    return(list("de_results" = de_results, "volcano_plots" = volcano_plots))
  }
  return(de_results)
}
