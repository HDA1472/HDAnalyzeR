utils::globalVariables(c("adj.P.Val", "P.Value", "logFC", "sig", "sig.label", "Sex", "Disease", "case", "control", "is_normal"))
#' Differential expression analysis with limma
#'
#' This function performs differential expression analysis using limma.
#' It can correct the results with sex, age, and BMI.
#' The output dataframe includes the logFC, the p-values, as well as the adjusted p-values with FDR.
#'
#' @param join_data (tibble). A tibble with the Olink data in wide format joined with metadata.
#' @param disease (character). The disease of interest.
#' @param correct (character or vector). The variables to correct the results with. Default c("Sex", "Age", "BMI").
#' @param correct_type (character or vector). The type of the variables to correct the results with. Default c("factor", "numeric", "numeric").
#' @param only_female (character or vector). The female specific diseases. Default is NULL.
#' @param only_male (character or vector). The male specific diseases. Default is NULL.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#'
#' @return de_res (tibble). A tibble with the differential expression results.
#' @keywords internal
do_limma_de <- function(join_data,
                        disease,
                        correct = c("Sex", "Age", "BMI"),
                        correct_type = c("factor", "numeric", "numeric"),
                        only_female = NULL,
                        only_male = NULL,
                        pval_lim = 0.05,
                        logfc_lim = 0) {

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

  if (!is.null(correct)) {
    for (i in 1:length(correct)) {
      if (correct_type[i] == "factor") {
        cofactor = paste("as.factor(", correct[i], ")")
      } else {
        cofactor = correct[i]
      }
      formula <- paste(formula, "+", cofactor)
    }
  }

  design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
  cols <- c("control", "case", correct)
  cols <- cols[!is.null(cols)]
  colnames(design) <- paste(cols)

  contrast <- limma::makeContrasts(Diff = case - control, levels = design)

  # Fit linear model to each protein assay
  data_fit <- join_data |>
    dplyr::select(-dplyr::any_of(c("Disease", "Sex", "Age", "BMI", correct))) |>
    tibble::column_to_rownames("DAid") |>
    t()

  fit <- limma::lmFit(data_fit, design = design, method = "robust", maxit = 10000)
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
    dplyr::mutate(sig = dplyr::case_when(
      adj.P.Val < pval_lim & logFC < -logfc_lim ~ "significant down",
      adj.P.Val < pval_lim & logFC > logfc_lim ~ "significant up",
      T ~ "not significant")
    ) |>
    dplyr::arrange(adj.P.Val)

  return(de_res)
}


#' Differential expression analysis with t-test
#'
#' This function performs differential expression analysis using t-test.
#'
#' @param long_data (tibble). A tibble with the Olink data in long format and the Sex column from metadata.
#' @param disease (character). The disease of interest.
#' @param assays (character or vector). The assays to run the differential expression analysis on.
#' @param only_female (character or vector). The female specific diseases. Default is NULL.
#' @param only_male (character or vector). The male specific diseases. Default is NULL.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#'
#' @return de_res (tibble). A tibble with the differential expression results.
#' @keywords internal
do_ttest_de <- function(long_data,
                        disease,
                        assays,
                        normality_res,
                        only_female = NULL,
                        only_male = NULL,
                        pval_lim = 0.05,
                        logfc_lim = 0) {

  de_res <- matrix(nrow=0, ncol=4)
  colnames(de_res) <- c("Assay", "P.Value", "logFC", "Disease")

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

  # Run statistical test for each assay
  i <- 0
  de_res_list <- lapply(assays, function(assay) {
    i <<- i + 1
    case_group <- long_data |>
      dplyr::filter(Disease == disease, Assay == assay) |>
      dplyr::pull(NPX)

    control_group <- long_data |>
      dplyr::filter(Disease != disease, Assay == assay) |>
      dplyr::pull(NPX)

    if (normality_res[i] == T) {
      test_res <- stats::t.test(case_group, control_group)
    } else {
      test_res <- stats::wilcox.test(case_group, control_group)
    }

    p.val <- test_res$p.value
    difference <- mean(case_group, na.rm = T) - mean(control_group, na.rm = T)

    de_res <- rbind(de_res, c(assay, p.val, difference, disease))

  })

  combined_de_res <- do.call(rbind, de_res_list)
  combined_de_res <- as.data.frame(combined_de_res) |>
    dplyr::mutate(P.Value = as.numeric(P.Value),
                  logFC = as.numeric(logFC))
  combined_de_res$adj.P.Val <- stats::p.adjust(combined_de_res$P.Value, method = "fdr")
  de_res <- combined_de_res |>
    dplyr::mutate(sig = dplyr::case_when(
      adj.P.Val < pval_lim & logFC < -logfc_lim ~ "significant down",
      adj.P.Val < pval_lim & logFC > logfc_lim ~ "significant up",
      TRUE ~ "not significant"
    )) |>
    dplyr::arrange(adj.P.Val)

  de_res <- tibble::as_tibble(de_res)

  return(de_res)
}


#' Create volcano plots
#'
#' This function creates volcano plots for the differential expression results.
#'
#' @param disease (character). The disease of interest.
#' @param de_result (tibble). The differential expression results.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#' @param top_up_prot (numeric). The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot (numeric). The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is "diff_exp".
#' @param subtitle (logical). If the subtitle should be displayed. Default is TRUE.
#'
#' @return p (plot). A ggplot object with the volcano plot.
#' @keywords internal
create_volcano <- function(disease,
                           de_result,
                           pval_lim = 0.05,
                           logfc_lim = 0,
                           top_up_prot = 40,
                           top_down_prot = 10,
                           palette = "diff_exp",
                           subtitle = T) {

  top.sig.down <- de_result |>
    dplyr::filter(adj.P.Val < pval_lim & logFC < -logfc_lim) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::pull(Assay)

  top.sig.up <- de_result |>
    dplyr::filter(adj.P.Val < pval_lim & logFC > logfc_lim) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::pull(Assay)

  top.sig.prot <- c(top.sig.up[1:top_up_prot], top.sig.down[1:top_down_prot])

  tab <- de_result |>
    dplyr::mutate(sig.label = ifelse(Assay %in% top.sig.prot, "top significance", 0))

  num.sig.up <- length(top.sig.up)
  num.sig.down <- length(top.sig.down)

  p <- de_result |>
    ggplot2::ggplot(ggplot2::aes(x = logFC,
                                 y = -log10(adj.P.Val),
                                 color = sig,
                                 label = Assay
    )) +
    ggplot2::geom_point(size = 1, alpha = 0.4)+
    ggrepel::geom_text_repel(data = subset(tab, sig.label == "top significance")) +
    ggplot2::geom_hline(yintercept = -log10(pval_lim), linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = logfc_lim, linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = -logfc_lim, linetype = 'dashed') +
    ggplot2::ggtitle(label = paste0(disease, ""),
                     subtitle = paste0("Num significant up = ", num.sig.up,
                                       "\nNum significant down = ", num.sig.down)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.subtitle = ggplot2::element_text(size = 10, face = "italic"))

    if (isFALSE(subtitle)) {
      p <- p + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
    }

    if (is.null(names(palette))) {
      p <- p + scale_color_hpa(palette)
    } else {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }

  return(p)
}


#' Run differential expression analysis with limma
#'
#' This function runs differential expression analysis using limma.
#' It can generate and save volcano plots.
#'
#' @param olink_data (tibble). A tibble with the Olink data in wide format.
#' @param metadata (tibble). A tibble with the metadata.
#' @param correct (character or vector). The variables to correct the results with. Default c("Sex", "Age", "BMI").
#' @param correct_type (character or vector). The type of the variables to correct the results with. Default c("factor", "numeric", "numeric").
#' @param wide (logical). If the data is in wide format. Default is TRUE.
#' @param only_female (character or vector). The female specific diseases. Default is NULL.
#' @param only_male (character or vector). The male specific diseases. Default is NULL.
#' @param volcano (logical). Generate volcano plots. Default is TRUE.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#' @param top_up_prot (numeric). The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot (numeric). The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is "diff_exp".
#' @param subtitle (logical). If the subtitle should be displayed. Default is TRUE.
#' @param save (logical). Save the volcano plots. Default is FALSE.
#'
#' @return de_results (list). A list with the differential expression results and volcano plots.
#'   - de_results: A list with the differential expression results.
#'   - volcano_plots: A list with the volcano plots.
#' @export
#'
#' @examples
#' test_data <- example_data |> dplyr::select(DAid, Assay, NPX)
#' do_limma(test_data, example_metadata, wide = FALSE)
do_limma <- function(olink_data,
                     metadata,
                     correct = c("Sex", "Age", "BMI"),
                     correct_type = c("factor", "numeric", "numeric"),
                     wide = T,
                     only_female = NULL,
                     only_male = NULL,
                     volcano = T,
                     pval_lim = 0.05,
                     logfc_lim = 0,
                     top_up_prot = 40,
                     top_down_prot = 10,
                     palette = "diff_exp",
                     subtitle = T,
                     save = F) {

  # Prepare Olink data and merge them with metadata
  if (isFALSE(wide)) {
    join_data <- olink_data |>
      dplyr::select(DAid, Assay, NPX) |>
      tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
      dplyr::left_join(
        metadata |> dplyr::select(dplyr::any_of(c("DAid", "Disease", correct))),
        by = "DAid")
  } else {
    join_data <- olink_data |> dplyr::left_join(
      metadata |> dplyr::select(dplyr::any_of(c("DAid", "Disease", correct))),
      by = "DAid")
  }

  # Run differential expression analysis
  levels <- unique(join_data$Disease)

  de_results <- lapply(levels,
                       function(disease) do_limma_de(join_data,
                                                     disease,
                                                     correct,
                                                     correct_type,
                                                     only_female,
                                                     only_male,
                                                     pval_lim,
                                                     logfc_lim))


  names(de_results) <- levels

  # Generate (and save) volcano plots
  if (volcano) {

    volcano_plots <- lapply(levels,
                            function(disease) create_volcano(disease,
                                                             de_results[[disease]],
                                                             pval_lim,
                                                             logfc_lim,
                                                             top_up_prot,
                                                             top_down_prot,
                                                             palette,
                                                             subtitle))
    names(volcano_plots) <- levels

    if (isTRUE(save)) {
      dir_name <- create_dir("results/volcano_plots", date = T)
      for (i in 1:length(levels)) {
        ggplot2::ggsave(volcano_plots[[i]], filename = paste0(dir_name, "/", levels[i], "_volcano.png"), width = 10, height = 8)
      }
    }
    return(list("limma_results" = de_results, "volcano_plots" = volcano_plots))
  }
  return(de_results)
}


#' Run differential expression analysis with t-test
#'
#' This function runs differential expression analysis using t-test.
#' It can generate and save volcano plots.
#'
#' @param olink_data (tibble). A tibble with the Olink data in wide format.
#' @param metadata (tibble). A tibble with the metadata.
#' @param wide (logical). If the data is in wide format. Default is TRUE.
#' @param only_female (character or vector). The female specific diseases. Default is NULL.
#' @param only_male (character or vector). The male specific diseases. Default is NULL.
#' @param volcano (logical). Generate volcano plots. Default is TRUE.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#' @param top_up_prot (numeric). The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot (numeric). The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is "diff_exp".
#' @param subtitle (logical). If the subtitle should be displayed. Default is TRUE.
#' @param save (logical). Save the volcano plots. Default is FALSE.
#'
#' @return de_results (list). A list with the differential expression results and volcano plots.
#'   - de_results: A list with the differential expression results.
#'   - volcano_plots: A list with the volcano plots.
#' @export
#'
#' @examples
#' test_data <- example_data |> dplyr::select(DAid, Assay, NPX)
#' do_ttest(test_data, example_metadata, wide = FALSE)
do_ttest <- function(olink_data,
                     metadata,
                     wide = T,
                     only_female = NULL,
                     only_male = NULL,
                     volcano = T,
                     pval_lim = 0.05,
                     logfc_lim = 0,
                     top_up_prot = 40,
                     top_down_prot = 10,
                     palette = "diff_exp",
                     subtitle = T,
                     save = F) {

  # Prepare Olink data and merge them with metadata
  if (isFALSE(wide)) {
    join_data <- olink_data |>
      dplyr::select(DAid, Assay, NPX) |>
      tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
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

  long_data <- join_data |>
    dplyr::select(-dplyr::any_of(c("Age", "BMI"))) |>
    tidyr::pivot_longer(!c("DAid", "Disease", "Sex"), names_to = "Assay", values_to = "NPX")

  assays <- unique(long_data$Assay)

  normality_res <- check_normality(
    join_data |>
      dplyr::select(-dplyr::any_of(c("DAid", "Disease", "Sex", "Age", "BMI")))
  ) |>
    dplyr::pull(is_normal)

  de_results <- lapply(levels, function(disease) {

    de_res <- do_ttest_de(long_data,
                          disease,
                          assays,
                          normality_res,
                          only_female,
                          only_male,
                          pval_lim,
                          logfc_lim)

  })

  names(de_results) <- levels

  # Generate (and save) volcano plots
  if (volcano) {

    volcano_plots <- lapply(levels,
                            function(disease) create_volcano(disease,
                                                             de_results[[disease]],
                                                             pval_lim,
                                                             logfc_lim,
                                                             top_up_prot,
                                                             top_down_prot,
                                                             palette,
                                                             subtitle))
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
