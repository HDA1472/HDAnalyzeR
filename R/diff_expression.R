utils::globalVariables(c("adj.P.Val", "P.Value", "logFC", "sig", "sig.label", "Sex",
                         "Disease", "case", "control", "is_normal", "Count"))
#' Differential expression analysis with limma
#'
#' `do_limma_de()` performs differential expression analysis using limma package.
#' It can correct the results for metadata columns like Sex, Age, or BMI.
#' The output tibble includes the logFC, p-values, as well as the FDR adjusted p-values.
#' The function removes the NAs from the columns that are used to correct for.
#'
#' @param join_data A tibble with the Olink data in wide format joined with metadata.
#' @param disease The disease of interest.
#' @param correct The variables to correct the results with. Default c("Sex", "Age").
#' @param correct_type The type of the variables to correct the results with. Default c("factor", "numeric").
#' @param only_female The female specific diseases. Default is NULL.
#' @param only_male The male specific diseases. Default is NULL.
#' @param pval_lim The p-value limit of significance. Default is 0.05.
#' @param logfc_lim The logFC limit of significance. Default is 0.
#'
#' @return A tibble with the differential expression results.
#' @keywords internal
do_limma_de <- function(join_data,
                        disease,
                        correct = c("Sex", "Age"),
                        correct_type = c("factor", "numeric", "numeric"),
                        only_female = NULL,
                        only_male = NULL,
                        pval_lim = 0.05,
                        logfc_lim = 0) {

  # Filter for Sex if disease is Sex specific
  if(!is.null(only_female) & disease %in% only_female) {
    join_data <- join_data |>
      dplyr::filter(Sex == "F")
  } else {
    join_data <- join_data
  }

  if(!is.null(only_male) & disease %in% only_male) {
    join_data <- join_data |>
      dplyr::filter(Sex == "M")
  } else {
    join_data <- join_data
  }

  join_data <- join_data |>
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c("Disease", correct)), is.na)) |>  # Remove NAs from columns in formula
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


#' Differential expression analysis with limma for continuous variable
#'
#' This function performs differential expression analysis using limma for a continuous variable.
#' The output dataframe includes the logFC, the p-values, as well as the adjusted p-values with FDR.
#' The function removes the NAs from the columns that are used to correct for.
#'
#' @param join_data (tibble). A tibble with the Olink data in wide format joined with metadata.
#' @param variable (character). The variable of interest.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#'
#' @return de_res (tibble). A tibble with the differential expression results.
#' @keywords internal
do_limma_continuous_de <- function(join_data,
                                   variable,
                                   pval_lim = 0.05,
                                   logfc_lim = 0) {

  join_data <- join_data |>
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c(variable)), is.na))  # Remove NAs from columns in formula

  # Design a model
  design <- stats::model.matrix(~0 + join_data[[variable]])

  # Fit linear model to each protein assay
  data_fit <- join_data |>
    dplyr::select(-dplyr::any_of(c("Disease", "Sex", "Age", "BMI"))) |>
    tibble::column_to_rownames("DAid") |>
    t()

  fit <- limma::lmFit(data_fit, design = design, method = "robust", maxit = 10000)


  # Apply empirical Bayes smoothing to the SE
  ebays_fit <- limma::eBayes(fit)

  # Extract DE results
  de_results <- limma::topTable(ebays_fit,
                                n = nrow(ebays_fit$p.value),
                                adjust.method = "fdr",
                                confint = TRUE)

  de_res <- de_results |>
    tibble::as_tibble(rownames = "Assay") |>
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
#' `do_ttest_de()` performs differential expression analysis using t-test.
#' It separates the data in case-control groups, checks for data normality and
#' perform a t-test or Wilcoxon test respectively. It also performs p value FDR adjustment.
#'
#' @param long_data A tibble with the Olink data in long format and the Sex column from metadata.
#' @param disease The disease of interest.
#' @param assays The assays to run the differential expression analysis on.
#' @param only_female The female specific diseases. Default is NULL.
#' @param only_male The male specific diseases. Default is NULL.
#' @param pval_lim The p-value limit of significance. Default is 0.05.
#' @param logfc_lim The logFC limit of significance. Default is 0.
#'
#' @return A tibble with the differential expression results.
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
      dplyr::filter(Sex == "F")
  } else if(!is.null(only_male) & disease %in% only_male) {
    long_data <- long_data |>
      dplyr::filter(Sex == "M")
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
#' `plot_volcano()` creates volcano plots for the differential expression results.
#' It colors and labels the top up and down regulated proteins.
#'
#' @param de_result The differential expression results.
#' @param pval_lim The p-value limit for significance. Default is 0.05.
#' @param logfc_lim The logFC limit for significance. Default is 0.
#' @param top_up_prot The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette The color palette for the plot. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is "diff_exp".
#' @param title The title of the plot or NULL for no title.
#' @param subtitle If the subtitle should be displayed. Default is TRUE.
#'
#' @return A ggplot object with the volcano plot.
#' @keywords internal
plot_volcano <- function(de_result,
                         pval_lim = 0.05,
                         logfc_lim = 0,
                         top_up_prot = 40,
                         top_down_prot = 10,
                         palette = "diff_exp",
                         title = NULL,
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
                                 label = Assay)) +
    ggplot2::geom_point(size = 1, alpha = 0.4)+
    ggrepel::geom_text_repel(data = subset(tab, sig.label == "top significance"), show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = -log10(pval_lim), linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = logfc_lim, linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = -logfc_lim, linetype = 'dashed') +
    ggplot2::labs(color = "Significance")

    if (!is.null(title)) {
      p <- p + ggplot2::ggtitle(label = paste0(title, ""))
    }
    if (isTRUE(subtitle)) {
      p <- p + ggplot2::ggtitle(label = paste0(title, ""),
                                subtitle = paste0("Num significant up = ", num.sig.up,
                                                  "\nNum significant down = ", num.sig.down))
    }
    if (is.null(names(palette))) {
      p <- p + scale_color_hpa(palette)
    } else {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }

  return(p + theme_hpa() + ggplot2::theme(legend.position = "none"))
}


#' Run differential expression analysis with limma
#'
#' `do_limma()` performs differential expression analysis using limma package.
#' It can correct the results for metadata columns like Sex, Age, or BMI.
#' The output tibble includes the logFC, p-values, as well as the FDR adjusted p-values.
#' The function removes the NAs from the columns that are used to correct for.
#' It can generate and save volcano plots.
#'
#' @param olink_data A tibble with the Olink data in wide format.
#' @param metadata A tibble with the metadata.
#' @param correct The variables to correct the results with. Default c("Sex", "Age", "BMI").
#' @param correct_type The type of the variables to correct the results with. Default c("factor", "numeric", "numeric").
#' @param wide If the data is in wide format. Default is TRUE.
#' @param only_female The female specific diseases. Default is NULL.
#' @param only_male The male specific diseases. Default is NULL.
#' @param volcano Generate volcano plots. Default is TRUE.
#' @param pval_lim The p-value limit of significance. Default is 0.05.
#' @param logfc_lim The logFC limit of significance. Default is 0.
#' @param top_up_prot The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette The color palette for the plot. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is "diff_exp".
#' @param subtitle If the subtitle should be displayed. Default is TRUE.
#' @param save Save the volcano plots. Default is FALSE.
#'
#' @return A list with the differential expression results and volcano plots.
#'   - de_results: A list with the differential expression results.
#'   - volcano_plots: A list with the volcano plots.
#' @export
#'
#' @examples
#' de_results <- do_limma(example_data, example_metadata, wide = FALSE)
#'
#' # Results for AML
#' de_results$de_results$AML
#'
#' # Volcano plot for AML
#' de_results$volcano_plots$AML
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
                            function(disease) plot_volcano(de_results[[disease]],
                                                           pval_lim,
                                                           logfc_lim,
                                                           top_up_prot,
                                                           top_down_prot,
                                                           palette,
                                                           disease,
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


#' Run differential expression analysis with limma for continuous variable
#'
#' This function runs differential expression analysis using limma for a continuous variable.
#' It can generate and save volcano plots.
#'
#' @param olink_data (tibble). A tibble with the Olink data in wide format.
#' @param metadata (tibble). A tibble with the metadata.
#' @param variable (character). The variable of interest.
#' @param wide (logical). If the data is in wide format. Default is TRUE.
#' @param volcano (logical). Generate volcano plots. Default is TRUE.
#' @param pval_lim (numeric). The p-value limit for significance. Default is 0.05.
#' @param logfc_lim (numeric). The logFC limit for significance. Default is 0.
#' @param top_up_prot (numeric). The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot (numeric). The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is "diff_exp".
#' @param subtitle (logical). If the subtitle should be displayed. Default is TRUE.
#' @param save (logical). Save the volcano plots. Default is FALSE.
#'
#' @return de_results (list). A list with the differential expression results and volcano plots.
#'   - de_results: A list with the differential expression results.
#'   - volcano_plots: A list with the volcano plots.
#' @export
#'
#' @examples
#' do_limma_continuous(example_data, example_metadata, "Age", wide = FALSE)
do_limma_continuous <- function(olink_data,
                                metadata,
                                variable,
                                wide = T,
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
        metadata |> dplyr::select(dplyr::any_of(c("DAid", variable))),
        by = "DAid")
  } else {
    join_data <- olink_data |> dplyr::left_join(
      metadata |> dplyr::select(dplyr::any_of(c("DAid", variable))),
      by = "DAid")
  }

  # Run differential expression analysis
  de_results <- do_limma_continuous_de(join_data,
                                       variable,
                                       pval_lim,
                                       logfc_lim)

  # Generate (and save) volcano plots
  if (volcano) {
    volcano_plots <- plot_volcano(de_results,
                                  pval_lim,
                                  logfc_lim,
                                  top_up_prot,
                                  top_down_prot,
                                  palette,
                                  variable,
                                  subtitle)

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


#' Run differential expression analysis with t-test
#'
#' `do_ttest()` performs differential expression analysis using t-test.
#' It separates the data in case-control groups, checks for data normality and
#' perform a t-test or Wilcoxon test respectively. It also performs p value FDR adjustment.
#' It can generate and save volcano plots.
#'
#' @param olink_data A tibble with the Olink data in wide format.
#' @param metadata A tibble with the metadata.
#' @param wide If the data is in wide format. Default is TRUE.
#' @param only_female The female specific diseases. Default is NULL.
#' @param only_male The male specific diseases. Default is NULL.
#' @param volcano Generate volcano plots. Default is TRUE.
#' @param pval_lim The p-value limit of significance. Default is 0.05.
#' @param logfc_lim The logFC limit of significance. Default is 0.
#' @param top_up_prot The number of top up regulated proteins to label on the plot. Default is 40.
#' @param top_down_prot The number of top down regulated proteins to label on the plot. Default is 10.
#' @param palette The color palette for the plot. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is "diff_exp".
#' @param subtitle If the subtitle should be displayed. Default is TRUE.
#' @param save Save the volcano plots. Default is FALSE.
#'
#' @return A list with the differential expression results and volcano plots.
#'   - de_results: A list with the differential expression results.
#'   - volcano_plots: A list with the volcano plots.
#' @export
#'
#' @examples
#' de_results <- do_ttest(example_data, example_metadata, wide = FALSE)
#'
#' # Results for AML
#' de_results$de_results$AML
#'
#' # Volcano plot for AML
#' de_results$volcano_plots$AML
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
                            function(disease) plot_volcano(de_results[[disease]],
                                                           pval_lim,
                                                           logfc_lim,
                                                           top_up_prot,
                                                           top_down_prot,
                                                           palette,
                                                           disease,
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


#' Plot summary visualizations for the differential expression results
#'
#' This function creates summary visualizations for the differential expression results.
#' It plots a barplot with the number of significant proteins for each disease.
#' It also creates upset plots both for the significant up and down regulated proteins for each disease.
#'
#' @param de_results A list with the differential expression results.
#' @param disease_palette The color palette for the disease. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is NULL.
#' @param diff_exp_palette The color palette for the differential expression. If it is a character, it should be one of the palettes from `get_hpa_palettes()`. Default is "diff_exp".
#'
#' @return A list containing the following plots:
#'   - de_barplot: A barplot with the number of significant proteins for each disease.
#'   - upset_plot_up: An upset plot with the significant up regulated proteins for each disease.
#'   - upset_plot_down: An upset plot with the significant down regulated proteins for each disease.
#' @export
#'
#' @examples
#' # Run differential expression analysis
#' de_results <- do_limma(example_data, example_metadata, wide = FALSE)
#'
#' # Plot summary visualizations
#' plot_de_summary(de_results)
plot_de_summary <- function(de_results, disease_palette = NULL, diff_exp_palette = "diff_exp") {
  barplot_data <- de_results$de_results |>
    dplyr::bind_rows() |>
    dplyr::mutate(sig = factor(sig, levels = c("not significant", "significant down", "significant up"))) |>
    dplyr::group_by(Disease, sig) |>
    dplyr::summarise(Count = dplyr::n()) |>
    dplyr::ungroup()

  de_barplot <- barplot_data |>
    ggplot2::ggplot(ggplot2::aes(x = Disease, y = Count, fill = sig)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(x = "", y = "Number of proteins", fill = "Significance") +
    theme_hpa(angled = T) +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_text(face = "bold"))

  if (is.null(names(diff_exp_palette)) && !is.null(diff_exp_palette)) {
    de_barplot <- de_barplot + scale_fill_hpa(diff_exp_palette)
  } else if (!is.null(diff_exp_palette)) {
    de_barplot <- de_barplot + ggplot2::scale_fill_manual(values = diff_exp_palette)
  }

  significant_proteins_up <- lapply(names(de_results$de_results), function(disease) {
    significant_proteins_up <- de_results$de_results[[disease]] |>
      dplyr::filter(sig == "significant up") |>
      dplyr::pull(Assay)

  })
  names(significant_proteins_up) <- names(de_results$de_results)

  significant_proteins_down <- lapply(names(de_results$de_results), function(disease) {

    significant_proteins_down <- de_results$de_results[[disease]] |>
      dplyr::filter(sig == "significant down") |>
      dplyr::pull(Assay)

  })
  names(significant_proteins_down) <- names(de_results$de_results)

  significant_proteins <- list("up" = significant_proteins_up, "down" = significant_proteins_down)

  if (is.null(names(disease_palette)) && !is.null(disease_palette)) {
    pal <- get_hpa_palettes()[[disease_palette]]
  } else if (!is.null(disease_palette)) {
    pal <- disease_palette
  } else {
    pal <- "black"
  }

  upset_plot_up <- UpSetR::upset(UpSetR::fromList(significant_proteins$up),
                                 order.by = "freq",
                                 nsets = length(names(de_results$de_results)),
                                 sets.bar.color = pal)

  upset_plot_down <- UpSetR::upset(UpSetR::fromList(significant_proteins$down),
                                   order.by = "freq",
                                   nsets = length(names(de_results$de_results)),
                                   sets.bar.color = pal)

  return(list("de_barplot" = de_barplot,
              "upset_plot_up" = upset_plot_up,
              "upset_plot_down" = upset_plot_down))
}
