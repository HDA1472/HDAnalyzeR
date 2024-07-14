utils::globalVariables(c("terms", "value", "component", "positive"))
#' Create PCA loadings plot
#'
#' This function creates a plot of the PCA loadings for the top-8 features and first 4 PCs.
#'
#' @param tidied_res (tibble). A tibble with the results of the PCA analysis
#'
#' @return (plot). A ggplot object
#' @keywords internal
plot_loadings <- function(tidied_res) {

  p <- tidied_res |>
    dplyr::filter(component %in% paste0("PC", 1:4)) |>
    dplyr::group_by(component) |>
    dplyr::top_n(8, abs(value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(terms = tidytext::reorder_within(terms, abs(value), component)) |>
    dplyr::mutate(positive = factor((value > 0), levels = c(TRUE, FALSE))) |>
    ggplot2::ggplot(ggplot2::aes(abs(value), terms, fill = positive)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap( ~ component, scales = "free_y") +
    tidytext::scale_y_reordered() +
    ggplot2::labs(x = "Absolute value of contribution",
                  y = NULL, fill = "Positive?") +
    ggplot2::theme_classic()

  if (is.null(names(palette))) {
    p <- p + scale_color_hpa(palette)
  } else {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }
  return(p)
}


#' Create a plot of the dimensionality reduction results
#'
#' This function creates a plot of the dimensionality reduction results.
#' The function can be used to plot the results of PCA (1st vs 2nd PC) or UMAP analyses.
#'
#' @param res (tibble). A tibble with the results of the dimensionality reduction analysis
#' @param x (string). The name of the column in `res` that contains the x-axis values
#' @param y (string). The name of the column in `res` that contains the y-axis values
#' @param metadata (tibble). A tibble with metadata information to be added to the plot
#' @param variable (string). The name of the column in `metadata` that contains the variable to be used to plot the points color
#' @param palette (vector). A vector with the colors to be used in the plot
#'
#' @return (plot). A ggplot object
#' @keywords internal
plot_dim_reduction <- function(res, x, y, metadata, variable, palette) {

  if (!is.null(metadata)) {
    p <- res |>
      dplyr::left_join(metadata, by = "DAid") |>
      ggplot2::ggplot(ggplot2::aes(!!rlang::sym(x), !!rlang::sym(y))) +
      ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(variable)),
                          alpha = 0.7,
                          size = 2) +
      ggplot2::labs(color = variable) +
      ggplot2::theme_classic()
  } else {
    p <- res |>
      ggplot2::ggplot(ggplot2::aes(!!rlang::sym(x), !!rlang::sym(y))) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::theme_classic()
  }

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }

  return(p)
}


#' Perform PCA analysis
#'
#' This function performs a PCA analysis on the data provided.
#' The function can also create plots of the PCA results and save them in the results directory.
#'
#' @param olink_data (tibble). A tibble with the data to be used in the PCA analysis
#' @param metadata (tibble). A tibble with metadata information to be used in the PCA plots. Default is NULL
#' @param variable (string). The name of the column in `metadata` that contains the variable
#' to be used to plot the points color. Default is "Disease"
#' @param palette (vector). A vector with the colors to be used in the PCA plots. Default is NULL
#' @param wide (logical). If TRUE, the data is assumed to be in wide format. Default is TRUE
#' @param impute (logical). If TRUE, the data is imputed before the PCA analysis. Default is TRUE
#' @param plots (logical). If TRUE, the function creates plots of the PCA results. Default is FALSE
#' @param save (logical). If TRUE, the plots are saved in the results directory. Default is FALSE
#'
#' @return (list). A list with the PCA results and, if requested, the PCA plots
#'   - pca_res (tibble). A tibble with the PCA results
#'   - loadings (tibble). A tibble with the PCA loadings
#'   - pca_plot (plot). A ggplot object with the PCA plot
#'   - loadings_plot (plot). A ggplot object with the PCA loadings plot
#' @export
#'
#' @examples
#' test_data <- example_data |> dplyr::select(DAid, Assay, NPX)
#' do_pca(test_data,
#'        metadata = example_metadata,
#'        wide = FALSE,
#'        variable = "Disease",
#'        plots = TRUE)
do_pca <- function(olink_data,
                   metadata = NULL,
                   variable = "Disease",
                   palette = NULL,
                   wide = T,
                   impute = T,
                   plots = F,
                   save = F) {

  wide_data <- widen_data(olink_data, wide)

  set.seed(123)
  if (isTRUE(impute)) {
    pca_rec <- recipes::recipe( ~ ., data = wide_data) |>
      recipes::update_role(DAid, new_role = "id")  |>
      recipes::step_normalize(recipes::all_predictors()) |>
      recipes::step_impute_knn(recipes::all_predictors(), neighbors = 5) |>
      recipes::step_pca(recipes::all_predictors())

    pca_prep <- recipes::prep(pca_rec)

    tidied_pca <- broom::tidy(pca_prep, 3)

  } else {
    pca_rec <- recipes::recipe( ~ ., data = wide_data) |>
      recipes::update_role(DAid, new_role = "id")  |>
      recipes::step_normalize(recipes::all_predictors()) |>
      recipes::step_pca(recipes::all_predictors())

    pca_prep <- recipes::prep(pca_rec)

    tidied_pca <- broom::tidy(pca_prep, 2)
  }

  loadings_data <- tidied_pca |>
    dplyr::rename(Assay = terms, Value = value, PC = component)

  pca_res <-  recipes::juice(pca_prep)

  if (isTRUE(plots)) {
    pca_plot <- plot_dim_reduction(pca_res, "PC1", "PC2", metadata, variable, palette)
    loadings_plot <- plot_loadings(tidied_pca)

    if (isTRUE(save)) {
      dir_name <- create_dir("results/pca_plots", date = T)
      ggplot2::ggsave(pca_plot, filename = paste0(dir_name, "/pca_plot.png"), width = 10, height = 8)
      ggplot2::ggsave(loadings_plot, filename = paste0(dir_name, "/loadings_plot.png"), width = 10, height = 8)
    }
    return(
      list(
        "pca_res" = pca_res,
        "loadings" = loadings_data,
        "pca_plot" = pca_plot,
        "loadings_plot" = loadings_plot
      )
    )
  } else {
    return(list("pca_res" = pca_res, "loadings" = loadings_data))
  }
}


#' Perform UMAP analysis
#'
#' This function performs a UMAP analysis on the data provided.
#' The function can also create plots of the UMAP results and save them in the results directory.
#'
#' @param olink_data (tibble). A tibble with the data to be used in the UMAP analysis
#' @param metadata (tibble). A tibble with metadata information to be used in the UMAP plots. Default is NULL
#' @param variable (string). The name of the column in `metadata` that contains the variable
#' to be used to plot the points color. Default is "Disease"
#' @param palette (vector). A vector with the colors to be used in the UMAP plots. Default is NULL
#' @param wide (logical). If TRUE, the data is assumed to be in wide format. Default is TRUE
#' @param impute (logical). If TRUE, the data is imputed before the UMAP analysis. Default is TRUE
#' @param plots (logical). If TRUE, the function creates plots of the UMAP results. Default is FALSE
#' @param save (logical). If TRUE, the plots are saved in the results directory. Default is FALSE
#'
#' @return (list). A list with the UMAP results and, if requested, the UMAP plots
#'   - umap_res (tibble). A tibble with the UMAP results
#'   - umap_plot (plot). A ggplot object with the UMAP plot
#' @export
#'
#' @examples
#' test_data <- example_data |> dplyr::select(DAid, Assay, NPX)
#' do_umap(test_data,
#'        metadata = example_metadata,
#'        wide = FALSE,
#'        variable = "Disease",
#'        plots = TRUE)
do_umap <- function(olink_data,
                    metadata = NULL,
                    variable = "Disease",
                    palette = NULL,
                    wide = T,
                    impute = T,
                    plots = F,
                    save = F) {

  wide_data <- widen_data(olink_data, wide)

  set.seed(123)
  if (isTRUE(impute)) {
    umap_rec <- recipes::recipe( ~ ., data = wide_data) |>
      recipes::update_role(DAid, new_role = "id")  |>
      recipes::step_normalize(recipes::all_predictors()) |>
      recipes::step_impute_knn(recipes::all_predictors(), neighbors = 5) |>
      embed::step_umap(recipes::all_predictors())

    umap_prep <- recipes::prep(umap_rec)

  } else {
    umap_rec <- recipes::recipe( ~ ., data = wide_data) |>
      recipes::update_role(DAid, new_role = "id")  |>
      recipes::step_normalize(recipes::all_predictors()) |>
      embed::step_umap(recipes::all_predictors())

    umap_prep <- recipes::prep(umap_rec)

  }

  umap_res <-  recipes::juice(umap_prep)

  if (isTRUE(plots)) {
    umap_plot <- plot_dim_reduction(umap_res, "UMAP1", "UMAP2", metadata, variable, palette)

    if (isTRUE(save)) {
      dir_name <- create_dir("results/umap_plots", date = T)
      ggplot2::ggsave(umap_plot, filename = paste0(dir_name, "/umap_plot.png"), width = 10, height = 8)
    }
    return(
      list(
        "umap_res" = umap_res,
        "umap_plot" = umap_plot
      )
    )
  } else {
    return("umap_res" = umap_res)
  }
}
