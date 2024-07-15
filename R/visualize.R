utils::globalVariables(c("Value"))
#' Create boxplots for proteins
#'
#' This function creates a boxplot for the top n proteins in the dataset.
#' It annotates the boxplot with color for the selected disease.
#' It is also possible to add points to the boxplot.
#'
#' @param join_data (tibble). The dataset with the wide Olink data joined with the metadata.
#' @param proteins (vector). The proteins to include in the boxplot.
#' @param disease (character). The disease to annotate.
#' @param points (logical). Whether to add points to the boxplot.
#' @param palette (character). The color palette to use. Default is red3.
#'
#' @return boxplot_panel (plot). The boxplot panel with the selected proteins.
#' @export
#'
#' @examples
#' wide_data <- widen_data(example_data, FALSE)
#' join_data <- wide_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Disease, Sex))
#' plot_protein_boxplot(join_data, c("AARSD1", "ABL1"), "AML", palette = "cancers12")
plot_protein_boxplot <- function(join_data,
                                 proteins,
                                 disease,
                                 points = T,
                                 palette = NULL) {

  # Prepare palettes
  pals <- get_hpa_palettes()
  if (!is.null(palette) && is.null(names(palette))) {
    pal <- pals[palette]
    pal <- unlist(pals[[palette]])
  } else if (!is.null(palette)) {
    pal <- palette
  } else {
    pal <- "red3"
  }

  long_data <- join_data |>
    dplyr::select(Disease, dplyr::all_of(proteins)) |>
    tidyr::pivot_longer(cols = !Disease, names_to = "Protein", values_to = "Value")

  long_data$Protein <- factor(long_data$Protein, levels = proteins, labels = proteins)

  # Create boxplot
  boxplot <- long_data |>
    ggplot2::ggplot(ggplot2::aes(x = Disease, y = Value)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_boxplot(data = subset(long_data, Disease == disease),
                          ggplot2::aes(fill = Disease),
                          alpha = 0.5,
                          show.legend = FALSE,
                          outlier.shape = NA)

  if (isTRUE(points)) {
    boxplot <- boxplot +
      ggplot2::geom_point(data = subset(long_data, Disease != disease),
                          position = ggplot2::position_jitter(width = 0.1),
                          color = 'grey',
                          alpha = 0.3)

    if (!is.null(palette)) {
      boxplot <- boxplot +
        ggplot2::geom_point(data = subset(long_data, Disease == disease),
                            ggplot2::aes(fill = Disease),
                            position = ggplot2::position_jitter(width = 0.1),
                            color = pal[disease],
                            alpha = 0.5,
                            show.legend = FALSE)
    } else {
      boxplot <- boxplot +
        ggplot2::geom_point(data = subset(long_data, Disease == disease),
                            ggplot2::aes(fill = Disease),
                            position = ggplot2::position_jitter(width = 0.1),
                            color = pal,
                            alpha = 0.5,
                            show.legend = FALSE)
    }
  }

  boxplot_panel <- boxplot +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::xlab('') +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~ Protein, scale="free_y")

  return(boxplot_panel)
}


#' Create a scatter plot with regression line
#'
#' This function creates a scatter plot with a linear regression line.
#' It is possible to add the standard error of the regression line, as well as the R-squared and p-value.
#'
#' @param plot_data (tibble). The wide dataset containing the data to plot as cols.
#' @param x (character). The column name of the x-axis variable.
#' @param y (character). The column name of the y-axis variable.
#' @param se (logical). Whether to add the standard error of the regression line.
#' @param line_color (character). The color of the regression line.
#' @param r_2 (logical). Whether to add the R-squared and p-value to the plot.
#'
#' @return scatter (plot). The scatter plot with the regression line.
#' @export
#'
#' @examples
#' plot_scatter_with_regression(example_metadata, "Age", "BMI")
plot_scatter_with_regression <- function(plot_data,
                                         x,
                                         y,
                                         se = F,
                                         line_color = "black",
                                         r_2 = T) {
  # Fit the linear model
  formula <- stats::as.formula(paste(y, "~", x))
  model <- stats::lm(formula, data = plot_data)

  # Get the R-squared and p-value
  summary_model <- summary(model)
  r_squared <- summary_model$r.squared
  p_val <- stats::coef(summary_model)[2, 4]

  # Create the plot
  x <- rlang::sym(x)
  y <- rlang::sym(y)
  scatter <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = !!x, y = !!y)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = se, color = line_color)

  if (isTRUE(r_2)) {
    scatter <- scatter +
      ggplot2::annotate("text",
                        x = Inf,
                        y = Inf,
                        label = paste("R2 =", round(r_squared, 2), "\nP =", format.pval(round(p_val, 4))),
                        hjust = 1.1,
                        vjust = 2,
                        size = 5,
                        color = "black") +
    ggplot2::theme_classic()
  } else {
    scatter <- scatter +
      ggplot2::theme_classic()
  }

  return(scatter)
}
