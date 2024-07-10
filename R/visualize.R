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
#' create_protein_boxplot(join_data, c("A1BG", "A2M"), "AML", palette = "cancers12")
create_protein_boxplot <- function(join_data,
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
