create_protein_boxplot <- function(join_data,
                                   features,
                                   disease,
                                   nfeatures = 9,
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

  top_features <- features |>
    dplyr::arrange(desc(Scaled_Importance)) |>
    dplyr::select(Variable) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    head(nfeatures)
  top_features <- top_features[['Variable']]

  long_data <- join_data |>
    dplyr::select(Disease, dplyr::all_of(top_features)) |>
    tidyr::pivot_longer(cols = !Disease, names_to = "Protein", values_to = "Value")

  long_data$Protein <- factor(long_data$Protein, levels = top_features, labels = top_features)

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
