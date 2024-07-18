utils::globalVariables(c("DAid", "Assay", "NPX", "adj.P.Val", "Age", "BMI", "n"))
#' Check the column types of the dataframe
#'
#' The function checks the column types of the input dataframe and returns the counts of each class.
#'
#' @param df (tibble). The input dataframe.
#'
#' @return class_summary (table). A table with the counts of each class in the dataframe.
#' @keywords internal
check_col_types <- function(df) {
  # Get the classes of all columns
  col_classes <- lapply(df, class)

  # Summarize the counts of each class
  class_summary <- table(unlist(col_classes))

  return(class_summary)
}

#' Calculate the percentage of NAs in each column
#'
#' The function calculates the percentage of NAs in each column of the input dataframe.
#' It filters out the columns with 0% missing data and returns the rest in descending order.
#'
#' @param df (tibble). The input dataframe.
#'
#' @return na_percentage (tibble). A tibble with the column names and the percentage of NAs in each column.
#' @keywords internal
calc_na_percentage_col <- function(df) {

  na_percentage <- df |>
    dplyr::summarise_all(~ round(sum(is.na(.)) / dplyr::n() * 100, 1)) |>
    tidyr::gather(key = "column", value = "na_percentage") |>
    dplyr::filter(na_percentage > 0) |>  # Filter out columns with no NAs
    dplyr::arrange(dplyr::desc(na_percentage))

  return(na_percentage)
}


#' Calculate the percentage of NAs in each row
#'
#' The function calculates the percentage of NAs in each row of the input dataframe.
#' It filters out the rows with 0% missing data and returns the rest in descending order.
#'
#' @param df (tibble). The input dataframe.
#'
#' @return na_percentage (tibble). A tibble with the DAids and the percentage of NAs in each row.
#' @keywords internal
calc_na_percentage_row <- function(df) {

  na_percentage <- df |>
    dplyr::rowwise() |>
    dplyr::mutate(na_percentage = round(sum(is.na(dplyr::across(dplyr::everything())))/ncol(df) * 100, 1)) |>
    dplyr::ungroup() |>
    dplyr::filter(na_percentage > 0) |>
    dplyr::arrange(dplyr::desc(na_percentage)) |>
    dplyr::select(DAid, na_percentage)

  return(na_percentage)
}


#' Check normality of the data
#'
#' @param df (tibble). The input dataframe.
#'
#' @return normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status.
#' @keywords internal
check_normality <- function(df) {

  future::plan(future::multicore)

  # Perform Shapiro-Wilk test in parallel
  p_values <- future.apply::future_lapply(df |> dplyr::select(-dplyr::any_of(c("DAid"))), function(column) {
    stats::shapiro.test(column)$p.value
  })

  adj.P.Val <- stats::p.adjust(unlist(p_values), method = "BH")

  # Determine normality based on adjusted p-values
  normality <- adj.P.Val > 0.05

  normality_results <- tibble::tibble(
    Protein = names(p_values),
    p_value = unlist(p_values),
    adj.P.Val = adj.P.Val,
    is_normal = normality
  ) |>
    dplyr::arrange((adj.P.Val))

  return(normality_results)
}


#' Print the summary of the quality control results
#'
#' The function prints the summary of the quality control results of the input dataframe.
#'
#' @param sample_n (numeric). The number of samples.
#' @param var_n (numeric). The number of variables.
#' @param class_summary (table). A table with the counts of each class in the dataframe.
#' @param na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column.
#' @param na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row.
#' @param normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status.
#' @param cor_results (tibble). A tibble with the filtered protein pairs and their correlation values.
#' @param heatmap (plot). A heatmap of protein-protein correlations.
#' @param threshold (numeric). The reporting protein-protein correlation threshold.
#'
#' @return NULL.
#' @keywords internal
print_summary <- function(sample_n, var_n, class_summary, na_percentage_col, na_percentage_row,
                          normality_results = F, cor_results = F, heatmap = F,  threshold = F) {

  print("Summary:")
  print("Note: In case of long output, only the first 10 rows are shown. To see the rest display the object with view()")
  print(paste0("Number of samples: ", sample_n))
  print(paste0("Number of variables: ", var_n))
  print("--------------------------------------")
  for (class_name in names(class_summary)) {
    print(paste(class_name, ":", class_summary[class_name]))
  }
  print("--------------------------------------")
  print("NA percentage in each column:")
  print(na_percentage_col)
  print("--------------------------------------")
  print("NA percentage in each row:")
  print(na_percentage_row)
  print("--------------------------------------")
  if (!isFALSE(normality_results)) {
    print("Normality test results:")
    print(normality_results)
    print("--------------------------------------")
  }
  if (!isFALSE(cor_results)) {
    print(paste0("Protein-protein correlations above ", threshold, ":"))
    print(cor_results)
    print("--------------------------------------")
  }
  if (!isFALSE(heatmap)) {
    print("Correlation heatmap:")
    print(heatmap)
    print("--------------------------------------")
  }

  invisible(NULL)
}


#' Create summary plots for age and sex
#'
#' The function creates three plots: two ridge plots for the age and BMI distributions and a bar plot for the number of samples per sex.
#'
#' @param metadata (tibble). The metadata dataframe.
#' @param disease_palette (character). The name of the palette to use for the disease variable. Default is NULL.
#' @param sex_palette (character). The name of the palette. Default is "sex_hpa".
#'
#' @return A list containing the following elements:
#'  - sex_barplot (plot). A bar plot for the number of samples per sex.
#'  - age_dist_plot (plot). A ridge plot for the age distributions.
#'  - bmi_dist_plot (plot). A ridge plot for the BMI distributions.
#' @keywords internal
plot_metadata_summary <- function(metadata, disease_palette = NULL, sex_palette = "sex_hpa") {

  plot_list <- list()

  # Ridge plot for age distributions
  if ("Age" %in% names(metadata)) {
    age_dist_plot <- metadata |>
      ggplot2::ggplot(ggplot2::aes(x = Age, y = Disease, fill = Disease)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
      ggplot2::labs(x = "Age", y = "Disease") +
      theme_hpa() +
      ggplot2::theme(legend.position = "none")

    if (is.null(names(disease_palette)) && !is.null(disease_palette)) {
      age_dist_plot <- age_dist_plot + scale_fill_hpa(disease_palette)
    } else if (!is.null(disease_palette)) {
      age_dist_plot <- age_dist_plot + ggplot2::scale_fill_manual(values = disease_palette)
    }

    plot_list$age_dist_plot <- age_dist_plot
  }

  # Ridge plot for BMI distributions
  if ("BMI" %in% names(metadata)) {
    bmi_dist_plot <- metadata |>
      ggplot2::ggplot(ggplot2::aes(x = BMI, y = Disease, fill = Disease)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
      ggplot2::labs(x = "BMI", y = "Disease") +
      theme_hpa() +
      ggplot2::theme(legend.position = "none")

    if (is.null(names(disease_palette)) && !is.null(disease_palette)) {
      bmi_dist_plot <- bmi_dist_plot + scale_fill_hpa(disease_palette)
    } else if (!is.null(disease_palette)) {
      bmi_dist_plot <- bmi_dist_plot + ggplot2::scale_fill_manual(values = disease_palette)
    }

    plot_list$bmi_dist_plot <- bmi_dist_plot
  }

  # Bar plot for the number of samples
  if ("Sex" %in% names(metadata)) {
    sex_barplot <- metadata |>
      dplyr::count(Disease, Sex) |>
      ggplot2::ggplot(ggplot2::aes(x = n, y = Disease, fill = Sex)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::labs(x = "Number of samples", y = "") +
      theme_hpa()

    if (is.null(names(sex_palette)) && !is.null(sex_palette)) {
      sex_barplot <- sex_barplot + scale_fill_hpa(sex_palette)
    } else if (!is.null(sex_palette)) {
      sex_barplot <- sex_barplot + ggplot2::scale_fill_manual(values = sex_palette)
    }

    plot_list$sex_barplot <- sex_barplot
  }
  return(plot_list)
}


#' Summarize the quality control results of Olink data
#'
#' The function summarizes the quality control results of the input dataframe.
#' It can handles both long and wide dataframes.
#'
#' @param df (tibble). The input dataframe.
#' @param wide (logical). Whether the input dataframe is in wide format. Default is TRUE.
#' @param threshold (numeric). The reporting protein-protein correlation threshold. Default is 0.8.
#' @param report (logical). Whether to print the summary. Default is TRUE.
#'
#' @return A list containing the following elements:
#'   - na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column.
#'   - na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row.
#'   - normality_results (tibble). A tibble with the protein names, p-values, adjusted p-values, and normality status.
#'   - cor_matrix (matrix). A matrix of protein-protein correlations.
#'   - cor_results (tibble). A tibble with the filtered protein pairs and their correlation values.
#'   - heatmap (ggplot). A heatmap of protein-protein correlations.
#' @export
#'
#' @examples
#' qc_summary_data(example_data, wide = FALSE, threshold = 0.7)
qc_summary_data <- function(df, wide = T, threshold = 0.8, report = T) {

  if (isFALSE(wide)) {
    wide_data <- widen_data(df)
  } else {
    wide_data <- df
  }
  sample_n <- nrow(wide_data)
  protein_n <- ncol(wide_data) - 1
  class_summary <- check_col_types(wide_data)
  na_percentage_col <- calc_na_percentage_col(wide_data)
  na_percentage_row <- calc_na_percentage_row(wide_data)
  normality_results <- check_normality(wide_data)
  cor <- create_corr_heatmap(wide_data |> dplyr::select(-dplyr::any_of(c("DAid"))),
                             threshold = threshold)
  cor_matrix <- cor$cor_matrix
  cor_results <- cor$cor_results
  p <- cor$p

  if (isTRUE(report)) {
    print_summary(sample_n, protein_n, class_summary, na_percentage_col, na_percentage_row,
                  normality_results, cor_results, p, threshold)
  }

  return(
    list(
      "na_percentage_col" = na_percentage_col,
      "na_percentage_row" = na_percentage_row,
      "normality_results" = normality_results,
      "cor_matrix" = cor_matrix,
      "cor_results" = cor_results,
      "heatmap" = p
      )
    )
}


#' Summarize the quality control results of metadata
#'
#' The function summarizes the quality control results of the input dataframe.
#'
#' @param metadata (tibble). The metadata dataframe.
#' @param disease_palette (character). The name of the palette to use for the disease variable. Default is NULL.
#' @param sex_palette (character). The name of the palette. Default is "sex_hpa".
#' @param report (logical). Whether to print the summary. Default is TRUE.
#'
#' @return A list containing the following elements:
#'   - na_percentage_col (tibble). A tibble with the column names and the percentage of NAs in each column.
#'   - na_percentage_row (tibble). A tibble with the DAids and the percentage of NAs in each row.
#'   - sex_barplot (plot). A bar plot for the number of samples per sex.
#'   - age_dist_plot (plot). A ridge plot for the age distributions.
#'   - bmi_dist_plot (plot). A ridge plot for the BMI distributions.
#' @export
#'
#' @examples
#' qc_summary_metadata(example_metadata)
qc_summary_metadata <- function(metadata, disease_palette = NULL, sex_palette = "sex_hpa", report = T) {

  sample_n <- nrow(metadata)
  var_n <- ncol(metadata) - 1
  class_summary <- check_col_types(metadata)
  na_percentage_col <- calc_na_percentage_col(metadata)
  na_percentage_row <- calc_na_percentage_row(metadata)

  if (isTRUE(report)) {
    print_summary(sample_n, var_n, class_summary, na_percentage_col, na_percentage_row)
  }

  metadata_plot <- plot_metadata_summary(metadata, disease_palette, sex_palette)

  return(
    list(
      "na_percentage_col" = na_percentage_col,
      "na_percentage_row" = na_percentage_row,
      "sex_barplot" = metadata_plot$sex_barplot,
      "age_dist_plot" = metadata_plot$age_dist_plot,
      "bmi_dist_plot" = metadata_plot$bmi_dist_plot
    )
  )
}
