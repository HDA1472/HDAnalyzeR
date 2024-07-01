#' Remove batch effects
#'
#' This function removes batch effects from the data using the limma package.
#' It converts the dataframe into matrix and transposes it to get it ready for limma.
#' It removes the batch effects and then converts the data back to normal.
#'
#' @param wide_data (tibble). A dataframe containing the data to be normalized. The data should be in wide format.
#' @param batch (vector). A character vector containing the batch information.
#' @param batch2 (vector). A character vector containing the second batch information. Default is NULL.
#'
#' @return no_batch_effects (tibble). A dataframe containing the data without batch effects.
#' @export
#'
#' @examples
#' test_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX)
#'
#' batch <- test_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Cohort), by = "DAid") |>
#'   dplyr::pull("Cohort")
#'
#' input_data <- test_data |> dplyr::select(-DAid)
#' remove_batch_effects(input_data, batch)
remove_batch_effects <- function(wide_data,
                                 batch,
                                 batch2 = NULL
                                 ) {

  # Prepare the data for limma
  mat_data <- as.matrix(wide_data)
  transposed_mat <- t(mat_data)

  # Remove batch effects
  no_batch_effects_res <- limma::removeBatchEffect(transposed_mat, batch=batch, batch2=batch2)
  transposed_no_batch_effects <- t(no_batch_effects_res)
  no_batch_effects <- tibble::as_tibble(transposed_no_batch_effects)

  return(no_batch_effects)
}


#' Normalize data
#'
#' This function normalizes the data by removing batch effects and scaling the data.
#' It first converts the data to wide format if it is not already in wide format.
#' It then removes the batch effects and scales or centers the data.
#' To remove batch effects, it uses the remove_batch_effects function, that utilizes limma package.
#' For scaling, it uses the scale function from base R.
#'
#' @param olink_data (tibble). A dataframe containing Olink data to be normalized.
#' @param metadata (tibble). A dataframe containing the metadata information.
#' @param wide (logical). A logical value indicating whether the data is in wide format. Default is TRUE.
#' @param center (logical). A logical value indicating whether to center the data. Default is TRUE.
#' @param scale (logical). A logical value indicating whether to scale the data. Default is TRUE.
#' @param batch (character). The metadata column containing the batch information. Default is NULL.
#' @param batch2 (character). The metadata column containing the second batch information. Default is NULL.
#' @param return_long (logical). A logical value indicating whether to return the data in long format. Default is FALSE.
#' @param save (logical). A logical value indicating whether to save the data. Default is FALSE.
#' @param file_name (character). The name of the file to be saved. Default is "normalized_data".
#'
#' @return scaled_data (tibble). A dataframe containing the normalized data.
#' @export
#'
#' @examples
#' scaled_data <- normalize_data(example_data, example_metadata, wide = FALSE, batch = "Cohort")
normalize_data <- function(olink_data,
                           metadata = NULL,
                           wide = T,
                           center = T,
                           scale = T,
                           batch = NULL,
                           batch2 = NULL,
                           return_long = F,
                           save = F,
                           file_name = "normalized_data"
                           ) {

  if (isFALSE(wide)) {
    wide_data <- olink_data |>
      dplyr::select(DAid, Assay, NPX) |>
      tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  } else {
    wide_data <- olink_data
  }

  # Prepare the data for scaling
  id_col <- wide_data |> dplyr::pull(DAid)
  input_data <- wide_data |> dplyr::select(-DAid)

  # Remove batch effects
  if (!is.null(batch)) {
    batch <- wide_data |>
      dplyr::left_join(metadata |>
                         dplyr::select(dplyr::any_of(c("DAid", batch))),
                       by = "DAid") |>
      dplyr::pull(batch)
    if (!is.null(batch2)) {
      batch2 <- wide_data |>
        dplyr::left_join(metadata |>
                           dplyr::select(DAid, batch2),
                         by = "DAid") |>
        dplyr::pull(batch2)
    }
    data_wo_batch_effects <- remove_batch_effects(input_data, batch = batch, batch2 = batch2)
  } else {
    data_wo_batch_effects <- input_data
  }

  # Scale the data
  scaled_data <- as.data.frame(scale(data_wo_batch_effects, center = center, scale = scale))
  names(scaled_data) <- names(input_data)

  # Prepare data to be returned
  scaled_data$DAid <- id_col
  scaled_data <- scaled_data |> dplyr::select(DAid, dplyr::everything())

  if (isTRUE(return_long)) {
    scaled_data <- scaled_data |> tidyr::pivot_longer(cols = -DAid, names_to = "Assay", values_to = "NPX")
  }

  if (isTRUE(save)) {
    dir_name <- create_dir("data/processed/data_metadata", date = T)
    save_df(scaled_data, dir_name, file_name, "rda")
  }

  return(scaled_data)
}

