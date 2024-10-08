% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{normalize_data}
\alias{normalize_data}
\title{Normalize data and remove batch effects}
\usage{
normalize_data(
  olink_data,
  metadata = NULL,
  wide = TRUE,
  center = TRUE,
  scale = TRUE,
  batch = NULL,
  batch2 = NULL,
  return_long = FALSE,
  save = FALSE,
  file_name = "normalized_data"
)
}
\arguments{
\item{olink_data}{A dataset containing Olink data to be normalized.}

\item{metadata}{A dataset containing the metadata information.}

\item{wide}{A logical value indicating whether the data is in wide format. Default is TRUE.}

\item{center}{A logical value indicating whether to center the data. Default is TRUE.}

\item{scale}{A logical value indicating whether to scale the data. Default is TRUE.}

\item{batch}{The metadata column containing the batch information. In order to correct for batch effects, this parameter should be provided. Default is NULL.}

\item{batch2}{The metadata column containing the second batch information. Default is NULL.}

\item{return_long}{A logical value indicating whether to return the data in long format. Default is FALSE.}

\item{save}{A logical value indicating whether to save the data. Default is FALSE.}

\item{file_name}{The name of the file to be saved. Default is "normalized_data".}
}
\value{
A tibble containing the normalized data.
}
\description{
\code{normalize_data()} normalizes the data by scaling them and removing their batch effects.
It first converts the data to wide format if they are not already. It then removes
the batch effects and scales or centers the data. To remove batch effects, it uses the
\code{remove_batch_effects()}, that utilizes limma package. For scaling, it uses the \code{scale()}
from base R.
}
\examples{
# Non-normalized data
example_data |>
  dplyr::select(DAid, Assay, NPX) |>
  tidyr::pivot_wider(names_from = "Assay", values_from = "NPX")

# Center data
normalize_data(example_data, example_metadata, wide = FALSE, center = TRUE, scale = FALSE)

# Center and scale data (z-score scaling)
normalize_data(example_data, example_metadata, wide = FALSE, center = TRUE, scale = TRUE)

# Center, scale and remove batch effects
normalize_data(example_data, example_metadata, wide = FALSE, batch = "Cohort")
}
