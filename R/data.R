#' Example cancer data
#'
#' A subset of data from a synthetic cancer dataset. Only the first 100 Assays.
#' With DAid, Assay_Warning, QC_Warning, and PlateID as extra columns.
#' The original dataset was processed with the `process_example_data` script.
#'
#' @format ## `example_data`
#' A data frame with 56,142 rows and 10 columns:
#' \describe{
#'   \item{DAid}{The Disease Atlas sample ID}
#'   \item{Sample}{The Sample ID}
#'   \item{OlinkID}{The Olink ID of the Assay}
#'   \item{UniProt}{The UniProt ID of the Assay}
#'   \item{Assay}{The Assay name}
#'   \item{Panel}{The Olink Panel in which the Assay belongs}
#'   \item{NPX}{The NPX value for each Sample and Assay}
#'   \item{Assay_Warning}{The Assay warning status for the sample}
#'   \item{QC_Warning}{The QC warning status for the sample}
#'   \item{PlateID}{The ID of the plate where the sample was processed}
#' }
#' @source <https://github.com/buenoalvezm/Pan-cancer-profiling/blob/main/data/cancer_data_synthetic.rds>
"example_data"
