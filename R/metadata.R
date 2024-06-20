#' Example cancer metadata
#'
#' A subset of metadata from a synthetic cancer dataset. With DAid, Age, BMI, and Cohort  as extra
#' columns. The original dataset was processed with the `process_example_metadata` script.
#'
#' @format ## `example_metadata`
#' A data frame with 586 rows and 9 columns:
#' \describe{
#'   \item{DAid}{The Disease Atlas sample ID}
#'   \item{Sample}{The Sample ID}
#'   \item{GROUP}{The cancer type}
#'   \item{Stage}{The cancer stage}
#'   \item{Grade}{The cancer grade}
#'   \item{Sex}{The patient sex}
#'   \item{Age}{The patient age}
#'   \item{BMI}{The patient BMI}
#'   \item{Cohort}{The cohort the patient belongs to}
#' }
#' @source <https://github.com/buenoalvezm/Pan-cancer-profiling/blob/main/data/cancer_metadata_synthetic.rds>
"example_metadata"
