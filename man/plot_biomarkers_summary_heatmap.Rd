% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{plot_biomarkers_summary_heatmap}
\alias{plot_biomarkers_summary_heatmap}
\title{Plot a summary heatmap of the combined differential expression and classification models results}
\usage{
plot_biomarkers_summary_heatmap(
  de_results,
  ml_results,
  order_by = NULL,
  pval_lim = 0.05,
  logfc_lim = 0
)
}
\arguments{
\item{de_results}{A list of differential expression results.}

\item{ml_results}{A list of classification models results.}

\item{order_by}{The control group to order the heatmap.}

\item{pval_lim}{The p-value limit to filter the differential expression results of the \code{order_by} group.}

\item{logfc_lim}{The log2 fold change limit to filter the differential expression results  of the \code{order_by} group.}
}
\value{
The summary heatmap of the combined differential expression and classification models results.
}
\description{
\code{plot_biomarkers_summary_heatmap} plots a summary heatmap of the combined differential
expression and classification models results. The heatmap shows the log2 fold change
and adjusted p-value of the differential expression results, and the scaled importance
and sign of the classification models results. The heatmap is ordered and the selected
assays are based on the specified control group.
}
\details{
It is very important the de_results and ml_results are in the same order
and in the right format (see examples).
}
\examples{
# Prepare differential expression results
de_results_amlbrc <- do_limma(example_data,
                              example_metadata,
                              case = "AML",
                              control = c("BRC"),
                              wide = FALSE,
                              only_female = "BRC")

de_results_amlcll <- do_limma(example_data,
                              example_metadata,
                              case = "AML",
                              control = c("CLL"),
                              wide = FALSE,
                              only_female = "BRC")

de_results_amllungc <- do_limma(example_data,
                                example_metadata,
                                case = "AML",
                                control = c("LUNGC"),
                                wide = FALSE,
                                only_female = "BRC")

# Combine the results
res_de <- list("BRC" = de_results_amlbrc,
               "CLL" = de_results_amlcll,
               "LUNGC" = de_results_amllungc)

res_amlbrc <- do_rreg(example_data,
                      example_metadata,
                      case = "AML",
                      type = 'elnet',
                      control = c("BRC"),
                      wide = FALSE,
                      only_female = "BRC",
                      cv_sets = 2,
                      grid_size = 1,
                      ncores = 1)

res_amlcll <- do_rreg(example_data,
                      example_metadata,
                      case = "AML",
                      type = 'elnet',
                      control = c("CLL"),
                      wide = FALSE,
                      only_female = "BRC",
                      cv_sets = 2,
                      grid_size = 1,
                      ncores = 1)

res_amllungc <- do_rreg(example_data,
                        example_metadata,
                        case = "AML",
                        control = c("LUNGC"),
                        type = 'elnet',
                        wide = FALSE,
                        only_female = "BRC",
                        cv_sets = 2,
                        grid_size = 1,
                        ncores = 1)

# Combine the results
res_ml <- list("BRC" = res_amlbrc,
               "CLL" = res_amlcll,
               "LUNGC" = res_amllungc)

# Create the summary heatmap
plot_biomarkers_summary_heatmap(res_de, res_ml, order_by = "BRC")
}
