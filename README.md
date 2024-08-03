# HDAnalyzeR <img src="man/figures/logo.png" align="right" height="200"/>

[![R-CMD-check](https://github.com/HDA1472/DA_RPackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/HDA1472/DA_RPackage/actions/workflows/R-CMD-check.yaml) [![Version](https://img.shields.io/badge/Version-0.1.0-blue)](https://github.com/HDA1472/DA_RPackage) [![License](https://img.shields.io/badge/license-Apache2.0-yellow)](https://github.com/HDA1472/DA_RPackage/blob/main/LICENSE.md)

HDAnalyzeR is an R package developed by the Human Disease Blood Atlas project, designed to facilitate proteomics analysis for biomarker selection from blood plasma samples. It is optimized to work with Olink proteomics data, but it can be adapted to other proteomics platforms. In order to use the package without issues the data should have these three necessary columns: `DAid` with the Sample IDs, `Assay` with the protein names, and `NPX` with the protein expression data. The metadata should contain the `DAid`, `Disease`, and `Sex` columns, where the `Disease` column should contain the different class names (Healthy, Disease, etc.), while in the `Sex` column the data should be encoded as M (males) and F (females).

HDAnalyzeR offers ready-to-use functions for common proteomics tasks such as protein differential expression analysis, classification models, imputation methods, dimensionality reduction, and data visualization, aiming to streamline workflows and enhance the standardization and efficiency of biomarker discovery in disease research.

## Installation

> ⚠️ We recommend to install it in a toy environment as it is currently under development.

You can install the latest version of HDAnalyzeR from GitHub:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install HDAnalyzeR
options(timeout = 1200)  # Set timeout to 20 minutes to avoid timeout errors
devtools::install_github("HDA1472/DA_RPackage")
```

## Example

The following example showcases how to perform a differential expression analysis. It is one of the many features of HDAnalyzeR. A complete guide is available through [package's documentation](https://hda1472.github.io/DA_RPackage/).

``` r
library(HDAnalyzeR)

# Prepare data
wide_data <- widen_data(example_data)

# Run differential expression analysis
de_results <- do_limma(wide_data, example_metadata, case = "AML", control = c("CLL", "MYEL"))

# DE results and volcano plot for AML
de_results$de_results
de_results$volcano_plot
```

## Issues and Support

If you encounter any bugs or you want to recommend new features and changes to existing ones, please open a new issue on our GitHub repository.

## Contact

For any questions or further information, please contact us at [konstantinos.antonopoulos\@scilifelab.se](mailto:konstantinos.antonopoulos@scilifelab.se){.email}.
