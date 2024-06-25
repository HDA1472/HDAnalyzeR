# Load package functions with load_all() first
example_df <- import_df("data-raw/cancer_metadata_synthetic.rds")

# Add DAid, Age, BMI, and Cohort columns
unique_samples_df <- example_df |>
  dplyr::distinct(Sample) |>
  dplyr::mutate(DAid = sprintf("DA%05d", dplyr::row_number()))

example_metadata <- example_df |>
  dplyr::left_join(unique_samples_df, by = "Sample") |>
  dplyr::select(DAid, everything()) |>
  dplyr::mutate(
    Age = round(runif(dplyr::n(), min = 40, max = 90)),
    BMI = round(runif(dplyr::n(), min = 20, max = 35), 1),
    Cohort = dplyr::case_when(
      GROUP %in% c("BRC", "CVX", "ENDC", "OVC", "PRC") ~ "Gender_specific",
      TRUE ~ "UCAN"
    )
  ) |>
  dplyr::rename(Disease = GROUP)
