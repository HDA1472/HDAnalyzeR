# Load package functions with load_all() first
example_df <- import_df("data-raw/cancer_data_synthetic.rds")

unique_assays <- unique(example_df$Assay)[1:100]
example_data <- example_df |> dplyr::filter(Assay %in% unique_assays)
example_data <- example_data |>
  dplyr::mutate(
    Assay_Warning = "PASS",
    QC_Warning = "PASS"
  )

# Randomly change some QC_Warning values to "WARN" or "MANUAL_WARN"
set.seed(42)
num_rows <- nrow(example_data)
num_change <- round(num_rows * 0.01)  # Number of rows to change (1% of total rows)

indices_assay <- sample(num_rows, num_change)
example_data$Assay_Warning[indices_assay] <- "WARN"

indices_qc <- sample(num_rows, num_change)
example_data$QC_Warning[indices_qc] <- sample(c("WARN", "MANUAL_WARN"), num_change, replace = TRUE)

# Add DAid, LOD and PlateID columns
unique_samples_df <- example_data |>
  dplyr::distinct(Sample) |>
  dplyr::mutate(DAid = sprintf("DA%05d", dplyr::row_number()))

example_data <- example_data |>
  dplyr::left_join(unique_samples_df, by = "Sample") |>
  dplyr::select(DAid, everything())

# example_data <- example_data %>%
#   dplyr::mutate(LOD = seq(-11.0, 5.0, length.out = dplyr::n()))

example_data <- example_data |>
  dplyr::mutate(
    PlateID = dplyr::case_when(
      grepl("_1$", Sample) | grepl("_2$", Sample) | grepl("_3$", Sample) |
        grepl("_4$", Sample) | grepl("_5$", Sample) |
        grepl("_6$", Sample) | grepl("_7$", Sample) | grepl("_8$", Sample) |
        grepl("_9$", Sample) | grepl("_10$", Sample) ~ "Run001",
      grepl("_11$", Sample) | grepl("_12$", Sample) | grepl("_13$", Sample) |
        grepl("_14$", Sample) | grepl("_15$", Sample) |
        grepl("_16$", Sample) | grepl("_17$", Sample) | grepl("_18$", Sample) |
        grepl("_19$", Sample) | grepl("_20$", Sample) ~ "Run002",
      grepl("_21$", Sample) | grepl("_22$", Sample) | grepl("_23$", Sample) |
        grepl("_24$", Sample) | grepl("_25$", Sample) |
        grepl("_26$", Sample) | grepl("_27$", Sample) | grepl("_28$", Sample) |
        grepl("_29$", Sample) | grepl("_30$", Sample) ~ "Run003",
      grepl("_31$", Sample) | grepl("_32$", Sample) | grepl("_33$", Sample) |
        grepl("_34$", Sample) | grepl("_35$", Sample) |
        grepl("_36$", Sample) | grepl("_37$", Sample) | grepl("_38$", Sample) |
        grepl("_39$", Sample) | grepl("_40$", Sample) ~ "Run004",
      grepl("_41$", Sample) | grepl("_42$", Sample) | grepl("_43$", Sample) |
        grepl("_44$", Sample) | grepl("_45$", Sample) |
        grepl("_46$", Sample) | grepl("_47$", Sample) | grepl("_48$", Sample) |
        grepl("_49$", Sample) | grepl("_50$", Sample) ~ "Run005"
    )
  )
