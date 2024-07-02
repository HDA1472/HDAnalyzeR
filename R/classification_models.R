#' Split data into training and test sets
#'
#' This function splits the data into training and test sets based on user defined ratio.
#'
#' @param join_data (tibble). Wide data joined with metadata.
#' @param ratio (numeric). Ratio of training data to test data. Default is 0.75.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with two elements:
#'  - train_set (tibble). The training set.
#'  - test_set (tibble). The test set.
#' @export
#'
#' @examples
#' wide_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX)
#' join_data <- wide_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Disease, Sex))
#'
#' data_splits <- split_data(join_data)
split_data <- function(join_data, ratio = 0.75, seed = 123) {

  set.seed(seed)
  data_split <- rsample::initial_split(join_data, prop = ratio, strata = "Disease")
  train_data <- rsample::training(data_split)
  test_data <- rsample::testing(data_split)

  return(
    list(
      "train_set" = train_data,
      "test_set" = test_data
    )
  )
}


#' Filter data in case of sex specific diseases
#'
#' This function filters the control data based on the disease.
#' It also filters the diseases vector to keep only the diseases that we will pick samples from.
#'
#' @param control_data (tibble). Control data to be filtered.
#' @param disease (character). Disease to be filtered.
#' @param diseases (vector). Vector of diseases.
#' @param only_female (vector). Vector of diseases that are female specific.
#' @param only_male (vector). Vector of diseases that are male specific.
#'
#' @return A list with two elements:
#'  - control_data (tibble). Filtered control data.
#'  - diseases_subset (vector). Filtered diseases vector.
#' @export
#'
#' @examples
#' wide_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX)
#' join_data <- wide_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Disease, Sex))
#'
#' diseases <- unique(example_metadata$Disease)
#' control_data <- join_data |> dplyr::filter(Disease != "BRC")
#'
#' filter_sex_specific_disease(control_data,
#'                             "BRC",
#'                             diseases,
#'                             only_female = c("BRC", "CVX", "ENDC", "OVC"),
#'                             only_male = "PRC")
filter_sex_specific_disease <- function(control_data,
                                        disease,
                                        diseases,
                                        only_female = NULL,
                                        only_male = NULL
                                        ) {

  if(!is.null(only_female) & disease %in% only_female) {
    control_data <- control_data |>
      dplyr::filter(Sex == "Female")
    diseases_subset <- diseases[!diseases %in% only_male]
  } else if(!is.null(only_male) & disease %in% only_male) {
    control_data <- control_data |>
      dplyr::filter(Sex == "Male")
    diseases_subset <- diseases[!diseases %in% only_female]
  } else {
    control_data <- control_data
    diseases_subset <- diseases
  }

  return(
    list(
      "control_data" = control_data,
      "diseases_subset" = diseases_subset
    )
  )
}


#' Make groups for classification models
#'
#' This function creates class-balanced groups for classification models.
#' It separates the data into control and case groups.
#' It also filters the control data based on the disease in case of sex specific diseases.
#'
#' @param join_data (tibble). Wide data joined with metadata.
#' @param diseases (vector). Vector of diseases.
#' @param only_female (vector). Vector of diseases that are female specific.
#' @param only_male (vector). Vector of diseases that are male specific.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with combined, balanced control-case groups for each disease
#' @export
#'
#' @examples
#' wide_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX)
#' join_data <- wide_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Disease, Sex))
#'
#' diseases <- unique(example_metadata$Disease)
#'
#' group_list <- make_groups(join_data,
#'                           diseases,
#'                           only_female = c("BRC", "CVX", "ENDC", "OVC"),
#'                           only_male = "PRC")
make_groups <- function(join_data,
                        diseases,
                        only_female = NULL,
                        only_male = NULL,
                        seed = 123
                        ) {

  set.seed(seed)
  group_list <- list()

  # Separate data in control and case groups
  for (disease in diseases) {
    case_data <- join_data |> dplyr::filter(Disease == disease)
    control_data <- join_data |> dplyr::filter(Disease != disease)

    # Filter for gender if disease is gender specific
    filtered_results <- filter_sex_specific_disease(control_data,
                                                    disease,
                                                    diseases,
                                                    only_female,
                                                    only_male)

    control_data <- filtered_results$control_data
    diseases_subset <- filtered_results$diseases_subset  # Keep only diseases that we will pick samples from

    # Calculate amount of data from each disease in control group
    case_sample_num <- nrow(case_data)
    samples_per_disease <- ceiling(case_sample_num/(length(unique(diseases_subset))-1))

    # Loop over control group, randomly select the number of each control class samples
    group <- case_data
    for (control_class in diseases_subset[diseases_subset != disease]) {
      control_class_data <- control_data |>
        dplyr::filter(Disease == control_class) |>
        dplyr::sample_n(size = samples_per_disease, replace = TRUE)
      group <- rbind(group, control_class_data)
    }

    group_list[[disease]] <- group
    print(disease)
    print(case_sample_num)
    print(samples_per_disease)
  }
  return(group_list)
}


