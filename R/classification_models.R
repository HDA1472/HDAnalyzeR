utils::globalVariables(c("roc_auc"))
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
  }

  return(group_list)
}


#' Hyperparameter optimization for elastic net
#'
#' This function performs hyperparameter optimization for elastic net models.
#' It uses the glmnet engine for logistic regression and tunes either only penalty or both penalty and mixture.
#'
#'
#' @param train_data (list). List of training data sets from make_groups().
#' @param exclude_cols (vector). Columns to exclude from the model.
#' @param disease (character). Disease to predict.
#' @param type (character). Type of regularization. Default is "lasso". Other options are "ridge" and "elnet".
#' @param metric (function). Metric to optimize. Default is roc_auc.
#' @param cv_sets (numeric). Number of cross-validation sets. Default is 5.
#' @param grid_size (numeric). Size of the grid for hyperparameter optimization. Default is 10.
#' @param ncores (numeric). Number of cores to use for parallel processing. Default is 4.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with three elements:
#'  - elnet_tune (tibble). Hyperparameter optimization results.
#'  - train_set (tibble). Training set.
#'  - wf (workflow). Workflow object.
#' @export
#'
#' @examples
#' wide_data <- example_data |>
#'   dplyr::select(DAid, Assay, NPX) |>
#'   tidyr::pivot_wider(names_from = Assay, values_from = NPX)
#' join_data <- wide_data |>
#'   dplyr::left_join(example_metadata |> dplyr::select(DAid, Disease, Sex))
#' diseases <- unique(example_metadata$Disease)
#' group_list <- make_groups(join_data,
#'                           diseases,
#'                           only_female = c("BRC", "CVX", "ENDC", "OVC"),
#'                           only_male = "PRC")
#'
#' hypopt_res <- elnet_hypopt(group_list, "Sex", "AML", type = "elnet", grid_size = 5, ncores = 1)
elnet_hypopt <- function(train_data,
                         exclude_cols,
                         disease,
                         type = "lasso",
                         metric = roc_auc,
                         cv_sets = 5,
                         grid_size = 10,
                         ncores = 4,
                         seed = 123
) {

  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  }

  # Prepare train data and create cross-validation sets with binary classifier
  train_set <- train_data[[disease]] |>
    dplyr::mutate(Disease = ifelse(Disease == disease, 1, 0)) |>
    dplyr::mutate(Disease = as.factor(Disease)) |>
    dplyr::select(-dplyr::any_of(exclude_cols))
  train_folds <- rsample::vfold_cv(train_set, v = cv_sets, strata = Disease)

  elnet_rec <- recipes::recipe(Disease ~ ., data = train_set) |>
    recipes::update_role(DAid, new_role = "id") |>
    recipes::step_normalize(recipes::all_numeric()) |>
    recipes::step_nzv(recipes::all_numeric()) |>
    recipes::step_corr(recipes::all_numeric()) |>
    recipes::step_impute_knn(recipes::all_numeric())

  mixture_value <- switch(type,
                          "lasso" = 1,
                          "ridge" = 0,
                          "elnet" = NULL,  # allow tuning for elastic net
                          stop("Invalid type specified. Choose 'lasso', 'ridge', or 'elnet'."))

  if (is.null(mixture_value)) {
    elnet_spec <- parsnip::logistic_reg(
      penalty = tune::tune(),  # lambda
      mixture = tune::tune()  # alpha
    ) |>
      parsnip::set_engine("glmnet")
  } else {
    elnet_spec <- parsnip::logistic_reg(
      penalty = tune::tune(),
      mixture = mixture_value
    ) |>
      parsnip::set_engine("glmnet")
  }

  elnet_wf <- workflows::workflow() |>
    workflows::add_model(elnet_spec) |>
    workflows::add_recipe(elnet_rec)

  elnet_grid <- elnet_wf |>
    workflows::extract_parameter_set_dials() |>
    dials::grid_latin_hypercube(size = grid_size)

  roc_res <- yardstick::metric_set(yardstick::roc_auc)

  set.seed(seed)
  ctrl <- tune::control_grid(save_pred = TRUE, parallel_over = "everything")
  elnet_tune <- elnet_wf |>
    tune::tune_grid(
      train_folds,
      grid = elnet_grid,
      control = ctrl,
      metrics = roc_res
    )

  return(list(
    "elnet_tune" = elnet_tune,
    "train_set" = train_set,
    "wf" = elnet_wf
  ))
}
