utils::globalVariables(c("roc_auc", ".config", ".pred_class", ".pred_0", "Scaled_Importance",
                         "Importance", "Variable", "std_err", "Type"))
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
#'  - data_split (list). The data split object.
#' @keywords internal
split_data <- function(join_data, ratio = 0.75, seed = 123) {

  set.seed(seed)
  data_split <- rsample::initial_split(join_data, prop = ratio, strata = "Disease")
  train_data <- rsample::training(data_split)
  test_data <- rsample::testing(data_split)

  return(list("train_set" = train_data,
              "test_set" = test_data,
              "data_split" = data_split))
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
#' @keywords internal
filter_sex_specific_disease <- function(control_data,
                                        disease,
                                        diseases,
                                        only_female = NULL,
                                        only_male = NULL) {

  if(!is.null(only_female) & disease %in% only_female) {
    control_data <- control_data |>
      dplyr::filter(Sex == "F")
    diseases_subset <- diseases[!diseases %in% only_male]
  } else if(!is.null(only_male) & disease %in% only_male) {
    control_data <- control_data |>
      dplyr::filter(Sex == "M")
    diseases_subset <- diseases[!diseases %in% only_female]
  } else {
    control_data <- control_data
    diseases_subset <- diseases
  }

  return(list("control_data" = control_data,
              "diseases_subset" = diseases_subset))
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
#' @keywords internal
make_groups <- function(join_data,
                        diseases,
                        only_female = NULL,
                        only_male = NULL,
                        seed = 123) {

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


#' Visualize hyperparameter optimization results
#'
#' This function visualizes the hyperparameter optimization results.
#'
#' @param tune_res (tibble). Hyperparameter optimization results.
#' @param x (character). X-axis variable of the plot.
#' @param color (character). Color variable of the plot. Default is NULL.
#' @param disease (character). Disease to predict.
#'
#' @return hypopt_plot (plot). Hyperparameter optimization plot.
#' @keywords internal
vis_hypopt <- function(tune_res,
                       x,
                       color,
                       disease) {

  hypopt_res <- tune_res |>
    tune::collect_metrics()

  x <- rlang::sym(x)
  mean <- rlang::sym('mean')
  if (!is.null(color)) {
    color <- rlang::sym(color)
  }
  hypopt_plot <- ggplot2::ggplot(hypopt_res, ggplot2::aes(x = !!x, y = !!mean, color = !!color)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - std_err,
                                        ymax = mean + std_err),
                           alpha = 0.5) +
    ggplot2::ggtitle(label = paste0(disease,'')) +
    viridis::scale_color_viridis() +
    ggplot2::labs(y = "metric_mean") +
    ggplot2::theme_classic()

  return(hypopt_plot)
}


#' Hyperparameter optimization for elastic net
#'
#' This function performs hyperparameter optimization for elastic net models.
#' It uses the glmnet engine for logistic regression and tunes either only penalty or both penalty and mixture.
#'
#'
#' @param train_data (list). List of training data sets from make_groups().
#' @param test_data (list). List of testing data sets from make_groups().
#' @param disease (character). Disease to predict.
#' @param type (character). Type of regularization. Default is "lasso". Other options are "ridge" and "elnet".
#' @param cv_sets (numeric). Number of cross-validation sets. Default is 5.
#' @param grid_size (numeric). Size of the grid for hyperparameter optimization. Default is 10.
#' @param ncores (numeric). Number of cores to use for parallel processing. Default is 4.
#' @param hypopt_vis (logical). Whether to visualize hyperparameter optimization results. Default is TRUE.
#' @param exclude_cols (vector). Columns to exclude from the model. Default is NULL.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with three elements:
#'  - elnet_tune (tibble). Hyperparameter optimization results.
#'  - wf (workflow). Workflow object.
#'  - train_set (tibble). Training set.
#'  - test_set (tibble). Testing set.
#' @keywords internal
elnet_hypopt <- function(train_data,
                         test_data,
                         disease,
                         type = "lasso",
                         cv_sets = 5,
                         grid_size = 10,
                         ncores = 4,
                         hypopt_vis = TRUE,
                         exclude_cols = NULL,
                         seed = 123) {

  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  }

  # Prepare train data and create cross-validation sets with binary classifier
  train_set <- train_data[[disease]] |>
    dplyr::mutate(Disease = ifelse(Disease == disease, 1, 0)) |>
    dplyr::mutate(Disease = as.factor(Disease)) |>
    dplyr::select(-dplyr::any_of(exclude_cols))
  train_folds <- rsample::vfold_cv(train_set, v = cv_sets, strata = Disease)

  test_set <- test_data[[disease]] |>
    dplyr::mutate(Disease = ifelse(Disease == disease, 1, 0)) |>
    dplyr::mutate(Disease = as.factor(Disease)) |>
    dplyr::select(-dplyr::any_of(exclude_cols))

  elnet_rec <- recipes::recipe(Disease ~ ., data = train_set) |>
    recipes::update_role(DAid, new_role = "id") |>
    recipes::step_normalize(recipes::all_numeric()) |>
    recipes::step_nzv(recipes::all_numeric()) |>
    recipes::step_corr(recipes::all_numeric()) |>
    recipes::step_impute_knn(recipes::all_numeric())

  if (type == "elnet") {
    elnet_spec <- parsnip::logistic_reg(
      penalty = tune::tune(),  # lambda
      mixture = tune::tune()  # alpha
    ) |>
      parsnip::set_engine("glmnet")
  } else if (type == "lasso") {
    elnet_spec <- parsnip::logistic_reg(
      penalty = tune::tune(),
      mixture = 1
    ) |>
      parsnip::set_engine("glmnet")
  } else if (type == "ridge") {
    elnet_spec <- parsnip::logistic_reg(
      penalty = tune::tune(),
      mixture = 0
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

  if (hypopt_vis) {
    if (type == "elnet") {
      hypopt_plot <- vis_hypopt(elnet_tune, "penalty", "mixture", disease)
    } else {
      hypopt_plot <- vis_hypopt(elnet_tune, "penalty", NULL, disease)
    }
    return(list("elnet_tune" = elnet_tune,
                "elnet_wf" = elnet_wf,
                "train_set" = train_set,
                "test_set" = test_set,
                "hyperopt_vis" = hypopt_plot))
  }

  return(list("elnet_tune" = elnet_tune,
              "elnet_wf" = elnet_wf,
              "train_set" = train_set,
              "test_set" = test_set))
}


#' Hyperparameter optimization for random forest
#'
#' This function performs hyperparameter optimization for random forest models.
#' It uses the ranger engine for logistic regression and tunes the number of trees and
#' the number of variables randomly sampled at each split.
#'
#' @param train_data (list). List of training data sets from make_groups().
#' @param test_data (list). List of testing data sets from make_groups().
#' @param disease (character). Disease to predict.
#' @param cv_sets (numeric). Number of cross-validation sets. Default is 5.
#' @param grid_size (numeric). Size of the grid for hyperparameter optimization. Default is 10.
#' @param ncores (numeric). Number of cores to use for parallel processing. Default is 4.
#' @param hypopt_vis (logical). Whether to visualize hyperparameter optimization results. Default is TRUE.
#' @param exclude_cols (vector). Columns to exclude from the model. Default is NULL.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with four elements:
#'  - rf_tune (tibble). Hyperparameter optimization results.
#'  - rf_wf (workflow). Workflow object.
#'  - train_set (tibble). Training set.
#'  - test_set (tibble). Testing set.
#'  - hyperopt_vis (plot). Hyperparameter optimization plot.
#' @keywords internal
rf_hypopt <- function(train_data,
                      test_data,
                      disease,
                      cv_sets = 5,
                      grid_size = 10,
                      ncores = 4,
                      hypopt_vis = TRUE,
                      exclude_cols = NULL,
                      seed = 123) {

  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  }

  # Prepare train data and create cross-validation sets with binary classifier
  train_set <- train_data[[disease]] |>
    dplyr::mutate(Disease = ifelse(Disease == disease, 1, 0)) |>
    dplyr::mutate(Disease = as.factor(Disease)) |>
    dplyr::select(-dplyr::any_of(exclude_cols)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character) & !dplyr::all_of("DAid"), as.factor))  # Solve bug of Gower distance function

  train_folds <- rsample::vfold_cv(train_set, v = cv_sets, strata = Disease)

  test_set <- test_data[[disease]] |>
    dplyr::mutate(Disease = ifelse(Disease == disease, 1, 0)) |>
    dplyr::mutate(Disease = as.factor(Disease)) |>
    dplyr::select(-dplyr::any_of(exclude_cols)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character) & !dplyr::all_of("DAid"), as.factor))

  rf_rec <- recipes::recipe(Disease ~ ., data = train_set) |>
    recipes::update_role(DAid, new_role = "id") |>
    recipes::step_normalize(recipes::all_numeric()) |>
    recipes::step_nzv(recipes::all_numeric()) |>
    recipes::step_corr(recipes::all_numeric()) |>
    recipes::step_impute_knn(recipes::all_numeric())

  rf_spec <- parsnip::rand_forest(
    trees = 1000,
    mtry = tune::tune(),
    min_n = tune::tune()
  ) |>
    parsnip::set_mode("classification") |>
    parsnip::set_engine("ranger", importance = "permutation")

  disease_pred <- train_set |> dplyr::select(-dplyr::any_of(c("Disease", "DAid", "Sex", "Age", "BMI")))

  rf_wf <- workflows::workflow() |>
    workflows::add_model(rf_spec) |>
    workflows::add_recipe(rf_rec)

  rf_grid <- rf_wf |>
    workflows::extract_parameter_set_dials() |>
    dials::finalize(disease_pred) |>
    dials::grid_latin_hypercube(size = grid_size)

  roc_res <- yardstick::metric_set(yardstick::roc_auc)

  set.seed(seed)
  ctrl <- tune::control_grid(save_pred = TRUE, parallel_over = "everything")
  rf_tune <- rf_wf |>
    tune::tune_grid(
      train_folds,
      grid = rf_grid,
      control = ctrl,
      metrics = roc_res
    )

  if (hypopt_vis) {
    hypopt_plot <- vis_hypopt(rf_tune, "min_n", "mtry", disease)

    return(list("rf_tune" = rf_tune,
                "rf_wf" = rf_wf,
                "train_set" = train_set,
                "test_set" = test_set,
                "hyperopt_vis" = hypopt_plot))
  }

  return(list("rf_tune" = rf_tune,
              "rf_wf" = rf_wf,
              "train_set" = train_set,
              "test_set" = test_set))
}


#' Fit the final model
#'
#' This function fits the final model using the best hyperparameters from hyperparameter optimization.
#'
#' @param train_set (tibble). Training set.
#' @param tune_res (tibble). Hyperparameter optimization results.
#' @param wf (workflow). Workflow object.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with three elements:
#'  - final_elnet (parsnip model). Final model.
#'  - best_elnet (tibble). Best hyperparameters from hyperparameter optimization.
#'  - final_wf (workflow). Final workflow object.
#' @keywords internal
finalfit <- function(train_set,
                     tune_res,
                     wf,
                     seed = 123) {

  best <- tune_res |>
    tune::select_best(metric = "roc_auc") |>
    dplyr::select(-.config)

  final_wf <- tune::finalize_workflow(wf, best)

  final <- final_wf |>
    parsnip::fit(train_set)

  return(list("final" = final,
              "best" = best,
              "final_wf" = final_wf))
}


#' Test the final model
#'
#' This function tests the final model on the test set.
#' It calculates the accuracy, sensitivity, specificity, AUC, and confusion matrix.
#'
#' @param train_set (tibble). Training set.
#' @param test_set (tibble). Testing set.
#' @param disease (character). Disease to predict.
#' @param finalfit_res (list). Results from elnet_finalfit().
#' @param exclude_cols (vector). Columns to exclude from the model. Default is NULL.
#' @param type (character). Type of regularization. Default is "lasso". Other options are "ridge" and "elnet".
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is NULL.
#'
#' @return A list with two elements:
#'  - metrics (list). A list with 5 metrics:
#'   - accuracy (numeric). Accuracy of the model.
#'   - sensitivity (numeric). Sensitivity of the model.
#'   - specificity (numeric). Specificity of the model.
#'   - auc (numeric). AUC of the model.
#'   - conf_matrix (tibble). Confusion matrix of the model.
#'   - roc_curve (tibble). ROC curve of the model.
#'  - mixture (numeric). Mixture of lasso and ridge regularization.
#' @keywords internal
testfit <- function(train_set,
                    test_set,
                    disease,
                    finalfit_res,
                    exclude_cols = NULL,
                    type = "lasso",
                    seed = 123,
                    palette = NULL) {

  set.seed(seed)
  splits <- rsample::make_splits(train_set, test_set)

  preds <- tune::last_fit(finalfit_res$final_wf,
                          splits,
                          metrics = yardstick::metric_set(yardstick::roc_auc))
  res <- stats::predict(finalfit_res$final, new_data = test_set)

  res <- dplyr::bind_cols(res, test_set |> dplyr::select(Disease))

  accuracy <- res |>
    yardstick::accuracy(Disease, .pred_class)

  sensitivity <- res |>
    yardstick::sensitivity(Disease, .pred_class)

  specificity <- res |>
    yardstick::specificity(Disease, .pred_class)

  if (is.null(names(palette)) && !is.null(palette)) {
    disease_color <- get_hpa_palettes()[[palette]][[disease]]
  } else if (!is.null(palette)) {
    disease_color <- palette
  } else {
    disease_color <- "black"
  }

  selected_point <- tibble::tibble(x = 1 - specificity$.estimate,
                                   y = sensitivity$.estimate)
  roc <- preds |>
    tune::collect_predictions(summarize = F) |>
    yardstick::roc_curve(truth = Disease, .pred_0) |>
    ggplot2::ggplot(ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
    ggplot2::geom_path(colour = disease_color, size = 2) +
    ggplot2::geom_point(data = selected_point, ggplot2::aes(x = x, y = y), size = 2, shape = 4, colour = "black") +
    ggplot2::geom_abline(lty = 3) +
    ggplot2::coord_equal() +
    theme_hpa()

  auc <- preds |>
    tune::collect_metrics()

  cm <- res |>
    yardstick::conf_mat(Disease, .pred_class)

  if (type == "elnet") {
    mixture <- finalfit_res$best_elnet$mixture
  } else if (type == "lasso") {
    mixture <- 1
  } else if (type == "ridge") {
    mixture <- 0
  } else {
    mixture <- NA
  }

  return(list("metrics" = list("accuracy" = round(accuracy$.estimate, 2),
                              "sensitivity" = round(sensitivity$.estimate, 2),
                              "specificity" = round(specificity$.estimate, 2),
                              "auc" = round(auc$.estimate, 2),
                              "conf_matrix" = cm,
                              "roc_curve" = roc),
              "mixture" = mixture))
}


#' Generate subtitle for variable importance plot
#'
#' This function generates a subtitle for the variable importance plot.
#'
#' @param features (tibble). Features with importance values.
#' @param accuracy (numeric). Accuracy of the model.
#' @param sensitivity (numeric). Sensitivity of the model.
#' @param specificity (numeric). Specificity of the model.
#' @param auc (numeric). AUC of the model.
#' @param mixture (numeric). Mixture of lasso and ridge regularization.
#' @param subtitle (vector). Vector of subtitles to include in the plot. Default is all.
#'
#' @return subtitle (character). Subtitle for the plot.
#' @keywords internal
generate_subtitle <- function(features, accuracy, sensitivity, specificity, auc, mixture, subtitle = NULL) {
  subtitle_parts <- c()

  if ("accuracy" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('accuracy = ', round(accuracy, 2), '    '))
  }

  if ("sensitivity" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('sensitivity = ', round(sensitivity, 2), '    '))
  }

  if ("specificity" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('specificity = ', round(specificity, 2), '    '))
  }

  if ("auc" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('AUC = ', round(auc, 2), '    '))
  }

  if (length(subtitle_parts) > 0) {
    subtitle_parts <- c(subtitle_parts, '\n')
  }

  if ("features" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('Features = ', nrow(features), '    '))
  }

  if ("top-features" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('top-features = ',
                                               nrow(features |> dplyr::filter(Scaled_Importance >= 50)),
                                               '    '))
  }

  if ("mixture" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('Lasso/Ridge ratio = ', round(mixture, 2), '    '))
  }

  subtitle <- paste(subtitle_parts, collapse = '')

  return(subtitle)
}


#' Plot variable importance
#'
#' This function collects the features and their importance.
#' It scales their importance and plots it against them.
#'
#' @param finalfit_res (list). Results from finalfit().
#' @param disease (character). Disease to predict.
#' @param accuracy (numeric). Accuracy of the model.
#' @param sensitivity (numeric). Sensitivity of the model.
#' @param specificity (numeric). Specificity of the model.
#' @param auc (numeric). AUC of the model.
#' @param mixture (numeric). Mixture of lasso and ridge regularization.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is NULL.
#' @param vline (logical). Whether to add a vertical line at 50% importance. Default is TRUE.
#' @param subtitle (vector). Vector of subtitles to include in the plot. Default is a list with all.
#'
#' @return A list with two elements:
#'  - features (tibble). Features with importance values.
#'  - var_imp_plot (plot). Variable importance plot.
#' @keywords internal
plot_var_imp <- function (finalfit_res,
                          disease,
                          accuracy,
                          sensitivity,
                          specificity,
                          auc,
                          mixture,
                          palette = NULL,
                          vline = T,
                          subtitle = c("accuracy",
                                       "sensitivity",
                                       "specificity",
                                       "auc",
                                       "features",
                                       "top-features",
                                       "mixture")) {

  features <- finalfit_res$final |>
    workflows::extract_fit_parsnip() |>
    vip::vi(lambda = finalfit_res$best_model$penalty,
            alpha = finalfit_res$best_model$mixture
    ) |>
    dplyr::mutate(
      Importance = abs(Importance),
      Variable = forcats::fct_reorder(Variable, Importance)
    ) |>
    dplyr::arrange(dplyr::desc(Importance)) |>
    dplyr::mutate(Scaled_Importance = scales::rescale(Importance, to = c(0, 100))) |>
    dplyr::filter(Scaled_Importance > 0)

  subtitle_text <- generate_subtitle(features, accuracy, sensitivity, specificity, auc, mixture, subtitle)

  # Prepare palettes
  pals <- get_hpa_palettes()
  if (!is.null(palette) && is.null(names(palette))) {
    pal <- pals[palette]
    pal <- unlist(pals[[palette]])
  } else if (!is.null(palette)) {
    pal <- palette
  } else {
    pal <- "red3"
  }

  var_imp_plot <- features |>
    ggplot2::ggplot(ggplot2::aes(x = Scaled_Importance, y = Variable)) +
    ggplot2::geom_col(ggplot2::aes(fill = ifelse(Scaled_Importance > 50, disease, NA))) +
    ggplot2::labs(y = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +  # Keep x-axis tick labels at 0 and 100
    ggplot2::scale_fill_manual(values = pal, na.value = "grey50") +
    ggplot2::ggtitle(label = paste0(disease, ''),
                     subtitle = subtitle_text) +
    ggplot2::xlab('Importance') +
    ggplot2::ylab('Features') +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none")

  if (isTRUE(vline)) {
    var_imp_plot <- var_imp_plot +
      ggplot2::geom_vline(xintercept = 50, linetype = 'dashed', color = 'black')
  }

  return(list("features" = features,
              "var_imp_plot" = var_imp_plot))
}


#' Elastic net classification model pipeline
#'
#' This function runs the elastic net classification model pipeline.
#' It splits the data into training and test sets, creates class-balanced groups, and fits the model.
#' It also performs hyperparameter optimization, fits the final model, tests it, and plots useful visualizations.
#'
#' @param olink_data (tibble). Olink data.
#' @param metadata (tibble). Metadata.
#' @param wide (logical). Whether the data is wide format. Default is FALSE.
#' @param only_female (vector). Vector of diseases that are female specific. Default is NULL.
#' @param only_male (vector). Vector of diseases that are male specific. Default is NULL.
#' @param exclude_cols (vector). Columns to exclude from the model. Default is "Sex".
#' @param ratio (numeric). Ratio of training data to test data. Default is 0.75.
#' @param type (character). Type of regularization. Default is "lasso". Other options are "ridge" and "elnet".
#' @param cv_sets (numeric). Number of cross-validation sets. Default is 5.
#' @param grid_size (numeric). Size of the grid for hyperparameter optimization. Default is 10.
#' @param ncores (numeric). Number of cores to use for parallel processing. Default is 4.
#' @param hypopt_vis (logical). Whether to visualize hyperparameter optimization results. Default is TRUE.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is NULL.
#' @param vline (logical). Whether to add a vertical line at 50% importance. Default is TRUE.
#' @param subtitle (vector). Vector of subtitles to include in the plot. Default is a list with all.
#' @param nfeatures (numeric). Number of top features to include in the boxplot. Default is 9.
#' @param points (logical). Whether to add points to the boxplot. Default is TRUE.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with results for each disease. The list contains:
#'  - hypopt_res (list). Hyperparameter optimization results.
#'  - finalfit_res (list). Final model fitting results.
#'  - testfit_res (list). Test model fitting results.
#'  - var_imp_res (list). Variable importance results.
#' @export
#'
#' @examples
#' unique_samples <- unique(example_data$Sample)
#' filtered_data <- example_data |>
#'  dplyr::filter(Sample %in% unique_samples[1:148])
#'
#' res <- do_elnet(filtered_data,
#'                 example_metadata,
#'                 palette = "cancers12",
#'                 cv_sets = 2,
#'                 grid_size = 1,
#'                 ncores = 1)
do_elnet <- function(olink_data,
                     metadata,
                     wide = F,
                     only_female = NULL,
                     only_male = NULL,
                     exclude_cols = "Sex",
                     ratio = 0.75,
                     type = "lasso",
                     cv_sets = 5,
                     grid_size = 10,
                     ncores = 4,
                     hypopt_vis = TRUE,
                     palette = NULL,
                     vline = T,
                     subtitle = c("accuracy",
                                  "sensitivity",
                                  "specificity",
                                  "auc",
                                  "features",
                                  "top-features",
                                  "mixture"),
                     nfeatures = 9,
                     points = T,
                     seed = 123) {

  # Prepare datasets
  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
  join_data <- wide_data |>
    dplyr::left_join(metadata |> dplyr::select(DAid, Disease, Sex))
  diseases <- unique(metadata$Disease)

  # Prepare sets and groups
  data_split <- split_data(join_data, ratio, seed)
  diseases <- unique(join_data$Disease)
  train_list <- make_groups(data_split$train_set,
                            diseases,
                            only_female,
                            only_male,
                            seed)
  test_list <- make_groups(data_split$test_set,
                           diseases,
                           only_female,
                           only_male,
                           seed)

  message("Sets and groups are ready. Model fitting is starting...")

  # Run model
  elnet_results <- lapply(diseases, function(disease) {
    message(paste0("Classification model for ", disease, " is starting..."))
    hypopt_res <- elnet_hypopt(train_list,
                               test_list,
                               disease,
                               type,
                               cv_sets,
                               grid_size,
                               ncores,
                               hypopt_vis,
                               exclude_cols,
                               seed)

    finalfit_res <- finalfit(hypopt_res$train_set,
                             hypopt_res$elnet_tune,
                             hypopt_res$elnet_wf,
                             seed)

    testfit_res <- testfit(hypopt_res$train_set,
                           hypopt_res$test_set,
                           disease,
                           finalfit_res,
                           exclude_cols,
                           type,
                           seed,
                           palette)

    var_imp_res <- plot_var_imp(finalfit_res,
                                disease,
                                testfit_res$metrics$accuracy,
                                testfit_res$metrics$sensitivity,
                                testfit_res$metrics$specificity,
                                testfit_res$metrics$auc,
                                testfit_res$mixture,
                                palette = palette,
                                vline = vline,
                                subtitle)

    top_features <- var_imp_res$features |>
      dplyr::arrange(dplyr::desc(Scaled_Importance)) |>
      dplyr::select(Variable) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
      utils::head(nfeatures)
    proteins <- top_features[['Variable']]

    boxplot_res <- plot_protein_boxplot(join_data,
                                        proteins,
                                        disease,
                                        points,
                                        palette)

    return(list("hypopt_res" = hypopt_res,
                "finalfit_res" = finalfit_res,
                "testfit_res" = testfit_res,
                "var_imp_res" = var_imp_res,
                "boxplot_res" = boxplot_res))
  })

  names(elnet_results) <- diseases

  return(elnet_results)
}


#' Random forest classification model pipeline
#'
#' This function runs the random forest classification model pipeline.
#' It splits the data into training and test sets, creates class-balanced groups, and fits the model.
#' It also performs hyperparameter optimization, fits the final model, tests it, and plots useful visualizations.
#'
#' @param olink_data (tibble). Olink data.
#' @param metadata (tibble). Metadata.
#' @param wide (logical). Whether the data is wide format. Default is FALSE.
#' @param only_female (vector). Vector of diseases that are female specific. Default is NULL.
#' @param only_male (vector). Vector of diseases that are male specific. Default is NULL.
#' @param exclude_cols (vector). Columns to exclude from the model. Default is "Sex".
#' @param ratio (numeric). Ratio of training data to test data. Default is 0.75.
#' @param cv_sets (numeric). Number of cross-validation sets. Default is 5.
#' @param grid_size (numeric). Size of the grid for hyperparameter optimization. Default is 10.
#' @param ncores (numeric). Number of cores to use for parallel processing. Default is 4.
#' @param hypopt_vis (logical). Whether to visualize hyperparameter optimization results. Default is TRUE.
#' @param palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is NULL.
#' @param vline (logical). Whether to add a vertical line at 50% importance. Default is TRUE.
#' @param subtitle (vector). Vector of subtitles to include in the plot. Default is a list with all except mixture.
#' @param nfeatures (numeric). Number of top features to include in the boxplot. Default is 9.
#' @param points (logical). Whether to add points to the boxplot. Default is TRUE.
#' @param seed (numeric). Seed for reproducibility. Default is 123.
#'
#' @return A list with results for each disease. The list contains:
#'  - hypopt_res (list). Hyperparameter optimization results.
#'  - finalfit_res (list). Final model fitting results.
#'  - testfit_res (list). Test model fitting results.
#'  - var_imp_res (list). Variable importance results.
#' @export
#'
#' @examples
#' unique_samples <- unique(example_data$Sample)
#' filtered_data <- example_data |>
#'  dplyr::filter(Sample %in% unique_samples[1:148])
#'
#' res <- do_rf(filtered_data,
#'              example_metadata,
#'              palette = "cancers12",
#'              cv_sets = 2,
#'              grid_size = 1,
#'              ncores = 1)
do_rf <- function(olink_data,
                  metadata,
                  wide = F,
                  only_female = NULL,
                  only_male = NULL,
                  exclude_cols = "Sex",
                  ratio = 0.75,
                  cv_sets = 5,
                  grid_size = 10,
                  ncores = 4,
                  hypopt_vis = TRUE,
                  palette = NULL,
                  vline = T,
                  subtitle = c("accuracy",
                               "sensitivity",
                               "specificity",
                               "auc",
                               "features",
                               "top-features"),
                  nfeatures = 9,
                  points = T,
                  seed = 123) {

  # Prepare datasets
  if (isFALSE(wide)) {
    wide_data <- widen_data(olink_data)
  } else {
    wide_data <- olink_data
  }
  join_data <- wide_data |>
    dplyr::left_join(metadata |> dplyr::select(DAid, Disease, Sex))
  diseases <- unique(metadata$Disease)

  # Prepare sets and groups
  data_split <- split_data(join_data, ratio, seed)
  diseases <- unique(join_data$Disease)
  train_list <- make_groups(data_split$train_set,
                            diseases,
                            only_female,
                            only_male,
                            seed)
  test_list <- make_groups(data_split$test_set,
                           diseases,
                           only_female,
                           only_male,
                           seed)

  message("Sets and groups are ready. Model fitting is starting...")

  # Run model
  rf_results <- lapply(diseases, function(disease) {
    message(paste0("Classification model for ", disease, " is starting..."))
    hypopt_res <- rf_hypopt(train_list,
                            test_list,
                            disease,
                            cv_sets,
                            grid_size,
                            ncores,
                            hypopt_vis,
                            exclude_cols,
                            seed)

    finalfit_res <- finalfit(hypopt_res$train_set,
                             hypopt_res$rf_tune,
                             hypopt_res$rf_wf,
                             seed)

    testfit_res <- testfit(hypopt_res$train_set,
                           hypopt_res$test_set,
                           disease,
                           finalfit_res,
                           exclude_cols,
                           type = "other",
                           seed,
                           palette)

    var_imp_res <- plot_var_imp(finalfit_res,
                                disease,
                                testfit_res$metrics$accuracy,
                                testfit_res$metrics$sensitivity,
                                testfit_res$metrics$specificity,
                                testfit_res$metrics$auc,
                                testfit_res$mixture,
                                palette = palette,
                                vline = vline,
                                subtitle)

    top_features <- var_imp_res$features |>
      dplyr::arrange(dplyr::desc(Scaled_Importance)) |>
      dplyr::select(Variable) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
      utils::head(nfeatures)
    proteins <- top_features[['Variable']]

    boxplot_res <- plot_protein_boxplot(join_data,
                                        proteins,
                                        disease,
                                        points,
                                        palette)

    return(list("hypopt_res" = hypopt_res,
                "finalfit_res" = finalfit_res,
                "testfit_res" = testfit_res,
                "var_imp_res" = var_imp_res,
                "boxplot_res" = boxplot_res))
  })

  names(rf_results) <- diseases

  return(rf_results)
}


#' Plot protein features summary
#'
#' This function plots the number of proteins and the number of top proteins for each disease.
#' It also plots the upset plot of the top or all proteins.
#'
#' @param ml_results (list). Results from do_elnet() or do_rf().
#' @param importance (numeric). Importance threshold for top features. Default is 50.
#' @param upset_top_features (logical). Whether to plot the upset plot for the top features. Default is FALSE.
#' @param disease_palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is NULL.
#' @param feature_type_palette (character or vector). The color palette for the plot. If it is a character, it should be one of the palettes from get_hpa_palettes(). Default is "all-features" = "pink" and "top-features" = "darkblue".
#'
#' @return A list with two elements:
#'   - features_barplot (plot). Barplot of the number of proteins and top proteins for each disease.
#'   - upset_plot_features (plot). Upset plot of the top or all proteins.
#' @export
#'
#' @examples
#' unique_samples <- unique(example_data$Sample)
#' filtered_data <- example_data |>
#'  dplyr::filter(Sample %in% unique_samples[1:148])
#' res <- do_elnet(filtered_data,
#'                 example_metadata,
#'                 cv_sets = 2,
#'                 grid_size = 1,
#'                 ncores = 1)
#'
#' plot <- plot_features_summary(res)
plot_features_summary <- function(ml_results,
                                  importance = 50,
                                  upset_top_features = F,
                                  disease_palette = NULL,
                                  feature_type_palette = c("all-features" = "pink", "top-features" = "darkblue")) {

  barplot_data <- lapply(names(ml_results), function(disease) {

    features <- ml_results[[disease]]$var_imp_res$features |>
      dplyr::mutate(Disease = disease) |>
      dplyr::select(Disease, Variable) |>
      dplyr::rename(Assay = Variable) |>
      dplyr::group_by(Disease) |>
      dplyr::summarise(Count = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(Type = "all-features")

    top_features <- ml_results[[disease]]$var_imp_res$features |>
      dplyr::mutate(Disease = disease) |>
      dplyr::filter(Scaled_Importance >= importance) |>
      dplyr::select(Disease, Variable) |>
      dplyr::rename(Assay = Variable) |>
      dplyr::group_by(Disease) |>
      dplyr::summarise(Count = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(Type = "top-features")

    barplot_data <- rbind(features, top_features)
  })

  barplot_data <- do.call(rbind, barplot_data)

  features_barplot <- barplot_data |>
    ggplot2::ggplot(ggplot2::aes(x = Disease, y = Count, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(x = "", y = "Number of protein", fill = "Feature type") +
    theme_hpa(angled = T) +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_text(face = "bold"))

  if (is.null(names(feature_type_palette)) && !is.null(feature_type_palette)) {
    features_barplot <- features_barplot + scale_fill_hpa(feature_type_palette)
  } else if (!is.null(feature_type_palette)) {
    features_barplot <- features_barplot + ggplot2::scale_fill_manual(values = feature_type_palette)
  }

  upset_features <- lapply(names(ml_results), function(disease) {

    if (upset_top_features == T) {
      upset_features <- ml_results[[disease]]$var_imp_res$features |>
        dplyr::filter(Scaled_Importance >= importance) |>
        dplyr::pull(Variable)
    } else {
      upset_features <- ml_results[[disease]]$var_imp_res$features |>
        dplyr::pull(Variable)
    }

  })
  names(upset_features) <- names(ml_results)

  if (is.null(names(disease_palette)) && !is.null(disease_palette)) {
    pal <- get_hpa_palettes()[[disease_palette]]
  } else if (!is.null(disease_palette)) {
    pal <- disease_palette
  } else {
    pal <- "black"
  }

  upset_plot_features <- UpSetR::upset(UpSetR::fromList(upset_features),
                                 order.by = "freq",
                                 nsets = length(names(ml_results)),
                                 sets.bar.color = pal)


  return(list("features_barplot" = features_barplot, "upset_plot_features" = upset_plot_features))
}
