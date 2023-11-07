#| code-fold: true
#| code-summary: "**Click for a potential solution**. Only click if you are struggling or are out of time."
#' Calculate descriptive statistics of each metabolite.
#'
#' @param data Lipidomics dataset.
#'
#' @return A data.frame/tibble.
#'
descriptive_stats <- function(data) {
  data %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(value, list(
      mean = mean,
      sd = sd,
      median = median,
      iqr = IQR
    ))) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~round(.x, digits = 1)))
}

#' Generate the model variation results using bootstrap on a single metabolite.
#'
#' @param data The lipidomics data.
#'
#' @return A nested tibble.
#'
generate_model_variation <- function(data) {
  create_model_workflow(
    parsnip::logistic_reg() %>%
      parsnip::set_engine("glm"),
    data %>%
      create_recipe_spec(tidyselect::starts_with("metabolite_"))
  ) %>%
    tune::fit_resamples(
      resamples = rsample::bootstraps(data, times = 10),
      control = tune::control_resamples(
        extract = tidy_model_output,
        save_pred = TRUE
      )
    )
}

#' Tidy up the bootstrap output.
#'
#' @param bootstrap_results The bootstrap object with model results.
#'
#' @return A data frame.
#'
tidy_bootstrap_output <- function(bootstrap_results) {
  bootstrap_results %>%
    dplyr::select(id, .extracts) %>%
    # Need to unnest twice since first `.extracts` is a nest of another two
    # columns of `.extracts` and `.config`.
    tidyr::unnest(cols = .extracts) %>%
    tidyr::unnest(cols = .extracts) %>%
    dplyr::filter(stringr::str_detect(term, "metabolite_")) %>%
    add_original_metabolite_names(lipidomics)
}

#| eval: false
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
#' Calculate the uncertainty in results.
#'
#' @param data The lipidomics data.
#'
#' @return A data frame (or file path)
#'
calculate_variation <- function(data) {
  data %>%
    split_by_metabolite() %>%
    purrr::map(generate_model_variation) %>%
    purrr::map(tidy_bootstrap_output) %>% 
    purrr::list_rbind()
}

#' Plot the uncertainty in the estimates of the models.
#'
#' @param model_results The model results with the variation.
#'
#' @return A ggplot2 image.
#'
plot_variation <- function(model_results) {
  model_results %>%
    ggplot2::ggplot(ggplot2::aes(x = estimate)) +
    ggplot2::geom_dotplot() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
}

#| eval: false
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
## This should be in the R/functions.R file.
#' Plot for basic distribution of metabolite data.
#'
#' @param data The lipidomics dataset.
#'
#' @return A ggplot2 graph.
#'
plot_distributions <- function(data) {
  data %>% 
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
}

#' Convert column value strings into snakecase.
#'
#' @param data Data with string columns.
#' @param cols The column to convert into snakecase.
#'
#' @return A data frame.
#'
column_values_to_snake_case <- function(data, cols) {
  data %>%
    dplyr::mutate(dplyr::across({{ cols }}, snakecase::to_snake_case))
}

#| eval: true
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
#' Convert the metabolite long format into a wider one.
#'
#' @param data The lipidomics dataset.
#'
#' @return A wide data frame.
#'
metabolites_to_wider <- function(data) {
  data %>%
    tidyr::pivot_wider(
      names_from = metabolite,
      values_from = value,
      values_fn = mean,
      names_prefix = "metabolite_"
    )
}

#' A transformation recipe to pre-process the data.
#'
#' @param data The lipidomics dataset.
#' @param metabolite_variable The column of the metabolite variable.
#'
#' @return
#'
create_recipe_spec <- function(data, metabolite_variable) {
  recipes::recipe(data) %>%
    recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") %>%
    recipes::update_role(class, new_role = "outcome") %>%
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}

#' Create a workflow object of the model and transformations.
#'
#' @param model_specs The model specs
#' @param recipe_specs The recipe specs
#'
#' @return A workflow object
#'
create_model_workflow <- function(model_specs, recipe_specs) {
  workflows::workflow() %>%
    workflows::add_model(model_specs) %>%
    workflows::add_recipe(recipe_specs)
}

#' Create a tidy output of the model results.
#'
#' @param workflow_fitted_model The model workflow object that has been fitted.
#'
#' @return A data frame.
#'
tidy_model_output <- function(workflow_fitted_model) {
  workflow_fitted_model %>%
    workflows::extract_fit_parsnip() %>%
    broom::tidy(exponentiate = TRUE)
}

#' Convert the long form dataset into a list of wide form data frames.
#'
#' @param data The lipidomics dataset.
#'
#' @return A list of data frames.
#'
split_by_metabolite <- function(data) {
  data %>%
    column_values_to_snake_case(metabolite) %>%
    dplyr::group_split(metabolite) %>%
    purrr::map(metabolites_to_wider)
}

#' Generate the results of a model
#'
#' @param data The lipidomics dataset.
#'
#' @return A data frame.
#'
generate_model_results <- function(data) {
  create_model_workflow(
    parsnip::logistic_reg() %>%
      parsnip::set_engine("glm"),
    data %>%
      create_recipe_spec(tidyselect::starts_with("metabolite_"))
  ) %>%
    parsnip::fit(data) %>%
    tidy_model_output()
}

#| eval: false
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
#' Add the original metabolite names (not as snakecase) to the model results.
#'
#' @param model_results The data frame with the model results.
#' @param data The original lipidomics dataset.
#'
#' @return A data frame.
#'
add_original_metabolite_names <- function(model_results, data) {
  data %>%
    dplyr::mutate(term = metabolite) %>%
    column_values_to_snake_case(term) %>%
    dplyr::mutate(term = stringr::str_c("metabolite_", term)) %>%
    dplyr::distinct(term, metabolite) %>%
    dplyr::right_join(model_results, by = "term")
}

#| eval: false
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
#' Calculate the estimates for the model for each metabolite.
#'
#' @param data The lipidomics dataset.
#'
#' @return A data frame.
#'
calculate_estimates <- function(data) {
  data %>%
    column_values_to_snake_case(metabolite) %>%
    dplyr::group_split(metabolite) %>%
    purrr::map(metabolites_to_wider) %>%
    purrr::map(generate_model_results) %>%
    purrr::list_rbind() %>% 
    dplyr::filter(stringr::str_detect(term, "metabolite_")) %>%
    add_original_metabolite_names(data)
}

#| eval: false
#| code-fold: true
#| code-summary: "**Click for the solution**. Only click if you are struggling or are out of time."
#' Plot the estimates and standard errors of the model results.
#'
#' @param results The model estimate results.
#'
#' @return A ggplot2 figure.
#'
plot_estimates <- function(results) {
  results %>%
    ggplot2::ggplot(ggplot2::aes(
      x = estimate, y = metabolite,
      xmin = estimate - std.error,
      xmax = estimate + std.error
    )) +
    ggplot2::geom_pointrange() +
    ggplot2::coord_fixed(xlim = c(0, 5))
}

