#' Calculate descriptive statistics of each metabolite.
#'
#' @param data Lipidomics dataset.
#'
#' @return A data.frame/tibble.
#'
descriptive_stats <- function(data) {
  data %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, digits = 1)))
}

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


#' Convert column value string into snakecase.
#'
#' @param data Data with string columns
#' @param cols The column to convert into snakecase
#'
#' @return A data frame.
#'

column_values_to_snake_case <- function(data, cols) {
  data %>%
    dplyr::mutate(dplyr::across({{ cols }}, snakecase::to_snake_case))
}

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


#' A transformation recipe to pre-process the data
#'
#' @param data The lipidomics dataset.
#' @param metabolite_variable The column of the metabolite variable.
#'
#' @return Recipe with specifications

create_recipe_spec <- function(data, metabolite_variable) {
    recipes::recipe(data) %>%
        recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") %>%
        recipes::update_role(class, new_role = "outcome") %>%
        recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}





