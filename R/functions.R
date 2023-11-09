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
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~round(.x, digits = 1)))
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
