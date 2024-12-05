#' Add Fitted Distribution to Data Frame
#'
#' This function fits a specified distribution to the data in a given data frame and adds the fitted distribution as a new column.
#'
#' @param cost_df A data frame containing the data to which the distribution will be fitted. The data frame should have columns named "mean", "lb_95", and "ub_95".
#' @param dist A character string specifying the distribution to fit. The distribution should be one of the distributions supported by the `q` function in R (e.g., "norm" for normal distribution, "tri" for triangular distribution).
#'
#' @return A data frame with an additional column named `fitt_dist` containing the fitted distribution for each row.
#'
#' @examples
#' \dontrun{
#' cost_df <- data.frame(mean = c(10, 20), lb_95 = c(5, 15), ub_95 = c(15, 25))
#' add_fitted_dist(cost_df, "norm")
#' }
#' 
#' @importFrom stats optim
#' @importFrom purrr map_dbl
#' @importFrom dplyr %>%
#' @export
add_fitted_dist <- function(cost_df, dist) {
    fitdist_ci <- function(pars, data, dist) {
        qs <- c(0.025, 0.5, 0.975); 1:3 %>% map_dbl(~ (dist(qs[.x], pars[1], pars[2]) - data[.x])^2) %>% sum
    }
    dist_q <- str2expression(paste0("q", dist))
    init <- c(10, 10)
    fitt_dist <- vector()
    for(i in 1:nrow(cost_df)) {
        vals <- cost_df[i, ]
        if (is.na(vals$lb_95) | is.na(vals$ub_95)) {
            fitt_dist[i] <- "NA"
        } else {
            if (dist == "tri") {
                vals_temp <- as.numeric(vals[c("mean", "lb_95", "ub_95")])
                fitt_dist[i] <- paste0("r", dist, "(1000,", round(vals_temp[2], 10), ", ", round(vals_temp[3], 10),  ",", vals_temp[1], ")")
            } else {
                fitvals <- optim(init, fitdist_ci, data = as.numeric(vals[c("mean", "lb_95", "ub_95")]), dist = eval(dist_q))$par
                fitt_dist[i] <- paste0("r", dist, "(1000,", round(fitvals[1], 10), ", ", round(fitvals[2], 10), ")")
            }
        }
    }
    cost_df$fitt_dist <- fitt_dist
    cost_df
}


#' Generate and output samples for economic and risk data
#'
#' This function processes raw economic and risk data, fits distributions to the data, 
#' generates samples, and outputs the results to a specified file.
#'
#' @param econ_raw A data frame containing raw economic data with columns for age group, outcome, reference QALYs, mean QALYs, lower bound QALYs, and upper bound QALYs.
#' @param risks_raw A data frame containing raw risk data.
#' @param file.name A string specifying the name of the output file.
#'
#' @return A list containing sampled outcomes for risks, QALYs, and costs.
#' 
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Fits triangular distributions to the risk, QALY, and cost data.
#'   \item Generates samples for each outcome using the fitted distributions.
#'   \item Outputs the sampled outcomes to the specified file.
#' }
#'
#' @examples
#' \dontrun{
#' econ_raw <- read.csv("economic_data.csv")
#' risks_raw <- read.csv("risk_data.csv")
#' output_samples(econ_raw, risks_raw, "output_samples.csv")
#' }
#'
#' @export
output_samples <- function(econ_raw, risks_raw, file.name) {

    require(EnvStats)

    risks_uncert <- add_fitted_dist(risks_raw, "tri")
    qalys_raw <- econ_raw %>% select(age_group, outcome, ref_qaly, mean_qaly, lb_95_qaly, ub_95_qaly) %>%
        rename_with(~ stringr::str_remove(., "_qaly"))
    costs_raw <- econ_raw %>% select(age_group, outcome, ref_costs, mean_costs, lb_95_costs, ub_95_costs) %>%
        rename_with(~ stringr::str_remove(., "_costs"))

    costs_uncert <- add_fitted_dist(costs_raw, "tri")
    qalys_uncert <- add_fitted_dist(qalys_raw, "tri")

    run_samples_outcomes <- function(df_uncert, type = "risk") {
        outcomes <- df_uncert$outcome %>% unique
        risk_samples <- list()
        for (j in 1:length(outcomes)) {
            dist_sample <- matrix(0, 1000, 25)
            dists <- df_uncert %>% dplyr::filter(outcome == outcomes[j]) %>% select(mean, fitt_dist)

            for (i in 1:25) {
                if (!is.na(dists$mean[i])) {
                    if (dists$fitt_dist[i] == "NA" | dists$mean[i] == 0) {
                        dist_sample[, i] <- rep(dists$mean[i], 1000)
                    } else {
                        dist_sample[, i] <- eval(parse(text = dists$fitt_dist[i]))
                    }
                }
            }
            risk_samples[[paste0(outcomes[j])]] <- dist_sample
        }
        risk_samples
    }

    risk_samples <- run_samples_outcomes(risks_uncert, "risks")
    qaly_samples <- run_samples_outcomes(qalys_uncert, "qalys")
    costs_samples <- run_samples_outcomes(costs_uncert, "costs")

    samples_outcomes <- list(risks = risk_samples, qalys = qaly_samples, costs = costs_samples)

}

#' @title Plot Economics Data
#' @description This function generates and saves plots for economic data, risks data, and very high-risk data from the given object.
#' @param object An object containing economic data frames (`econ_df`, `risks_df`, `risks_vhr_df`) and an economic name (`econ_name`).
#' @details The function performs the following steps:
#' \itemize{
#'   \item Checks if the directory for saving plots exists, and creates it if it does not.
#'   \item Processes the economic data frame (`econ_df`) to create a combined data frame for QALY loss and Costs.
#'   \item Generates and saves a plot for the economic data with uncertainty ribbons and lines.
#'   \item Generates and saves a plot for the risks data with uncertainty ribbons and lines.
#'   \item Generates and saves a plot for the very high-risk data with lines.
#' }
#' @import ggplot2
#' @import dplyr
#' @import here
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_ribbon geom_line facet_grid scale_y_continuous labs theme_bw ggsave
#' @importFrom dplyr bind_rows rename_with mutate filter
#' @importFrom here here
#' @export
plot_economics <- function(object) {
    
    if (!dir.exists(here::here("outputs", object@econ_name, "figs"))) {
        # create the folder if it doesn't exist
        dir.create(here::here("outputs", object@econ_name, "figs"))
    }

    age_groups <- object@econ_df$age_group %>% unique

    econ_df <- bind_rows(
        object@econ_df[, 1:6] %>%
            rename_with(~gsub("_qaly$", "", .), c("ref_qaly", "mean_qaly", "lb_95_qaly", "ub_95_qaly")) %>%
            mutate(econ_metric = "QALY loss"),
        object@econ_df[, c(1:2, 7:10)] %>%
            rename_with(~gsub("_costs$", "", .), c("ref_costs", "mean_costs", "lb_95_costs", "ub_95_costs")) %>%
            mutate(econ_metric = "Costs (GBP)")
    )


    econ_df %>% filter(mean > 0) %>%
        mutate(age_group = factor(age_group, levels = age_groups)) %>% 
        ggplot() + 
            geom_ribbon(aes(x = age_group, ymin = lb_95, ymax = ub_95,
                fill = outcome, group = outcome), alpha = 0.5) + 
            geom_line(aes(x = age_group, y = mean, color = outcome, group = outcome), size = 2) + 
            facet_grid(vars(econ_metric), scales = "free_y") +
            scale_y_continuous(trans = "log10") + labs(x = "Age group", y = "Value") + 
            theme_bw()
    ggsave(here::here("outputs", object@econ_name, "figs", "economics_raw_uncertainty.png"))


    object@risks_df %>% 
        mutate(age_group = factor(age_group, levels = age_groups)) %>% 
        ggplot() + 
            geom_ribbon(aes(x = age_group, ymin = lb_95, ymax = ub_95,
                fill = outcome, group = outcome), alpha = 0.5) + 
            geom_line(aes(x = age_group, y = mean, color = outcome, group = outcome), size = 2) + 
            labs(x = "Age group", y = "Value") + 
            theme_bw() + facet_grid(vars(outcome), scales = "free_y")
    ggsave(here::here("outputs", object@econ_name, "figs", "risks_raw_uncertainty.png"))

    object@risks_vhr_df %>% 
        mutate(age_group = factor(age_group, levels = age_groups)) %>% 
        ggplot() + 
            geom_line(aes(x = age_group, y = mean, color = outcome, group = outcome), size = 2) + 
            labs(x = "Age group", y = "Value") + 
            theme_bw() + facet_grid(vars(outcome), scales = "free_y")
    ggsave(here::here("outputs",  object@econ_name, "figs", "risks_raw_vhr_uncertainty.png"))

}



#' @title Plot Economic Fit
#' 
#' @description This function generates and saves plots comparing fitted distributions to reference data for risks, QALYs, and costs across different age groups.
#' 
#' @param object An S4 object containing economic data and model outcomes.
#' 
#' @details The function performs the following steps:
#' \itemize{
#'   \item Checks if the directory for saving plots exists, and creates it if it doesn't.
#'   \item Extracts unique age groups from the economic data.
#'   \item Defines a helper function to extract mean and 95% uncertainty intervals for model outcomes.
#'   \item Extracts fitted model outcomes for risks, QALYs, and costs.
#'   \item Prepares data frames for risks, QALYs, and costs by merging reference data with fitted model outcomes.
#'   \item Defines a helper function to convert data frames to long format for plotting.
#'   \item Converts the prepared data frames to long format.
#'   \item Generates and saves plots comparing fitted model outcomes to reference data for risks, QALYs, and costs.
#' }
#' 
#' @return None. The function saves the generated plots to the specified directory.
#' 
#' @import dplyr
#' @import ggplot2
#' @import purrr
#' @import here
#' @importFrom stats quantile
#' @importFrom dplyr bind_cols
#' 
#' @examples
#' \dontrun{
#'   plot_economics_fit(my_economic_model)
#' }
plot_economics_fit <- function(object) {
    
    if (!dir.exists(here::here("outputs", object@econ_name, "figs"))) {
        # create the folder if it doesn't exist
        dir.create(here::here("outputs", object@econ_name, "figs"))
    }
    age_groups <- object@econ_df$age_group %>% unique

    extract_uncert_fit <- function(metric) {
        object@outcomes_vec %>% map_df(~
            bind_cols(
                object@samples_outcomes[[metric]][[.x]] %>% apply(2, mean),
                object@samples_outcomes[[metric]][[.x]] %>% apply(2, function(x) quantile(x, c(0.025, 0.975))) %>% t
            ) %>% setNames(c("mean_fit", "lb_95_fit", "ub_95_fit")) %>% mutate(outcome = .x, age_group = age_groups)
        )
    }
    model_fit_risks <- extract_uncert_fit("risks")
    model_fit_qalys <- extract_uncert_fit("qalys")
    model_fit_costs <- extract_uncert_fit("costs")

    risks_df <- object@risks_df %>% select(age_group, outcome, ref, mean, lb_95, ub_95)
    qaly_df <- object@econ_df %>% select(age_group, outcome, ref_qaly, mean_qaly, lb_95_qaly, ub_95_qaly) %>% 
        rename_with(~gsub("_qaly$", "", .), c("ref_qaly", "mean_qaly", "lb_95_qaly", "ub_95_qaly"))
    costs_df <- object@econ_df %>% select(age_group, outcome, ref_costs, mean_costs, lb_95_costs, ub_95_costs) %>% 
        rename_with(~gsub("_costs$", "", .), c("ref_costs", "mean_costs", "lb_95_costs", "ub_95_costs"))

    risk_df_short_risks <- risks_df %>% left_join(model_fit_risks, by = c("outcome", "age_group"))
    risk_df_short_qalys <- qaly_df %>% left_join(model_fit_qalys, by = c("outcome", "age_group"))
    risk_df_short_costs <- costs_df %>% left_join(model_fit_costs, by = c("outcome", "age_group"))

    make_long <- function(risk_df_short) {
        risk_df <- bind_rows(
            risk_df_short[, 1:6] %>%
                mutate(econ_metric = "Data"),
            risk_df_short[, c(1:3, 7:9)] %>%
                rename_with(~gsub("_fit$", "", .), c("mean_fit", "lb_95_fit", "ub_95_fit")) %>%
                mutate(econ_metric = "Fitted")
        )
    }

    risk_df_risks <- make_long(risk_df_short_risks)
    risk_df_qalys <- make_long(risk_df_short_qalys)
    risk_df_costs <- make_long(risk_df_short_costs)

    risk_df_risks %>% mutate(age_group = factor(age_group, levels = age_groups)) %>%
        ggplot() +             
            geom_linerange(aes(x = age_group, y = mean, ymin = lb_95, ymax = ub_95, color = econ_metric),
                size = 2, position = position_dodge(0.75)) +
            geom_point(aes(x = age_group, y = mean, color = econ_metric),
                size = 5, position = position_dodge(0.75)) + theme_bw() +
            facet_grid(vars(outcome), scales = "free_y")
      ggsave(here::here("outputs", object@econ_name, "figs", "compare_fit_risks.png"))
  

    risk_df_qalys %>% mutate(age_group = factor(age_group, levels = age_groups)) %>%
        ggplot() +             
            geom_linerange(aes(x = age_group, y = mean, ymin = lb_95, ymax = ub_95, color = econ_metric),
                size = 2, position = position_dodge(0.75)) +
            geom_point(aes(x = age_group, y = mean, color = econ_metric),
                size = 5, position = position_dodge(0.75)) + theme_bw() +
            facet_grid(vars(outcome), scales = "free_y")
      ggsave(here::here("outputs", object@econ_name, "figs", "compare_fit_qalys.png"))

    risk_df_costs %>% mutate(age_group = factor(age_group, levels = age_groups)) %>%
        ggplot() +             
            geom_linerange(aes(x = age_group, y = mean, ymin = lb_95, ymax = ub_95, color = econ_metric),
                size = 2, position = position_dodge(0.75)) +
            geom_point(aes(x = age_group, y = mean, color = econ_metric),
                size = 5, position = position_dodge(0.75)) + theme_bw() +
            facet_grid(vars(outcome), scales = "free_y")
      ggsave(here::here("outputs", object@econ_name, "figs", "compare_fit_costs.png"))

}


#' Convert Raw Outcomes to Risk
#'
#' This function converts raw outcome data to risk metrics based on population data and model cases sample mean.
#'
#' @param object An object containing UK population data.
#' @param outcomes_raw A data frame containing raw outcome data.
#' @param model_cases_sample_mean A data frame containing model cases sample mean.
#'
#' @return A data frame with calculated risk metrics for each age group and outcome.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts population and age group data from the input object.
#'   \item Defines helper functions to check the completeness of data and calculate risk per infection.
#'   \item Checks that each outcome has data for all age groups.
#'   \item Combines raw outcomes data with model cases data and population data.
#'   \item Calculates the risk per infection for each outcome and age group.
#' }
#'
#' @examples
#' \dontrun{
#'   result <- covert_raw_to_risk(object, outcomes_raw, model_cases_sample_mean)
#'   print(result)
#' }
#' 
#' @import dplyr
#' @import tidyr
#' @export
covert_raw_to_risk <- function(object, outcomes_raw, model_cases_sample_mean) {

    population <- object@uk_data$population
    pop_age <- object@uk_data$populationAgeGroup
    outcomes <- unique(outcomes_raw$outcome)
    age_groups <- unique(outcomes_raw$age_group)

    ## Functions used
    # Create function to get the population size of each age group
    check_n <- function(df, outcomes) {
        df %>% filter(outcome == outcomes) %>% nrow == length(age_groups)
    }

    calc_outcome_per_inf <- function(df) {
        df %>% mutate(risk = 
            case_when(
                metric == "inc_inf"~(value / 100000),
                metric == "inc"~(value / 100000 * population_age) / incidence,
                metric == "n"~value / incidence,
            )
        )
    }

    # Create function to get the population size of each age group
    df_pop <- data.frame(
        age_group = age_groups,
        population_age = pop_age
    )

    ## Check that each metric values has each value in age_groups
    for (vals in outcomes){
        if(!check_n(outcomes_raw, vals)) {
            print(paste0("Error: not all age groups are represented in the data for outcome, ", vals))
        }
    }

    # Combine the raw outcomes data with the model cases data and population data
    outcomes_raw_full <- outcomes_raw %>%
        pivot_longer(c(mean:ub_95), names_to = "stat_summary",
            values_to = "value") %>% 
        left_join(rename(model_cases_sample_mean,
            incidence = value), by = "age_group" ) %>%
        left_join(df_pop, by = "age_group")

    # Calculate the risk per infection
    outcomes_risks_full <- outcomes_raw_full %>% calc_outcome_per_inf %>%
        select(!c(value, incidence)) %>% tidyr::pivot_wider(names_from = stat_summary, values_from = risk) %>%
        mutate(metric = "risk")
    outcomes_risks_full
}