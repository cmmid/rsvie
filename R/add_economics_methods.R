## Add uncertainty to costs and qalys
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

#' @title Convert mean and ci into samples from a fitted distribution
#'
#' @description Convert mean and ci into samples from a fitted distribution.
#' @param econ_raw a dataframe with the economic parameters per age group.
#' @param risks_raw a dataframe with the risk of outcomes per age group.
#'
#' @import stringr
#' @import dplyr
#' @import ggplot2
#' @import purrr
#' @importFrom EnvStats rtri
#' @export
output_samples <- function(econ_raw, risks_raw, file.name) {

    risks_uncert <- add_fitted_dist(risks_raw, "tri")
    qalys_raw <- econ_raw %>% select(age_group, outcome, ref_qaly, mean_qaly, lb_95_qaly, ub_95_qaly) %>%
        rename_with(~ str_remove(., "_qaly"))
    costs_raw <- econ_raw %>% select(age_group, outcome, ref_costs, mean_costs, lb_95_costs, ub_95_costs) %>%
        rename_with(~ str_remove(., "_costs"))

    costs_uncert <- add_fitted_dist(costs_raw, "tri")
    qalys_uncert <- add_fitted_dist(qalys_raw, "tri")

    # inner workings
   # write.csv(risks_uncert, file = here::here("data", "inner_workings", paste0("outcome_risks_", file.name,"_dist.csv")))
   # write.csv(qalys_uncert, file = here::here("data", "inner_workings", "outcome_qalys_dist.csv"))
   # write.csv(costs_uncert, file = here::here("data", "inner_workings", "outcome_costs_dist.csv"))

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
            mutate(econ_metric = "Costs (£)")
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
    ggsave(here::here("outputs", object@econ_name, "figs", "risks_raw_vhr_uncertainty.png"))

}

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
            incidence = value) ) %>%
        left_join(df_pop)

    # Calculate the risk per infection
    outcomes_risks_full <- outcomes_raw_full %>% calc_outcome_per_inf %>%
        select(!c(value, incidence)) %>% tidyr::pivot_wider(names_from = stat_summary, values_from = risk) %>%
        mutate(metric = "risk")
    outcomes_risks_full
}

plot_risks <- function(outcomes_risks_full) {
    outcomes_risks_full %>% mutate(age_group = factor(age_group, levels = age_groups)) %>% 
        mutate(outcome = factor(outcome, levels = outcomes)) %>% 
        ggplot() + 
            geom_ribbon(aes(x = age_group, ymin = lb_95, ymax = ub_95, group = outcome), alpha = 0.5) +
            geom_line(aes(x = age_group, y = mean, group = outcome)) +
            facet_grid(vars(outcome), scales = "free_y") + 
            labs(y = "Risk per RSV infection", x = "Age group") + theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(here::here("outputs", "figs", "risks_per_infection.png"), width = 10, height = 10)
}