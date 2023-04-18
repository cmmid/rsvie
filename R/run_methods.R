convert_to_outcomes <- function(object) {
    S <- object@S
    outcomes_vec <- object@outcomes_vec %>% unique
    discount_rate <- object@econ_par$discount_rate

    VHR_inci <- LR_inci <- map(1:S, ~ matrix(0, 521, 25))

    for (s in 1:S) {
        raw_inci_s <- object@raw_inci[[s]]
        for (t_w in 1:nrow(raw_inci_s$inci)) {
            if (object@eff_1season) {
                t_w_mod <- (t_w - 1) %% 52 + 1 + 104
            } else {
                t_w_mod <- t_w
            }
            inci_tw <- raw_inci_s$inci[t_w_mod, ]
            LR_inci[[s]][t_w, ] <- (1:25 %>% map(~sum(inci_tw[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist)
            VHR_inci[[s]][t_w, ] <- (1:25 %>% map(~sum(inci_tw[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist)
        }
    } 

    # Get the incidence per age group and week for each outome
    outcome_per_sample <- list()
    for (s in 1:S) {
        outcome_per_sample[[s]] <- outcomes_vec %>%
        map_df(~ data.frame(
            age_group = rep(1:25, each = 52),
            outcome = .x,
            week_no = rep(1:52, 25),
            cases = c(t(t(LR_inci[[s]][1:52, ]) * object@samples_outcomes$risks[[.x]][s, ])), 
            cases_VHR = c(t(t(VHR_inci[[s]][1:52, ]) * object@samples_outcomes_VHR$risks[[.x]][s, ])),
            s = s
        ))
    }
    outcomes_week_age <- outcome_per_sample %>% bind_rows %>% mutate(cases_total = cases + cases_VHR)

    # Get the total cost and qaly loss per outcome over 10 years

    get_discounted_qalys <- function(s, outcome) {
        (((LR_inci[[s]] + VHR_inci[[s]]) * (object@samples_outcomes$risks[[outcome]][s, ] * object@samples_outcomes$qalys[[outcome]][s, ])) %>%
        apply(1, sum) * exp(-(1:521 - 1) * discount_rate / 52.0)) %>% sum
    }

    get_discounted_costs <- function(s, outcome) {
        (((LR_inci[[s]] + VHR_inci[[s]]) * (object@samples_outcomes$risks[[outcome]][s, ] * object@samples_outcomes$costs[[outcome]][s, ])) %>%
        apply(1, sum) * exp(-(1:521 - 1) * discount_rate / 52.0)) %>% sum
    }

    dis_qaly_df <- 1:S %>% map_df(~
        data.frame(
            outcomes = outcomes_vec,
            qaly = lapply(outcomes_vec, function(x) get_discounted_qalys(.x, x)) %>% unlist,
            s = .x
        )
    )
    dis_qaly_df <- dis_qaly_df %>% left_join(dis_qaly_df %>% group_by(s) %>% summarise(total = sum(qaly)), 
        by = "s")

    dis_cost_df <- 1:S %>% map_df(~
        data.frame(
            outcomes = outcomes_vec,
            cost = lapply(outcomes_vec, function(x) get_discounted_costs(.x, x)) %>% unlist,
            s = .x
        )
    )
    dis_cost_df <- dis_cost_df %>% left_join(dis_cost_df %>% group_by(s) %>%
        summarise(total = sum(cost)), 
        by = "s")

    outcomes <- list(
        outcomes = outcomes_week_age,
        costs = dis_cost_df,
        qaly = dis_qaly_df,
        doses = object@raw_inci[[1]]$doses
    )
    outcomes
}

plot_outcomes <- function(object) {
    object_get <- load(file = here::here("outputs",
        object@prog_name, "run_outputs.RData"))
    object <- get(object_get)
    bind_rows(
        object@outcomes$costs %>% select(s, total) %>% mutate(costs = "costs"),
        object@outcomes$qaly %>% select(s, total) %>% mutate(costs = "qaly")
    ) %>%
        ggplot() +
        geom_boxplot(aes(x = s, y = total, fill = costs)) +
        scale_y_continuous(trans = "log10")

    ggsave(here::here("outputs", object@prog_name, "figs", "plot_qaly_cost.png"))

}


get_model_sample <- function(object, seed) {
    model_calendar_sample <- object@model_calendar


    if (object@immune_profile$vhr$product != "none") {
        vhr_eff <- if_else(length(object@immune_profile$vhr$a) > 1, object@immune_profile$vhr$a[seed], object@immune_profile$vhr$a[1])
        model_calendar_sample$mAB_VHR <- object@model_calendar$mAB_VHR * vhr_eff
    }
    if (object@immune_profile$mass$product != "none") {
        mass_eff <- if_else(length(object@immune_profile$mass$a) > 1, object@immune_profile$mass$a[seed], object@immune_profile$mass$a[1])
        if(object@immune_profile$mass$product == "mab") {
            model_calendar_sample$mAB_LR <- object@model_calendar$mAB_LR * mass_eff
        } else if (object@immune_profile$mass$product == "lav") {
            model_calendar_sample$LAV_LR <- object@model_calendar$LAV_LR * mass_eff
        } else if (object@immune_profile$mass$product == "mat") {
            model_calendar_sample$mat_LR[, 2:25] <- object@model_calendar$mat_LR[, 2:25] * mass_eff
            matt_eff <- if_else(length(object@immune_profile$mass$a_mat) > 1, object@immune_profile$mass$a_mat[seed], object@immune_profile$mass$a_mat[1])
            model_calendar_sample$mat_LR[, 1] <- object@model_calendar$mat_LR[, 1] * matt_eff
        }
    }
    model_calendar_sample
}

get_efficacies_sample <- function(object, seed) {
    efficacy_sample_list <- list(
        mab_vhr = 0,
        mab_mass = 0,
        lav_mass = 0,
        mat_mass = 0
    )
    for (eff_name in names(efficacy_sample_list) ) {
        if (length(object@efficacies[[eff_name]]) > 1) {
            efficacy_sample_list[[eff_name]] <- object@efficacies[[eff_name]][seed]
        } else {
            efficacy_sample_list[[eff_name]] <- object@efficacies[[eff_name]]
        }
    }
    efficacy_sample_list
}
