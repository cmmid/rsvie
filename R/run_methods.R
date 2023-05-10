#' @export
convert_to_outcomes <- function(object) {

    # Get some data
    S <- object@S
    outcomes_vec <- object@outcomes_vec %>% unique
    discount_rate <- object@econ_par$discount_rate
    if (is.null(object@immune_profile$disease_eff)) {
        stop("There is no disease specific efficacy at `immune_profile$disease_eff`")
    } else {
        disease_eff_probs <- object@immune_profile$disease_eff %>% 
            mutate(across(infection:death) / infection)
    }
    # Get states to figure out proportion that are in V at each time
    object_states <- checkout_states(object)
    epi_states <- c( "M", map(0:3, ~paste0(c("S", "E", "A", "I",  "R"), .x)) %>% unlist)
    df_estimate_novac <- object_states[model_type == "inter" & state_name %in% epi_states, .(notin_V = sum(incidence)), 
            keyby = .(t, age_group, risk_group, model_type)]

    # Lots of functions to help do this
    calculate_prop_inV <- function(df_estimate_novac, object_states, product) {
        df_estimate_vac_prod <- object_states[model_type == "inter" &  grepl(as.name(product), state_name), .(in_V = sum(incidence)), 
            keyby = .(t, age_group, risk_group, model_type)]
        df_estimated_both_prod <- left_join(df_estimate_novac, df_estimate_vac_prod,
            by = join_by(t, age_group, risk_group, model_type) ) %>%
            mutate(prop_in_V = case_when(
                    notin_V + in_V == 0~0,
                    notin_V + in_V != 0~in_V / (in_V + notin_V)
                )
            )
        df_estimated_both_prod %>% filter((t > 2 * 365)) %>%
            .[1:(7 * 3 * 25 * 521), ] %>%
            mutate(week = rep(1:(521), each = 7 * 3 * 25)) %>%
            select(!c(t, model_type)) %>% rename(t = week)
    }

    add_prop_disease <- function(df_inci, dis_string, product, disease_eff_probs) {
        prop <- disease_eff_probs[disease_eff_probs$product == product, ][[dis_string]]
        colname <- paste0("incidence_protected_", dis_string )

        df_inci %>% 
            mutate(!!colname :=
                case_when(
                    notin_V == 0~0,
                    notin_V != 0~pmin(1, prop_in_V * (prop - 1) / (1 - prop_in_V) )
                )
            )
    }

    add_healthout_protected <- function(df_estimated_both_in, disease_eff_probs, product) {
        df_estimated_both_in %>% 
            add_prop_disease("symptomatic", product, disease_eff_probs) %>%
            add_prop_disease("gp", product, disease_eff_probs) %>%
            add_prop_disease("hosp", product, disease_eff_probs) %>%
            add_prop_disease("a_e", product, disease_eff_probs) %>%
            add_prop_disease("icu", product, disease_eff_probs) %>%
            add_prop_disease("death", product, disease_eff_probs) %>% 
            summarise(across(incidence_protected_symptomatic:incidence_protected_death, ~mean(.x)), 
                .by = c(age_group, risk_group, t))
    }

    # combination of these
    df_estimated_both_mat <- calculate_prop_inV(df_estimate_novac, object_states, "mat") %>%
        add_healthout_protected(disease_eff_probs, "mat") 
    df_estimated_both_lav <- calculate_prop_inV(df_estimate_novac, object_states, "lav") %>%
        add_healthout_protected(disease_eff_probs, "lav")    
    df_estimated_both_mab <- df_estimated_both_mab <- calculate_prop_inV(df_estimate_novac, object_states, "mab") %>%
        add_healthout_protected(disease_eff_probs, "mab") 

    total_mat_prot <- 
        as.matrix( df_estimated_both_mat %>% select(incidence_protected_symptomatic:incidence_protected_death)) +
        as.matrix( df_estimated_both_lav %>% select(incidence_protected_symptomatic:incidence_protected_death)) +
        as.matrix( df_estimated_both_mab %>% select(incidence_protected_symptomatic:incidence_protected_death))
    df_estimated_both <- df_estimated_both_mat %>% select(age_group:t) %>% bind_cols(total_mat_prot)

    # Find incidence per age group
    VHR_inci <- LR_inci <- map(1:S, ~ matrix(0, 521, 25))
    for (s in 1:S) {
        raw_inci_s <- object@raw_inci[[s]]
        for (t_w in 1:nrow(raw_inci_s$inci)) {
            inci_tw <- raw_inci_s$inci[t_w, ]
            LR_inci[[s]][t_w, ] <- (1:25 %>% map(~sum(inci_tw[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist)
            VHR_inci[[s]][t_w, ] <- (1:25 %>% map(~sum(inci_tw[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist)
        }
    }

    # Get proportion fo incidence protected from severe outcomes
    VHR_dis_prot <- df_estimated_both %>% filter(risk_group == "VHR") %>% as.data.table
    LR_dis_prot <- df_estimated_both %>% filter(risk_group == "LR")  %>% as.data.table

    LR_prot <- outcomes_vec %>% map( ~LR_dis_prot[[paste0("incidence_protected_", .x)]] %>% 
                    matrix(521, 25, byrow = TRUE)) %>% setNames(outcomes_vec)

    VHR_prot <- outcomes_vec %>% map( ~VHR_dis_prot[[paste0("incidence_protected_", .x)]] %>% 
                    matrix(521, 25, byrow = TRUE)) %>% setNames(outcomes_vec)

    # Get the risk, incidence, costs and qalys
    outcome_per_sample <- list()
    for (s in 1:S) {
        outcome_per_sample[[s]] <- outcomes_vec %>%
        map_df( 
            function(x) {
                inci_risk_hc <- t(LR_inci[[s]] * (1 - LR_prot[[x]])) 
                inci_risk_hc_vhr <- t(VHR_inci[[s]] * (1 - VHR_prot[[x]])) 

                data.frame(
                    age_group = rep(1:25, each = 521),
                    outcome = x,
                    week_no = rep(1:521, 25),
                    cases = c(t(inci_risk_hc * object@samples_outcomes$risks[[x]][s, ])) + c(t(inci_risk_hc_vhr * object@samples_outcomes_VHR$risks[[x]][s, ])),
                    qaly = c(t(inci_risk_hc * object@samples_outcomes$risks[[x]][s, ] * object@samples_outcomes$qaly[[x]][s, ])) + c(t(inci_risk_hc_vhr * object@samples_outcomes_VHR$risks[[x]][s, ] * object@samples_outcomes$qaly[[x]][s, ])),
                    cost = c(t(inci_risk_hc * object@samples_outcomes$risks[[x]][s, ] * object@samples_outcomes$cost[[x]][s, ])) + c(t(inci_risk_hc_vhr * object@samples_outcomes_VHR$risks[[x]][s, ] * object@samples_outcomes$cost[[x]][s, ])),
                    s = s
                )
            }
        )
    } 

    outcomes_week_age <- outcome_per_sample %>% bind_rows 
    outcomes_week_age_trim <- outcomes_week_age %>% filter(week_no %in% (c(1:52) + 104))

    dis_qaly_df <- outcomes_week_age %>% mutate(qaly = qaly * exp(-(week_no * discount_rate / 52.0))) %>% group_by(outcome, age_group, s) %>% summarise(qaly = sum(qaly))
    dis_cost_df <- outcomes_week_age %>% mutate(cost = cost * exp(-(week_no * discount_rate / 52.0))) %>% group_by(outcome, age_group, s) %>% summarise(cost = sum(cost))

    dis_qaly_df <- dis_qaly_df %>% left_join(dis_qaly_df %>% group_by(s) %>% summarise(total = sum(qaly)),
        by = "s")
    dis_cost_df <- dis_cost_df %>% left_join(dis_cost_df %>% group_by(s) %>% summarise(total = sum(cost)),
        by = "s")

    outcomes <- list(
        outcomes = outcomes_week_age_trim,
        costs = dis_cost_df,
        qaly = dis_qaly_df,
        doses = object@raw_inci[[1]]$doses
    )
    outcomes
}

#' @export
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

#' @export
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
