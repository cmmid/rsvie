#' @title function to conver the incidence matrix from a run to a data.frame
#' 
#' @param object an RSVProgramme object
#' @importFrom tidyr pivot_longer
#' @return none
#' @title Checkout Incidence
#' @description This function processes the incidence data from an object.
#' @param object An object containing the raw incidence data.
#' @param rerun A logical value indicating whether to rerun the `rsvie::run` function on the object. Default is TRUE.
#' @return A data frame containing the processed incidence data with additional columns for time, age group, social group, and risk group.
#' @examples
#' \dontrun{
#'   result <- checkout_incidence(my_object)
#' }
#' @export
checkout_incidence <- function(object, rerun = FALSE) {
    if (rerun) {
        object <- rsvie::run(object)    
    }

    inci_df <- object@raw_inci[[1]]$inci %>% as.data.frame %>% mutate(t = 1:521) %>%
        pivot_longer(!t, names_to = "state", values_to = "incidence") %>%
        mutate(
            age_group = rep(rep(1:25, each = 9), 521),
            social_group = rep(rep(c("p", "c", "n"), each = 3), 521 * 25),
            risk_group = rep(c("VHR", "HR", "LR"), 3 * 25 * 521)
        ) 
    inci_df
}


#' @title Checkout States
#' @description This function processes the states of a given object and returns a data table with detailed information about each state.
#' @param object An object containing the states to be processed.
#' @return A data table with columns for time, incidence, age group, risk group, social group, model type, and state name.
#' @details The function generates a list of state names and categorizes them into risk groups, social groups, and model types. It then runs the state function on the input object and processes the resulting states into a data frame. The data frame is then transformed into a long format and additional columns are added to provide detailed information about each state.
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom data.table as.data.table
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @export
checkout_states <- function(object) {
    inner <- c(
        "M", 
        map(0:3, ~paste0(c("S", "E", "A", "I",  "R"), .x)) %>% unlist,
        "V_mat1", "V_mat2", "V_mat3",
        "V0_mab1", "V0_mab2", "V0_mab3",
        "V1_mab1", "V1_mab2", "V1_mab3",
        "V2_mab1", "V2_mab2", "V2_mab3",
        "V0_lav1",  "V0_lav2", "V0_lav3",
        "V1_lav1",  "V1_lav2", "V1_lav3",
        "V2_lav1",  "V2_lav2", "V2_lav3",
        "V3_lav1",  "V3_lav2", "V3_lav3"
    )
    state_names <- c(rep(
        inner, 18), paste0("sum", c(1:23)))

    L <- state_names %>% length
    SL <- inner %>% length

    risk_groups <- c(rep(rep(c("VHR", "HR", "LR"), each = SL), 6), rep("NA", 23))
    social_groups <- c(rep(rep(c("p", "c", "n"), each = SL * 3), 2), rep("NA", 23))
    model_types <- c(rep(c("no_inter", "inter"), each = SL * 9), rep("NA", 23))
    object <- 
    states <- run_state(object)    
    df_states <- states[[1]] %>% as.data.frame %>% mutate(t = 1:4745) %>%
        pivot_longer(!t, names_to = "time", values_to = "incidence") %>%
        mutate(
            age_group = rep(rep(1:25, each = L), 4745),
            risk_group = rep(risk_groups, 25 * 4745),
            social_group = rep(social_groups, 25 * 4745),
            model_type = rep(model_types, 25 * 4745),
            state_name = rep(state_names, 25 * 4745)
        )
    df_states %>% as.data.table
}



#' Compare Interventions
#'
#' This function compares the outcomes of a base intervention with multiple other interventions.
#'
#' @param obj_base An object containing the outcomes of the base intervention.
#' @param obj_inter A list of objects, each containing the outcomes of different interventions to be compared with the base intervention.
#'
#' @return A data frame summarizing the total cases, cases in very high-risk (VHR) groups, and total cases for each age group and outcome, for both the base and other interventions.
#'
#' @import dplyr
#' @import purrr
#' @export
compare_interventions <- function(obj_base, obj_inter) {

    base_sum <- obj_base@outcomes$outcomes %>% 
        group_by(age_group, outcome, s) %>% 
        summarise(cases = sum(cases), cases_VHR = sum(cases_VHR), cases_total = sum(cases_total)) %>%
        mutate(type = obj_base@prog_name)

    base_sum1 <- obj_inter %>% map_df(~
            .x@outcomes$outcomes %>% 
            group_by(age_group, outcome, s) %>% 
            summarise(cases = sum(cases), cases_VHR = sum(cases_VHR), cases_total = sum(cases_total)) %>%
            mutate(type = .x@prog_name)
    )

    base_compare <- base_sum %>% bind_rows(base_sum1)
    base_compare
}



#' Summarise Outcomes
#'
#' This function takes an object `x` and summarises the outcomes by grouping 
#' by `s`, `outcome`, and `age_group`, then calculates the total number of cases.
#'
#' @param x An object containing the outcomes data.
#' @return A data frame with the total number of cases grouped by `s`, `outcome`, and `age_group`.
#' @import dplyr
#' @export
summarise_outcomes <- function(x) {
    x@outcomes$outcomes %>% group_by(s, outcome, age_group) %>% summarise(cases_total = sum(cases) )
}


#' @title Calculate Averted Cases Data Frame
#' @description This function calculates the number of averted cases for different interventions compared to a base scenario.
#' @param base A data frame containing the base scenario data with columns `s`, `outcome`, `age_group`, and `cases_total`.
#' @param interventions A data frame containing the intervention scenario data with columns `s`, `outcome`, `age_group`, and `cases_total`.
#' @return A data frame with the calculated averted cases and proportions, with recoded age groups and outcomes.
#' @details The function performs the following steps:
#' \itemize{
#'   \item Joins the base and intervention data frames on `s`, `outcome`, and `age_group`.
#'   \item Calculates the total number of cases averted and the proportion of cases averted.
#'   \item Recodes the age groups and outcomes to more descriptive labels.
#'   \item Converts the recoded age groups and outcomes to factors with specific levels.
#' }
#' @examples
#' \dontrun{
#' base <- data.frame(s = 1:5, outcome = rep("symptomatic", 5), age_group = 1:5, cases_total = c(100, 200, 150, 120, 130))
#' interventions <- data.frame(s = 1:5, outcome = rep("symptomatic", 5), age_group = 1:5, cases_total = c(90, 180, 140, 110, 120))
#' get_averted_df(base, interventions)
#' }
#' @import dplyr
#' @export
get_averted_df <- function(base, interventions) {
    recode_age <- c("1" = "<1 month", "2" = "1 month", "3" = "2 month", "4" = "3 month", "5" = "4 month",
    "6" = "5 month", "7" = "6 month", "8" = "7 month", "9" = "8 month", "10" = "9 month", "11" = "10 month",
    "12" = "11 month", "13" = "1 year", "14" = "2 year", "15" = "3 years", "16" = "4 years", 
    "17" = "5-9 years", "18" = "10-14 years", "19" = "15-24 years", "20" = "25-34 years", "21" = "35-44 years",
    "22" = "45-54 years", "23" = "55-64 years", "24" = "65-74 years", "25" = "75+ years")
    relabel_outcomes <- c("symptomatic" = "Symptomatic cases", "gp" = "GP consultations", 
        "hosp" = "Hospital cases", "icu" = "ICU admissions", "a_e" = "A+E visits", "death" = "Deaths")
    RSV_impact <- interventions %>%
        left_join(base %>% rename(case_total_base = cases_total), by = c("s", "outcome", "age_group")) %>%
        mutate(total_case_averted = case_total_base - cases_total) %>% mutate(prop_cases_averted = (case_total_base - cases_total) / case_total_base)  %>% mutate(age_group = recode(age_group, !!!recode_age)) %>%
        mutate(age_group = factor(age_group, levels = unique(recode_age))) %>%
        mutate(outcome = recode(outcome, !!!relabel_outcomes)) %>%
        mutate(outcome = factor(outcome, levels = unique(relabel_outcomes)))
    RSV_impact
}