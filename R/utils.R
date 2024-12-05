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