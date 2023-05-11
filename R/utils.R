#' @title function to conver the incidence matrix from a run to a data.frame
#' 
#' @param object an RSVProgramme object
#' @importFrom tidyr pivot_longer
#' @return none
#' @export
checkout_incidence <- function(object, rerun = TRUE) {
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

#' @title function to conver the incidence matrix from a state run to a data.frame
#' 
#' @param object an RSVProgramme object
#' @importFrom tidyr pivot_longer
#' @import data.table
#' @return none
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