# Define the S4 class for the RSV intervention programme
#' @export
setClass(
    "RSVProgramme",
    representation(
        econ_name = "character",
        prog_name = "character",

        model = "refClass",
        uk_data = "list",
        model_par = "list",
        econ_par = "list",
        model_calendar = "list",
        dose_calendar = "list",
        immune_profile = "list",
        econ_df = "data.frame",
        risks_df = "data.frame",
        risks_vhr_df = "data.frame",
        outcomes_vec = "vector",
        post = "data.frame",
        seeds = "vector",
        S = "numeric",
        raw_inci = "list",
        wane_func_string = "character",
        wane_func_string_vhr = "character",
        outcomes = "list",
        samples_outcomes = "list",
        samples_outcomes_VHR = "list",
        efficacies = "vector",
        cov_mat = "numeric",
        full_output = "logical"
    )
)

#` The data used in the epidemic intervention model is as follows
#` * rsv_data_uk, demographic data from England and Wales
    # * numberDailyLiveBirths---number of daily live births in England and Wales
    # * population---total population
    # * ageGroupBoundary---age groups considered
#` data_inter_uk, data from England and Wales used in the intervention model
    #  * p_mat---proportion of each age group which are new mothers
    #  * u_p---proportion of each new mother age group
    #  * pVHR---proportion of each age groups which are considered very-high-risk
    #  * pHR---proportion of each age groups which are considered high-risk
    #  * pLR---proportion of each age groups which are considered low-risk
    #  * cnt_matrix_p---contact matrix showing mean number of total physical contacts between each age group
    #  * cnt_matrix_p_h---contact matrix showing mean number of household physical contacts between each age group
    #  * cnt_matrix_c---contact matrix showing mean number of total conversational contacts between each age group
    #  * cnt_matrix_c_h---contact matrix showing mean number of household conversational contacts between each age group
#` run_start, step number to start the model on (0)
#` run_burn, number of days to run the burn in model
#` run_full, number of days, after burn in, to run the model for


#' @title Initialise the RunInterventions class
#'
#' @param RunInterventions An empty RunInterventions class
#' @param uk_data_sum an .rda file with all the information needed to fit the UK model
#' @param run_burn Burn-in for the model (in days)
#' @param run_full Number of days to run the model for post burn-in
#' @return A classRunInterventions which is fully parameterised
#' @export
make_RunInterventions <- function(RunInterventions, uk_data_sum, run_burn, run_full) {
    classRunInterventions <- new(RunInterventions, uk_data_sum$numberDailyLiveBirths, uk_data_sum$population, uk_data_sum$ageGroupBoundary) # Calls class
    classRunInterventions$p_mat <- uk_data_sum$prop_mat
    classRunInterventions$pVHR <- uk_data_sum$pVHR
    classRunInterventions$pHR <- uk_data_sum$pHR
    classRunInterventions$pLR <- uk_data_sum$pLR

    classRunInterventions$cnt_matrix_p <- uk_data_sum$cnt_matrix_p
    classRunInterventions$cnt_matrix_p_h <- uk_data_sum$cnt_matrix_p_h
    classRunInterventions$pwp_p <- uk_data_sum$pwp_p
    classRunInterventions$pwn_p <- uk_data_sum$pwn_p
    classRunInterventions$nwp_p <- uk_data_sum$nwp_p
    classRunInterventions$nwn_p <- uk_data_sum$nwn_p

    classRunInterventions$cnt_matrix_c <- uk_data_sum$cnt_matrix_c
    classRunInterventions$cnt_matrix_c_h <- uk_data_sum$cnt_matrix_c_h
    classRunInterventions$pwp_c <- uk_data_sum$pwp_c
    classRunInterventions$pwn_c <- uk_data_sum$pwn_c
    classRunInterventions$nwp_c <- uk_data_sum$nwp_c
    classRunInterventions$nwn_c <- uk_data_sum$nwn_c

    classRunInterventions$u_p <- uk_data_sum$prop_mat[19:21] / sum(uk_data_sum$prop_mat[19:21])

    classRunInterventions$run_start <- 0
    classRunInterventions$run_burn <- run_burn
    classRunInterventions$run_full <- run_full
    classRunInterventions
 }


Rcpp::sourceCpp("./src/RunInterventions.cpp")

#' @title A constructor to make an RSVProgramme object 
#'
#' @param S an integer taking the number of Monte Carlo samples to run over
#' @return A RSVProgramme object
#' @export
make_rsv_programme <- function(S = 10) {
     # Uploads empty RunInterventions class
  #  Rcpp::sourceCpp("./src/RunInterventions.cpp")
    #Rcpp::loadModule("RunInterventionsModule", "rsvie")

    RunInterventionsModule <- Rcpp::Module("RunInterventionsModule", "rsvie")
    RunInterventions <- RunInterventionsModule$RunInterventions 

    load(file = system.file(package = "rsvie", "extdata", "data-raw", "uk_data_sum.RData"))
    load(file = system.file(package = "rsvie", "extdata", "data-raw", "posteriors.Rda"))
    seeds <- read.csv(file = system.file(package = "rsvie", "extdata", "data-raw", "seed_samples.csv"))[[1]] # loads uk_data_sum

    model_par <- list(burnin_yrs = 2, run_yrs = 10)
    econ_par <- list(time_horizon = 10, discount_rate = 0.035)

    model <- make_RunInterventions(RunInterventions, uk_data_sum,
        model_par$burnin_yrs * 365, model_par$run_yrs * 365)

    new("RSVProgramme", model = model, uk_data = uk_data_sum,
        model_par = model_par, econ_par = econ_par,
        post = post, seeds = seeds, S = S,
        full_output = FALSE,
        model_calendar = list(), dose_calendar = list(), immune_profile = list()
        )
}

#' @export
setGeneric("add_economics", function(object, econ_name, econ_raw, risks_raw, risks_vhr_raw) {
  standardGeneric("add_economics")
})

#' @title A setMethod for the `add_economics` function
#'
#' @param econ_name a string for the name gven the economic framework given
#' @param econ_raw a dataframe with the economic parameters per age group.
#' @param risks_raw a dataframe with the risk of outcomes per age group.
#' @param risks_vhr_raw a dataframe with the risk of outcomes per age group for very-high-risk infants
#' @return An RSVProgramme object
#' @export
setMethod("add_economics",  signature(object = "RSVProgramme"), 
    function(object, econ_name, econ_raw, risks_raw, risks_vhr_raw) {

        if (!dir.exists(here::here("outputs", econ_name))) {
            # create the folder if it doesn't exist
            dir.create(here::here("outputs", econ_name))
        }
        object@econ_name <- econ_name
        object@econ_df <- econ_raw
        object@risks_df <- risks_raw
        object@risks_vhr_df <- risks_vhr_raw
        object@outcomes_vec <- econ_raw$outcome %>% unique
        if(length(econ_raw$outcome %>% unique) == 0) {
            stop("ERROR: no outcomes detected.")
        }
        object@samples_outcomes <- output_samples(object@econ_df, object@risks_df, "")
        object@samples_outcomes_VHR <- output_samples(object@econ_df, object@risks_vhr_df, "VHR")
        object
    }
)

#' @export
setGeneric("add_programme", function(object, prog_name, cal, cal_vhr, immune_profile) {
  standardGeneric("add_programme")
})


#' @title A setMethod for the `add_programme` function
#'
#' @param prog_name an integer taking the number of Monte Carlo samples to run over
#' @param cal a dataframe with the implementation calendar for the population
#' @param cal_vhr dataframe with the implementation calendar for the very-high-risk population
#' @param immune_profile a list with information about the immune asumptions if the prophylatic in the programme
#' @return An RSVProgramme object
#' @export
setMethod("add_programme",  signature(object = "RSVProgramme"), 
    function(object, prog_name, cal, cal_vhr, immune_profile) {

        object@prog_name <- prog_name
        if (!is.null(immune_profile$full_output)) {
            object@full_output <- immune_profile$full_output
        }

        if (!dir.exists(here::here("outputs", "extra", prog_name))) {
            # create the folder if it doesn't exist
            dir.create(here::here("outputs", "extra", prog_name))
        }

        # These are all the programmes that are run
        all_dose <- list(
            mAB_VHR = matrix(0, 365, 25),
            mAB_HR = matrix(0, 365, 25),
            mAB_LR = matrix(0, 365, 25),
            LAV_HR = matrix(0, 365, 25),
            LAV_LR = matrix(0, 365, 25),
            mat_LR = matrix(0, 365, 25),
            mAB_VHR_custom_out = matrix(0, 365, 25),
            mAB_HR_custom_out = matrix(0, 365, 25),
            mAB_LR_custom_out = matrix(0, 365, 25)
        )
        all_sero <- all_dose

        # Extract the wane function in a string

        # Load relevant functions for getting string
        # Sort vhr
        daily_uptake_vhr <- calculate_daily_uptake(object, cal_vhr)
        daily_uptake <- calculate_daily_uptake(object, cal)
        all_both <- convert_mat_to_trans_cal(all_dose, all_sero, cal_vhr, daily_uptake_vhr, cal, daily_uptake, immune_profile)

        efficacy_estimates <- convert_efficacies(immune_profile)

        object@cov_mat <- get_coverage(immune_profile$mass, cal)
        object@model_calendar <- all_both[["sero"]]
        object@dose_calendar <- all_both[["dose"]]
        object@immune_profile <- immune_profile
        object@efficacies <- efficacy_estimates
        object
    }
)

#' @export
setGeneric("run", function(object, direct = FALSE, filename = NULL) {
  standardGeneric("run")
})

#' @title A setMethod for the `run` function
#'
#' @param object RSVProgramme object
#' @return An RSVProgramme object
#' @import furrr
#' @import future
#' @export
setMethod("run",  signature(object = "RSVProgramme"),

    function(object, direct = FALSE, filename = NULL) {
        cat("Running: Iterating through model simulations\n")
        future::plan(multisession, workers = 16)
        raw_inci <-
            map(seq_len(object@S),
                function(x) {

                cat(" Sample no: ", x, "    \r")
                model_calendar_eff <- get_model_sample(object, object@seeds[x])
                efficacies_sample <- get_efficacies_sample(object, object@seeds[x])
                efficacies_sample$direct <- direct

                object@model$Sample(
                    model_calendar_eff,
                    object@dose_calendar,
                    object@cov_mat,
                    efficacies_sample,
                    as.matrix(object@post)[object@seeds[x], ]
                )
                }
            )
        object@raw_inci <- raw_inci
        cat("Running: Calculating the outcomes\n")
        object@outcomes <- convert_to_outcomes(object)

        if (!dir.exists(here::here("outputs", "extra", object@prog_name))) {
            # create the folder if it doesn't exist
            dir.create(here::here("outputs", "extra", object@prog_name))
        }
        if (!is.null(filename)) {
            save(object, file = here::here("outputs", "extra", object@prog_name, "run_outputs.RData"))
        }
        object

    }
)

#' @export
setGeneric("run_state", function(object, S_i = 1, direct = FALSE) {
  standardGeneric("run_state")
})

#' @title A setMethod for the `run_state` function
#'
#' @param object RSVProgramme object
#' @return a matrix with all the values of the state variables in the dynamic transmission model
#' @export
setMethod("run_state",  signature(object = "RSVProgramme"),
    function(object, S_i = 1, direct = FALSE) {

        raw_state <-
            map(seq_len(S_i),
                function(x) {
                    model_calendar_eff <- get_model_sample(object, object@seeds[x])
                    efficacies_sample <- get_efficacies_sample(object, object@seeds[x])
                    efficacies_sample$direct <- direct
                    object@model$StatesValues(
                        model_calendar_eff,
                        object@dose_calendar,
                        object@cov_mat,
                        efficacies_sample,
                        as.matrix(object@post)[object@seeds[x], ]
                    )
                }
            )

        raw_state

    }
)
