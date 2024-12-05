#' @title Calculate Daily Uptake
#' @description This function calculates the daily uptake for different age groups based on the provided calibration data frame.
#' @param object An object containing the UK data with age group boundaries.
#' @param cal_df A data frame containing the start and end dates, and uptake types for each age group.
#' @return A matrix with 365 rows (days) and 25 columns (age groups) representing the daily uptake.
#' @details The function processes each age group, calculates the start and end positions within the year, and fills the uptake matrix accordingly. 
#' Different uptake types such as "linear", "birth_rate", "flu_oa", and "oa_2425" are handled.
#' @examples
#' \dontrun{
#' object <- your_object_with_uk_data
#' cal_df <- your_calibration_data_frame
#' daily_uptake_matrix <- calculate_daily_uptake(object, cal_df)
#' }
#' @importFrom lubridate dmy ymd year yday
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @export
calculate_daily_uptake <- function(object, cal_df) { 

    pop_size_days <- (object@uk_data$ageGroupBoundary * 365 ) %>% diff %>% round(0)

    cal_mat <- matrix(0, 365, 25)

    for (a in 1:25) {
        StartDate_format <- dmy(cal_df[a, ]$start_date)
        EndDate_format <- dmy(cal_df[a, ]$end_date)

        if (is.na(StartDate_format) | is.na(EndDate_format)) {
            next
        }        
        if (StartDate_format > EndDate_format) {
                stop("Error: Start date is after end date for age group: ", a)
        }
        # Sort out scaling for Start Date
        if ((StartDate_format <= ymd(paste0(year(StartDate_format), "-12-31"))) & (StartDate_format >= ymd(paste0(year(StartDate_format), "-07-01")))) {
            StartPos <- yday(StartDate_format) - 181
        }  else {
            StartPos <- yday(StartDate_format) - 181 + 365
        }

        if ((EndDate_format <= ymd(paste0(year(StartDate_format), "-12-31"))) & (EndDate_format >= ymd(paste0(year(StartDate_format), "-07-01")))) {
            EndPos <- yday(EndDate_format) - 181
        } else {
            EndPos <- yday(EndDate_format) + 365 - 181
        }

        if(!is.na(cal_df$uptake[a])) {
            if (cal_df$uptake[a] == "linear") {
                uptake <- 1 / pop_size_days[a]
            } else if (cal_df$uptake[a] == "birth_rate") {
                uptake <- 1 / 365
            } else if (cal_df$uptake[a] == "flu_oa") {
                flu_like <- read.table(system.file(package = "rsvie", "extdata", "uptakes_uk", "flu_oa.csv") , sep = ",") %>% as.numeric
                uptake <- get_daily_uptake(flu_like, StartPos)
            } else if (cal_df$uptake[a] == "oa_2425") {
                flu_like <- read.table(system.file(package = "rsvie", "extdata", "uptakes_uk", "oa_2425.csv") , sep = ",") %>% as.numeric
                uptake <- get_daily_uptake(flu_like, StartPos)
            }
            if (length(uptake) == 1) {
                if (StartPos < EndPos) {
                    cal_mat[c(StartPos:EndPos), a] <- uptake
                } else {
                    cal_mat[c(StartPos:365), a] <- uptake
                    cal_mat[c(1:EndPos), a] <- uptake
                }
            } else {
                cal_mat[, a] <- uptake
            }
        }

    }
    cal_mat
}


#' @title Calculate Daily Uptake from Weekly Data
#' @description This function calculates the daily uptake based on weekly uptake data.
#' @param up_week_raw A numeric vector containing the weekly uptake data.
#' @param start_time An integer representing the start time in days.
#' @return A numeric vector of length 365 representing the daily uptake.
#' @details The function first processes the weekly data to fill a weekly uptake vector, then calculates the daily uptake by interpolating the weekly data.
#' @examples
#' \dontrun{
#' up_week_raw <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' start_time <- 1
#' daily_uptake <- get_daily_uptake(up_week_raw, start_time)
#' }
#' @export
get_daily_uptake <- function(up_week_raw, start_time) {
    up_week <- vector("numeric", length = 52)
    up_day <- vector("numeric", length = 365)
    for (k in (start_time):(start_time + 20)) { 
        pos_week <- (k - 1) %% 52 + 1
        pos_ist <- (k - start_time) + 1
        up_week[pos_week] <- up_week_raw[pos_ist]
    }
            
    for (i in 1:365) {
        pos_day <- (i - 1 + start_time) %% 365 + 1
        up_day[pos_day] <- (up_week[((i / 7 + start_time) %% 52) + 1] - up_week[((i / 7 - 1 + start_time) %% 52) + 1]) / 7.0
        if (up_day[pos_day] < 0)
            up_day[pos_day] <- 0
    }
    up_day
}



#' @title Plot Dosing Calendar
#' @description This function generates and saves a plot of the dosing calendar for different age groups.
#' @param object An object containing the dosing calendar and other relevant data.
#' @details The function performs the following steps:
#' \itemize{
#'   \item Checks if the directory for saving the plot exists, and creates it if it does not.
#'   \item Extracts unique age groups from the economic data frame within the object.
#'   \item Updates the column names of the dose calendar matrices to match the age groups.
#'   \item Converts the dose calendar matrices into a data frame suitable for plotting.
#'   \item Generates a ggplot of the dosing calendar, with days of the year on the x-axis, age groups on the y-axis, and the proportion of doses as the fill color.
#'   \item Saves the plot as a PNG file in the specified directory.
#' }
#' @importFrom here here
#' @importFrom ggplot2 ggplot geom_tile facet_wrap scale_fill_gradient theme_bw theme labs ggsave
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map map_df
#' @importFrom magrittr %>%
#' @return None
#' @examples
#' \dontrun{
#' plot_calendar(my_object)
#' }
#' @export
plot_calendar <- function(object) {
    if (!dir.exists(here::here("outputs", "extra", object@prog_name, "figs"))) {
        dir.create(here::here("outputs", "extra", object@prog_name, "figs"))
    }
    age_groups <- object@econ_df$age_group %>% unique
    lapply(object@dose_calendar, function(x) colnames(x) <- age_groups)
    dose_calendar <- map(object@dose_calendar, ~{colnames(.x) <- age_groups; .x})
    names_dose <- dose_calendar %>% names

    df_calendar <- dose_calendar %>% map_df(
        ~as.data.frame(.x) %>% pivot_longer(everything(), names_to = "age_group", values_to = "prop") %>% 
            mutate(t = rep(1:365, each = 25), names(.x))
    ) %>% mutate(cal_name = rep(names_dose, each = 365 * 25)) %>%  mutate(age_group = factor(age_group, levels = age_groups))

    
    ggplot(df_calendar, aes(x = t, y = age_group, fill = prop)) + 
        geom_tile() + 
        facet_wrap(~cal_name, scales = "free") + 
        scale_fill_gradient(low = "blue", high = "red") + 
        theme_bw() + 
        theme(legend.position = "none") + 
        labs(x = "Day of year", y = "Age group")
    ggsave(here::here("outputs", object@prog_name, "figs", "dosing_calendar.png"))

}



#' Convert Matrices to Transmission Calendar
#'
#' This function converts dose and sero matrices into a transmission calendar based on the provided immune profile.
#'
#' @param all_dose A matrix representing all doses.
#' @param all_sero A matrix representing all sero values.
#' @param cal_vhr A calendar matrix for very high risk (VHR) individuals.
#' @param daily_uptake_vhr A matrix representing daily uptake for VHR individuals.
#' @param cal A calendar matrix for the general population.
#' @param daily_uptake A matrix representing daily uptake for the general population.
#' @param immune_profile A list containing immune profiles for VHR and general population.
#' 
#' @details
#' The function first creates a VHR calendar if the immune profile for VHR is not "none". It then calculates the sero and dose models using the `get_sero_dose_cals` helper function. The results are stored in the `all_sero` and `all_dose` matrices.
#' 
#' The `get_sero_dose_cals` function calculates the dose and sero models based on the provided calendar, daily uptake, and immune profile. It considers different sero delays and coverage values.
#' 
#' @return A list containing updated `dose` and `sero` matrices.
convert_mat_to_trans_cal <- function(all_dose, all_sero, cal_vhr, daily_uptake_vhr, cal, daily_uptake, immune_profile) { 


    # make vhr calendar
    ## vhr calendar must be a mab type, it can be palivuzumab or custom
    immune_profile_vhr <- immune_profile$vhr
    immune_profile_mass <- immune_profile$mass


    get_sero_dose_cals <- function(cal, daily_uptake, immune_profile) {
        dose_model <- matrix(0.0, 365, 25)
        sero_model <- matrix(0.0, 365, 25)
        ## make mass calendar 
        # get sero-delay


        if(immune_profile$sero_delay == "none") {
            sero_delay <-  c(1, rep(0, 29))
        } else if (immune_profile$sero_delay == "flu_like") {
            sero_delay <- dweibull(1:30, 2.42, 12.87)
        }
        # get coverage
        cov <- cal$cov
        A <- length(cov)
        # get efficacy
        eff <- immune_profile$a
        if(is.null(immune_profile[["gest_age_devilery_wks"]])) {
            for (t in 1:365) {
                for (a in 1:A) {
                    dose_model[t, a] <- daily_uptake[t, a] * cov[a]
                    for (k in 1:30) {
                        sero_model[(t + k - 1) %% 365 + 1, a] <- sero_model[(t + k - 1) %% 365 + 1, a] + daily_uptake[t, a] * cov[a] * sero_delay[k]
                    }
                }
            }
        } else {
            for (t in 1:365) {      
                for (a in 1:A) {
                    dose_model[t, a] <- daily_uptake[t, a]
                    for (k in 1:30) {
                        sero_model[(t + k - 1) %% 365 + 1, a] <- sero_model[(t + k - 1) %% 365 + 1, a] + daily_uptake[t, a] * sero_delay[k]
                    }
                }
            }
            for (t in 1:365) {
                gest_age_devilery_wks <- immune_profile$gest_age_devilery_wks
                numbers <- strsplit(gest_age_devilery_wks, " ")
                earliest_date <- 40 - as.numeric(numbers[[1]][1])
                latest_date <- 40 - as.numeric(numbers[[1]][2])
                eff_mat <- immune_profile$a_mat
                for(k in (latest_date * 7):(earliest_date * 7 - 1)) {
                    sero_model[(t + k - 1) %% 365 + 1, 1] <- sero_model[(t + k - 1) %% 365 + 1, 1] + dose_model[(t - 1) %% 365 + 1, 20] * 365 * dunif(k, latest_date * 7, earliest_date * 7)
                }
            }
            sero_model[sero_model > 1] <- 1
        }
        list(dose = dose_model, sero = sero_model)
    }

    if (immune_profile_vhr$product != "none") {
        cals_vhr <- get_sero_dose_cals(cal_vhr, daily_uptake_vhr, immune_profile_vhr)
        if(immune_profile_vhr$product == "mab" | immune_profile_vhr$product == "pal") {
            all_sero["mAB_VHR"] <- cals_vhr["sero"]
            all_dose["mAB_VHR"] <- cals_vhr["dose"]
        }
    }
    if (immune_profile_mass$product != "none") {
        cals <- get_sero_dose_cals(cal, daily_uptake, immune_profile_mass)
        if(immune_profile_mass$product == "mab") {
            all_sero["mAB_LR"] <- cals["sero"]
            all_dose["mAB_LR"] <- cals["dose"]
        } else if (immune_profile_mass$product == "lav") {
            all_sero["LAV_LR"] <- cals["sero"]
            all_dose["LAV_LR"] <- cals["dose"]
        } else if (immune_profile_mass$product == "mat") {
            all_sero["mat_LR"] <- cals["sero"]
            all_dose["mat_LR"] <- cals["dose"]
        }
    }

    list(dose = all_dose, sero = all_sero)
}



#' Convert Immune Profile to Efficacies
#'
#' This function takes an immune profile and converts it into a list of efficacies for different products.
#'
#' @param immune_profile A list containing immune profile information. It should have the following structure:
#' \describe{
#'   \item{vhr}{A list with element \code{b} representing the efficacy for very high risk (VHR) individuals.}
#'   \item{mass}{A list with elements \code{b} and \code{product} representing the efficacy and product type for mass immunization. The \code{product} can be "lav", "mab", or "mat".}
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{mab_vhr}{Efficacy for monoclonal antibodies (MAB) in very high risk individuals.}
#'   \item{mab_mass}{Efficacy for monoclonal antibodies (MAB) in mass immunization.}
#'   \item{lav_mass}{Efficacy for live attenuated vaccines (LAV) in mass immunization.}
#'   \item{mat_mass}{Efficacy for maternal antibodies (MAT) in mass immunization.}
#' }
#'
#' @details
#' The function initializes the efficacy list with zeros. It then checks if the immune profile contains efficacy values for VHR and mass immunization and updates the efficacy list accordingly.
#'
#' @examples
#' immune_profile <- list(
#'   vhr = list(b = 0.8),
#'   mass = list(b = 0.6, product = "lav")
#' )
#' convert_efficacies(immune_profile)
#'
#' @export
convert_efficacies <- function(immune_profile) {
    # add info on waning rates here
    efficacy_list <- list(
        mab_vhr = 0,
        mab_mass = 0,
        lav_mass = 0,
        mat_mass = 0
    )
    if(!is.null(immune_profile$vhr$b)) {
        efficacy_list$mab_vhr <- immune_profile$vhr$b
    }
    if(!is.null(immune_profile$mass$b)) {
        if (immune_profile$mass$product == "lav") {
            efficacy_list$lav_mass <- immune_profile$mass$b
        } else if (immune_profile$mass$product == "mab") {
            efficacy_list$mab_mass <- immune_profile$mass$b
        } else if (immune_profile$mass$product == "mat") {
            efficacy_list$mat_mass <- immune_profile$mass$b_mat
            efficacy_list$lav_mass <- immune_profile$mass$b
        }
    }
    efficacy_list
}


#' @title Get Coverage
#' @description This function calculates the coverage based on the immune profile mass and calibration data.
#' @param immune_profile_mass A data frame containing the immune profile mass information.
#' @param cal A data frame containing the calibration data.
#' @return The coverage value. If the product in the immune profile mass is "mat", it returns the mean of the positive coverage values from the calibration data. Otherwise, it returns 0.
#' @examples
#' immune_profile_mass <- data.frame(product = "mat")
#' cal <- data.frame(cov = c(0.1, 0.2, 0, -0.1))
#' get_coverage(immune_profile_mass, cal)
#' @export
get_coverage <- function(immune_profile_mass, cal) {
     if(immune_profile_mass$product == "mat") {
        cov <- mean(cal$cov[cal$cov > 0])
     } else {
        cov <- 0
     }
     cov
}