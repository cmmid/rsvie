

#' @title Convert daily uptake from calendar to model
#'
#' @description Convert mean and ci into samples from a fitted distribution.
#' @param object am S4 class of the defined RSV programme
#' @param cal_df a dataframe of a calendar being evaluated
#'
#' @import lubridate
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


#' Function to convert a custom weekly uptake vector into a daily uptake calendar
#'
#' @param up_week_raw A vector of weekly uptake proportions
#' @param start_time Start week when administration occurs
#' @return A vector of the daily proportion of a target group which are vaccinated.
get_daily_uptake <- function(up_week_raw, start_time)
{
    up_week <- vector("numeric", length = 52)
    up_day <- vector("numeric", length = 365)
    for (k in (start_time):(start_time + 20)) { 
        pos_week <- (k - 1) %% 52 + 1
        pos_ist <- (k - start_time) + 1
        up_week[pos_week] <- up_week_raw[pos_ist]
    }
            
    for (i in 1:365) {
        pos_day <- (i - 1 + start_time) %% 365 + 1;
        up_day[pos_day] <- (up_week[((i / 7 + start_time) %% 52) + 1] - up_week[((i / 7 - 1 + start_time) %% 52) + 1]) / 7.0;
        if (up_day[pos_day] < 0)
            up_day[pos_day] <- 0
    }
    up_day
}

#' @title function to plot the dosing calendar of the RSVprogramme]
#' 
#' @param object an RSVProgramme object
#' @import dplyr
#' @return none
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

get_coverage <- function(immune_profile_mass, cal) {
     if(immune_profile_mass$product == "mat") {
        cov <- mean(cal$cov[cal$cov > 0])
     } else {
        cov <- 0
     }
     cov
}
