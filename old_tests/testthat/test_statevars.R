

test_that("check state variables outputs", {

    RSVempty <- make_rsv_programme(S = 1)

    # Load the relevant csv files` 
    econ_raw <- read.csv(file = system.file(package = "rsvie", "extdata", "econ", "econ_pars.csv")) 
    risks_raw <- read.csv(file = system.file(package = "rsvie", "extdata", "econ", "outcome_risks.csv"))
    risks_vhr_raw <- read.csv(file = system.file(package = "rsvie", "extdata", "econ", "outcome_risks_vhr.csv"))
    cal_none <- read.csv(file = system.file(package = "rsvie", "extdata", "calendars", "cal_none.csv")) 
    immune_profile_none <- list(
        mass = list(product = "none"),
        vhr = list(
        product = "none"),
        direct = FALSE
    )

    RSVempty <- add_economics(RSVempty, econ_name = "E_W2023", econ_raw, risks_raw, risks_vhr_raw)
    RSV_none <- add_programme(RSVempty, prog_name = "none", cal_none, cal_none, immune_profile_none)

    load(file = system.file(package = "rsvie", "extdata", "tests", "state_none_180423.RData")) # none_test

    expect_equal(run_state(RSV_none), none_test)
})
