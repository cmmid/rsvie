load(file = here::here("data-raw", "rsv_data_uk.RData")) # loads ukdata
load(file = here::here("data-raw", "inter_data_uk.RData")) # loads data_inter_uk

get_population_age_group <- function(pop, ukdata) {
    pop_tot <- ukdata$population - 3e6
    ageStratification <- ukdata$ageGroupBoundary
    dailyBirthRate <- ukdata$numberDailyLiveBirths

    pop_age <- vector()
    for (i in 1:24){
        pop_age[i] <- dailyBirthRate*365*(ageStratification[i+1]-ageStratification[i])
    }
    pop_age[25] <- pop_tot - sum(pop_age)
    pop_age
}

pop_age <- get_population_age_group(population_tot, ukdata)

# Get the number per risk group 
nVHR <- c(2181.87522,  561.25445,  561.25445,  117.23218,  117.23218, 
    117.23218,   15.63978,   15.63978,   15.63978)

pVHR <- c( nVHR / pop_age[1:9], rep(0, 16))
pHR <- rep(0, 25)
pLR <- c(1 - nVHR / pop_age[1:9], rep(1, 16))

uk_data_sum <- list(
    ageGroupBoundary = ukdata$ageGroupBoundary,
    ageGroupNumber = ukdata$ageGroupNumber,
    numberDailyLiveBirths = ukdata$numberDailyLiveBirths,
    population = ukdata$population - 3e6,
    contactMatrixPhy = ukdata$contactMatrixPhy,
    contactMatrixCon = ukdata$contactMatrixCon,
    observationalData = ukdata$observationalData,
    populationAgeGroup = pop_age,
    prop_mat = data_inter_uk$prop_mat,
    pVHR = pVHR,
    pHR = pHR,
    pLR = pLR,
    nmat = data_inter_uk$prop_mat * pop_age,
    nVHR = pVHR * pop_age,
    nHR = pHR * pop_age,
    nLR = pLR * pop_age
)

uk_data_sum <- c(uk_data_sum, data_inter_uk[c(5:16)])
save(uk_data_sum, file = here::here("inst", "extdata", "data-raw", "uk_data_sum.RData"))