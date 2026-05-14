library(dplyr)
library(tidyr)

#### Data Loading & Clean-up ####
# let's load in the species bin counts by site x month
rawData <- read.csv("./data/stats_df.csv")
# now let's use lubridate to turn the month/year timestamp into actual time data
rawData$Month <- lubridate::parse_date_time(rawData$Month_Year, "b-y")
# this adds a column to the right end of the dataframe with a date for each
# sample. the date defaults to the 1st of the month, but that doesn't matter for
# our analysis

# we're going to do some filtering to remove all data before Feb '23 and after
# Oct '23, and we'll remove Hemiptera since we're not interested in this Order
rawData %>%
  filter(between(Month, as.Date("2023-02-01"), as.Date("2023-10-01"))) %>%
  filter_out(Order == "Hemiptera") -> filteredData
rm(rawData)
# now we'll expand our dataframe to include all combos of site x month x taxa
filteredData %>%
  expand(Exact.Site, Month, BIN) -> fullData

# now let's read in our sampling effort data frame
sampling_days <- read.csv("./data/sampling_days.csv")
# and get the correct time data
sampling_days$Month <- lubridate::parse_date_time(sampling_days$Month_Year, 
                                                  "y-b")
# and remove the same dates as before
sampling_days %>%
  filter(between(Month, as.Date("2023-02-01"), as.Date("2023-10-01"))) -> sampling_days

# Now let's join the sample days info with the expanded dataset
tmp <- left_join(fullData, sampling_days, by=c("Exact.Site" = "Site",
                                               "Month" = "Month"))
# and remove the "Month_Year" column that is no longer useful
tmp %>%
  select(-c(Month_Year)) -> tmp
# now we'll add the abundance data to the correct places
tmp2 <- right_join(filteredData, tmp, by=c("Exact.Site" = "Exact.Site",
                                      "Month" = "Month",
                                      "BIN" = "BIN"))
tmp2 %>%
  select(Exact.Site, Month, BIN, Abundnace, Days_on_Trap) -> tmp2
tmp2 %>%
  mutate(
    Abundnace = replace_when(Abundnace, is.na(Abundnace) & Days_on_Trap > 0 ~ 0)
  ) -> tmp1
# replace non-sampled months with 0's if necessary
# tmp1$Days_on_Trap[is.na(tmp1$Days_on_Trap)] <- 0

# let's put the data into the correct order. first sort by species bin, then by
# month/year, then by site
tmp1 <- tmp1[order(
  tmp1[,"BIN"],
  tmp1[,"Month"],
  tmp1[,"Exact.Site"]),]

# and before we go further let's make some reference keys to help us remember
# which number corresponds to which site / month / Order
site.key <- data.frame(Site = sort(unique(tmp1$Exact.Site)))
site.key$Site_ID <- as.numeric(factor(site.key$Site))

month.key <- data.frame(Month_Year = sort(unique(tmp1$Month)))
month.key$Month_ID <- as.numeric(factor(month.key$Month_Year))

# then replace values w/ numbers
tmp1$Exact.Site <- as.numeric(as.factor(tmp1$Exact.Site))
tmp1$Month <- as.numeric(as.factor(tmp1$Month))
tmp1$BIN <- as.numeric(as.factor(tmp1$BIN))

# data array creation
y <- array(data = NA, dim = c(max(tmp1$BIN),
                              max(tmp1$Exact.Site),
                              max(tmp1$Month)))
for(i in 1:nrow(tmp1)){
  y[
    tmp1$BIN[i],
    tmp1$Exact.Site[i],
    tmp1$Month[i]
  ] <- tmp1$Abundnace[i]
}

J <- array(data = NA, dim = c(max(tmp1$Exact.Site),
                              max(tmp1$Month)))
for(i in 1:nrow(tmp1)){
  J[
    tmp1$Exact.Site[i],
    tmp1$Month[i]
  ] <- tmp1$Days_on_Trap[i]
}
# getting order information
filteredData %>%
  select(Order, BIN) -> orders
order.key <- data.frame(Order = sort(unique(orders$Order)))
order.key$Order_ID <- as.numeric(factor(order.key$Order))
# we also need an index of which species bins belong to which Orders
orders$Order <- as.numeric(as.factor(orders$Order))
orders$BIN <- as.numeric(as.factor(orders$BIN))
orders <- unique(orders[,c('Order','BIN')])
rm(tmp, tmp2)

#### Covariate Processing ####
# our analysis includes 2 covariates on capture rates: average wind speed during
# the sampling period (wind), & number of trap-days (J)
# these covariates are also indexed by site/time, so it gets a similar array to 
# the species bin data
metData <- read.csv("./data/abiotic_var/gridMET_SDNHMsites_monthlyvalues_13may26.csv")
metData$month <- lubridate::parse_date_time(metData$month,"y-m-d")
metData %>%
  filter(between(month, as.Date("2023-02-01"), as.Date("2023-10-01"))) -> metData
metData$site.name <- as.numeric(as.factor(metData$site.name))
metData$month <- as.numeric(as.factor(metData$month))
metData$wind_speed_ms_mean <- scale(metData$wind_speed_ms_mean)
metData$max_relative_humidity_mean <- scale(metData$max_relative_humidity_mean)
metData$max_air_temperature_mean_K <- scale(metData$max_air_temperature_mean_K)
metData$precipitation_accumulation_mm <- scale(metData$precipitation_accumulation_mm)

filteredData %>%
  select(Exact.Site, Month, PM2.5, n_smoke) %>%
  unique() -> tmp
tmp$Exact.Site <- as.numeric(as.factor(tmp$Exact.Site))
tmp$Month <- as.numeric(as.factor(tmp$Month))
tmp$PM2.5 <- scale(tmp$PM2.5)
tmp$n_smoke <- scale(tmp$n_smoke)

covariates <- right_join(tmp, metData, by=c("Exact.Site" = "site.name",
                                           "Month" = "month"))

# making covariate arrays
# particulate matter (i.e., air pollution)
PM <- array(data = NA, dim = c(max(covariates$Exact.Site),
                               max(covariates$Month)))
for(i in 1:nrow(covariates)){
  PM[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$PM2.5[i]
}
PM[is.na(PM)] <- 0
# smokey days
smoke <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                  max(covariates$Month)))
for(i in 1:nrow(covariates)){
  smoke[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$n_smoke[i]
}
smoke[is.na(smoke)] <- 0
# PC1 (+ = wet/humid siteMonths | - = hot siteMonths)
PC1 <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                  max(covariates$Month)))
for(i in 1:nrow(covariates)){
  PC1[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$PC1[i]
}
PC1[is.na(PC1)] <- 0
# PC2 (+ = calm siteMonths | - = windy siteMonths)
PC2 <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                max(covariates$Month)))
for(i in 1:nrow(covariates)){
  PC2[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$PC2[i]
}
PC2[is.na(PC2)] <- 0
# wind speed
wind <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                 max(covariates$Month)))
for(i in 1:nrow(covariates)){
  wind[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$wind_speed_ms_mean[i]
}
wind[is.na(wind)] <- 0

#### Run Model ####
data_list <- list(
  nSite = max(site.key$Site_ID),
  nTaxa = max(tmp1$BIN),
  nOrder = max(order.key$Order_ID),
  nMonth = max(tmp1$Month),
  order = orders$Order,
  y = y,
  J = J,
  PM = PM,
  smoke = smoke,
  PC1 = PC1,
  PC2 = PC2
)

source("./functions/inits.R")

library(runjags)
# I guess my JAGS isn't stored where {runjags} expects it to be, so I have to 
# tell it where to look
runjags.options(jagspath = "C:/Users/rlarson/AppData/Local/Programs/JAGS/JAGS-4.3.2/x64/bin")
my_mod <- runjags::run.jags(
  model = "./model/JAGS_model.R",
  monitor = c(# hyper-hyperpriors for abundance
              "mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
              "mu.mu.beta2","mu.tau.beta2","mu.mu.beta3","mu.tau.beta3",
              "mu.phi","tau.phi",
              # hyperpriors for order-level responses in abundance
              "mu.beta0","mu.beta1","mu.beta2","mu.beta3",
              "tau.beta0","tau.beta1","tau.beta2","tau.beta3",
              # hyper-hyperpriors for capture/detection
              "mu.mu.alpha0","mu.tau.alpha0","mu.mu.alpha2","mu.tau.alpha2",
              "mu.mu.alpha3","mu.tau.alpha3",
              # hyperpriors for order-level responses in capture rates
              "mu.alpha0","tau.alpha0","alpha1","mu.alpha2","tau.alpha2",
              "mu.alpha3","tau.alpha3",
              # derived parameters
              "Ntotal", "Nsite"),
  data = data_list,
  n.chains = 3,
  inits = inits,
  burnin = 100,
  sample = 1000,
  adapt = 100,
  modules = "glm",
  thin = 3,
  method = "parallel",
  jags = runjags.getOption("jagspath")
)

varSum <- c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
            "mu.mu.beta2","mu.tau.beta2","mu.mu.beta3","mu.tau.beta3",
            "mu.phi","tau.phi",
            "mu.beta0","mu.beta1","mu.beta2","mu.beta3",
            "tau.beta0","tau.beta1","tau.beta2","tau.beta3",
            "mu.mu.alpha0","mu.tau.alpha0","mu.mu.alpha2","mu.tau.alpha2",
            "mu.mu.alpha3","mu.tau.alpha3",
            "mu.alpha0","tau.alpha0","alpha1","mu.alpha2","tau.alpha2",
            "mu.alpha3","tau.alpha3")
runjags::add.summary(my_mod, vars = varSum)
plot(my_mod, plot.type = "trace", vars = varSum)
