library(dplyr)
library(tidyr)

#### Data Loading & Clean-up ####
# let's load in the species bin counts by site x month
stats_df <- read.csv("./data/stats_df.csv")

# and remove Hemipterans; we're not interested in this order for this analysis
stats_df %>%
  filter(Order != "Hemiptera") -> rawData
rm(stats_df)
rawData %>%
  expand(Exact.Site, Month_Year, BIN) -> fullData
# now let's use lubridate to turn the month/year timestamp into actual time data
fullData$Month <- lubridate::parse_date_time(fullData$Month_Year, "b-y")
# this adds a column to the right end of the dataframe with a date for each
# sample. the date defaults to the 1st of the month, but that doesn't matter for
# our analysis

# now let's read in our sampling effort data frame
sampling_days <- read.csv("./data/sampling_days.csv")
# and get the correct time data
sampling_days$Month <- lubridate::parse_date_time(sampling_days$Month_Year, 
                                                  "y-b")
# Now let's joing the sample days info with the expanded dataset
tmp <- left_join(fullData, sampling_days, by=c("Exact.Site" = "Site",
                                               "Month" = "Month"))
# and remove the "Month_Year" columns that are no longer useful
tmp %>%
  select(-c(Month_Year.x, Month_Year.y)) -> tmp
# now we'll add the abundance data to the correct places
rawData$Month <- lubridate::parse_date_time(rawData$Month_Year, "b-y")
tmp2 <- right_join(rawData, tmp, by=c("Exact.Site" = "Exact.Site",
                                      "Month" = "Month",
                                      "BIN" = "BIN"))
tmp2 %>%
  select(Exact.Site, Month, BIN, Abundnace, Days_on_Trap) -> tmp2
tmp2 %>%
  mutate(
    Abundnace = replace_when(Abundnace, is.na(Abundnace) & Days_on_Trap > 0 ~ 0)
  ) -> tmp1
tmp1$Days_on_Trap[is.na(tmp1$Days_on_Trap)] <- 0

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
rawData %>%
  select(Order, BIN) -> orders
order.key <- data.frame(Order = sort(unique(orders$Order)))
order.key$Order_ID <- as.numeric(factor(order.key$Order))
# we also need an index of which species bins belong to which Orders
orders$Order <- as.numeric(as.factor(orders$Order))
orders$BIN <- as.numeric(as.factor(orders$BIN))
orders <- unique(orders[,c('Order','BIN')])
rm(tmp, tmp2)
# 20260417 - I might have to use nested indexing to handle model computation
# has_data <- which(!is.na(y), arr.ind = TRUE)
# y_long <- y[!is.na(y)]

#### Covariate Processing ####
# our analysis includes 2 covariates on capture rates: average wind speed during
# the sampling period (wind), & number of trap-days (J)
# these covariates are also indexed by site/time, so it gets a similar array to 
# the species bin data
rawData %>%
  select(Exact.Site, Month, wind_speed_ms_mean, PM2.5, n_smoke, precipitation_accumulation_mm, 
         max_relative_humidity_mean, max_air_temperature_mean_K) %>%
  unique() -> covariates
covariates$Exact.Site <- as.numeric(as.factor(covariates$Exact.Site))
covariates$Month <- as.numeric(as.factor(covariates$Month))
covariates$wind_speed_ms_mean <- scale(covariates$wind_speed_ms_mean)
covariates$PM2.5 <- scale(covariates$PM2.5)
covariates$n_smoke <- scale(covariates$n_smoke)
covariates$max_relative_humidity_mean <- scale(covariates$max_relative_humidity_mean)
covariates$max_air_temperature_mean_K <- scale(covariates$max_air_temperature_mean_K)
covariates$precipitation_accumulation_mm <- scale(covariates$precipitation_accumulation_mm)

# making covariate arrays
wind <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                 max(covariates$Month)))
for(i in 1:nrow(covariates)){
  wind[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$wind_speed_ms_mean[i]
}
wind[is.na(wind)] <- 0
# hasData_wind <- which(!is.na(wind), arr.ind = TRUE)
# ob_cov_long <- matrix(1,
#                       nrow(hasData_wind),
#                       ncol = 2)
# for(i in 1:nrow(hasData_wind)){
#   ob_cov_long[i,2] <- wind[
#     hasData_wind[i,1],
#     hasData_wind[i,2]
#   ]
# }
PM <- array(data = NA, dim = c(max(covariates$Exact.Site),
                               max(covariates$Month)))
for(i in 1:nrow(covariates)){
  PM[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$PM2.5[i]
}
PM[is.na(PM)] <- 0

smoke <- array(data = NA, dim = c(max(covariates$Exact.Site),
                               max(covariates$Month)))
for(i in 1:nrow(covariates)){
  smoke[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$n_smoke[i]
}
smoke[is.na(smoke)] <- 0
# temperature
temp <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                  max(covariates$Month)))
for(i in 1:nrow(covariates)){
  temp[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$max_air_temperature_mean_K[i]
}
temp[is.na(temp)] <- 0
# humidity
humid <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                 max(covariates$Month)))
for(i in 1:nrow(covariates)){
  humid[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$max_relative_humidity_mean[i]
}
humid[is.na(humid)] <- 0
# precip
precip <- array(data = NA, dim = c(max(covariates$Exact.Site),
                                  max(covariates$Month)))
for(i in 1:nrow(covariates)){
  precip[
    covariates$Exact.Site[i],
    covariates$Month[i]
  ] <- covariates$precipitation_accumulation_mm[i]
}
precip[is.na(precip)] <- 0

#### Run Model ####
data_list <- list(
  nSite = max(site.key$Site_ID),
  nTaxa = max(tmp1$BIN),
  nOrder = max(order.key$Order_ID),
  nMonth = max(tmp1$Month),
  order = orders$Order,
  y = y,
  J = J,
  wind = wind,
  PM = PM,
  smoke = smoke,
  temp = temp,
  humid = humid,
  precip = precip
)

source("./functions/inits.R")

library(runjags)
# I guess my JAGS isn't stored where {runjags} expects it to be, so I have to 
# tell it where to look
runjags.options(jagspath = "C:/Users/rlarson/AppData/Local/Programs/JAGS/JAGS-4.3.2/x64/bin")
my_mod <- runjags::run.jags(
  model = "./model/JAGS_model.R",
  monitor = c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
              "mu.mu.beta2","mu.tau.beta2",
              "mu.beta0","mu.beta1","mu.beta2","tau.beta0","tau.beta1","tau.beta2",
              "beta3","beta4","beta5",
              "mu.mu.alpha0","mu.tau.alpha0","mu.mu.alpha2","mu.tau.alpha2",
              "mu.alpha0","tau.alpha0","alpha1","mu.alpha2","tau.alpha2",
              "alpha3","alpha4",
              "Ntotal", "Nsite"),
  data = data_list,
  n.chains = 3,
  inits = inits,
  burnin = 100,
  sample = 1000,
  adapt = 100,
  modules = "glm",
  thin = 3,
  method = "parallel"#,
  # jags = runjags.getOption("jagspath")
)
system("say Calculations Complete.") # only works on Mac
varSum <- c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
            "mu.mu.beta2","mu.tau.beta2",
            "mu.beta0","mu.beta1","mu.beta2","tau.beta0","tau.beta1","tau.beta2",
            "beta3","beta4","beta5",
            "mu.mu.alpha0","mu.tau.alpha0","mu.mu.alpha2","mu.tau.alpha2",
            "mu.alpha0","tau.alpha0","alpha1","mu.alpha2","tau.alpha2",
            "alpha3","alpha4")
runjags::add.summary(my_mod, vars = varSum)
plot(my_mod, plot.type = "trace", vars = varSum)
