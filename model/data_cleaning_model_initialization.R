library(dplyr)
library(tidyr)

#### Data Loading & Clean-up ####
# let's load in the species bin counts by site x month
stats_df <- read.csv("./data/stats_df.csv")

# and remove Hemipterans; we're not interested in this order for this analysis
stats_df %>%
  filter(Order != "Hemiptera") -> rawData
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

#### Covariate Processing ####
# our analysis includes 2 covariates on capture rates: average wind speed during
# the sampling period (wind), & number of trap-days (J)
# these covariates are also indexed by site/time, so it gets a similar array to 
# the species bin data
wind <- array(data = NA, dim = c(max(tmp1$Exact.Site),
                                 max(tmp1$Month)))
for(i in 1:nrow(tmp1)){
  wind[
    tmp1$Exact.Site[i],
    tmp1$Month[i]
  ] <- tmp1$wind_speed_ms_mean[i]
}
wind[is.na(wind)]<-mean(wind, na.rm = TRUE)



#### Run Model ####
data_list <- list(
  nSite = max(site.key$Site_ID),
  nTaxa = max(tmp1$BIN),
  nOrder = max(order.key$Order_ID),
  nMonth = max(tmp1$Month),
  order = orders$Order,
  y = y
)

source("./functions/inits.R")

library(runjags)
runjags.options(jagspath = "C:/Users/rlarson/AppData/Local/Programs/JAGS/JAGS-4.3.1/x64/bin")
my_mod <- runjags::run.jags(
  model = "./model/JAGS_model.R",
  monitor = c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
              "mu.beta0","mu.beta1","tau.beta0","tau.beta1",
              "mu.phi","tau.phi",
              "mu.mu.alpha0","mu.tau.alpha0","mu.alpha0","tau.alpha0",
              "Nsite"),
  data = data_list,
  n.chains = 3,
  inits = inits,
  burnin = 10000,
  sample = 2000,
  adapt = 2000,
  modules = "glm",
  thin = 2,
  method = "parallel",
  jags = runjags.getOption("jagspath")
)
varSum <- c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1",
            "mu.beta0","mu.beta1","tau.beta0","tau.beta1",
            "mu.phi","tau.phi",
            "mu.mu.alpha0","mu.tau.alpha0",
            "mu.alpha0","tau.alpha0")
runjags::add.summary(my_mod, vars = varSum)
plot(my_mod, plot.type = "trace", vars = varSum)
