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

# these data.frames will help us connect the model output back to the actual
# species, season, or sites they are associated to. These data.frames will also
# give us the numeric index we will need for the JAGS model
taxa_map <- data.frame(
  Taxa = sort(unique(tmp1$BIN))
)
taxa_map$Taxa_ID <- as.numeric(as.factor(taxa_map$Taxa))

month_map <- data.frame(
  Month = unique(tmp1$Month)
)
month_map$Month_ID <- as.numeric(as.factor(month_map$Month))

site_map <- data.frame(
  Site = sort(unique(tmp1$Exact.Site))
)
site_map$Site_ID <- as.numeric(as.factor(site_map$Site))

# now that we have these we can make a temporary site_data and join the maps
tmp <- dplyr::left_join(
  tmp1,
  taxa_map,
  by = join_by(BIN == Taxa)
) %>%
  dplyr::left_join(
    .,
    month_map,
    by = "Month"
  ) %>%
  dplyr::left_join(
    .,
    site_map,
    by = join_by(Exact.Site == Site)
  )
# now we just need to drop the rows where Days_on_Trap == 0 (i.e., no sampling)
tmp <- tmp[tmp$Days_on_Trap > 0,]

# we have multiple months of sampling and therefore need to account for this in
# the model. The issue here, though, is that each site does not have data for
# each season. As such, we need to be able to index the correct seasonal
# parameter in our long-format capture array
month_density <- t(
  table(
    tmp$Month_ID,
    tmp$Site_ID
  )
)
# give row names that are site names
row.names(month_density) <- site_map$Site
# make this binary; 1 = site has data for a given month
month_density[month_density>0] <- 1
# this tells us which site and month have data
month_has_data <- which(
  month_density == 1,
  arr.ind = TRUE
)
# give useful headers to 'month_has_data'
colnames(month_has_data) <- c("site_id", "month_id")
# order by season they site (same ways a taxa capture data)
month_has_data <- month_has_data[
  order(
    month_has_data[,"month_id"],
    month_has_data[,"site_id"]
  ),
]
# month_has_data is our first step; it's effectively a map that lets us index 
# what sites and months have data. We do model the probability of each taxon at
# a site, but the months and sites that have data vary. As such, we are going to 
# make a unique identifier for each taxon, month, and site.
# It's important to note that 'tmp' is sorted by taxon, month, then site. Sites
# that don't have data for a given month are removed at this point in the script
# as well. As such, we need to make a unique identifier for every combination of
# taxon, month, and city.
combo <- tmp[
  ,
  c("Taxa_ID", "Month_ID", "Site_ID")
]
# make the unique identifier
combo$Combo_ID <- 1:nrow(combo)
# and then join this back to the taxon capture array
tmp <- dplyr::left_join(tmp,
                        combo,
                        by = c("Taxa_ID","Month_ID","Site_ID"))

# select covariates
rawData %>%
  select(Exact.Site, Month_Year, PM2.5, n_smoke, precipitation_accumulation_mm, max_relative_humidity_mean,
         max_air_temperature_mean_K, wind_speed_ms_mean) -> covs
covs <- covs[
  -which(duplicated(covs)),
]
covs$Month <- lubridate::parse_date_time(covs$Month_Year, "b-y")
covs_centered <- covs %>%
  dplyr::mutate_if(
    is.numeric,
    function(x) c(scale(x))
  )

# now we need to replicate this across our capture matrix
tmp_covs <- dplyr::inner_join(
  tmp[,c("Exact.Site", "Month")],
  covs_centered,
  by = c("Exact.Site","Month")
)
tmp_covs <- tmp_covs[,-c(3)]
# tmp_covs is ordered like tmp, so we should be good to go with this
# for lambda, we are going to use all the covariates except wind speed. We also 
# need to add a column for the intercept
lambda_covs <- matrix(
  1,
  ncol = ncol(tmp_covs)-3+1, # -3 columns (site, month, wind speed); +1 column (intercept)
  nrow = nrow(tmp_covs)
)
lambda_covs[,-1] <- as.matrix(
  tmp_covs[,-which(
    colnames(tmp_covs) %in% c("Exact.Site","Month","wind_speed_ms_mean")
  )]
)
# for capture rates (rho), we will use temp, humidity, and wind speed
rho_covs <- matrix(
  1,
  ncol = ncol(tmp_covs)-5+1,
  nrow = nrow(tmp_covs)
)
rho_covs[,-1] <- as.matrix(
  tmp_covs[,-which(
    colnames(tmp_covs) %in% c("Exact.Site","Month","PM2.5","n_smoke",
                              "precipitation_accumulation_mm")
  )]
)

data_list <- list(
  # capture data
  y = tmp$Abundnace,
  # J = tmp$Days_on_Trap,
  # indexes to fit the model in long format
  taxa_idx = tmp$Taxa_ID,
  site_idx = tmp$Site_ID,
  combo_taxa_idx = combo$Taxa_ID,
  combo_site_idx = combo$Site_ID,
  combo_idx = combo$Combo_ID,
  # covariates
  lambda_covs = lambda_covs,
  rho_covs = rho_covs,
  # number of species, parameters, etc.
  nTaxa = max(tmp$Taxa_ID),
  nSite = max(tmp$Site_ID),
  ncov_abundance = ncol(lambda_covs),
  nSamples = nrow(tmp),
  nMonth_Params = nrow(combo),
  ncov_capture = ncol(rho_covs)
)

source("./functions/inits_Indexed.R")

library(runjags)
# I guess my JAGS isn't stored where {runjags} expects it to be, so I have to 
# tell it where to look
runjags.options(jagspath = "C:/Users/rlarson/AppData/Local/Programs/JAGS/JAGS-4.3.2/x64/bin")
my_mod <- runjags::run.jags(
  model = "./model/nestedIndex_model.R",
  monitor = c(# abundance
    "b_abundance","tau_abundance","tau_shape","tau_rate","b_taxa","tau_taxa","b_taxa_site",
    # monthly variation
    "c_shape_lambda","c_rate_lambda","c_tau_lambda",
    # capture rates
    "c_cap","tau_cap","tau_shape_rho","tau_rate_rho","c_taxa_cap","tau_taxa_cap","c_taxa_site",
    # monthly variation capture
    "c_shape_rho","c_rate_rho","c_tau_rho"),
  data = data_list,
  n.chains = 3,
  inits = inits_Indexed,
  burnin = 1000,
  sample = 10000,
  adapt = 100,
  modules = "glm",
  thin = 3,
  method = "parallel",
  jags = runjags.getOption("jagspath")
)
#system("say Calculations Complete.") # only works on Mac
varSum <- c("b_abundance","tau_abundance","tau_shape","tau_rate","b_taxa","tau_taxa","b_taxa_site",
            "c_shape_lambda","c_rate_lambda","c_tau_lambda",
            "c_cap","tau_cap","tau_shape_rho","tau_rate_rho","c_taxa_cap","tau_taxa_cap","c_taxa_site",
            "c_shape_rho","c_rate_rho","c_tau_rho")
runjags::add.summary(my_mod, vars = varSum)
plot(my_mod, plot.type = "trace", vars = varSum)
