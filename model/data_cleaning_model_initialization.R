library(dplyr)
library(tidyr)

stats_df <- read.csv("~/GitHub/CIBI-Air-Quality-Project/data/stats_df.csv")
stats_df$Month <- lubridate::parse_date_time(stats_df$Month_Year, "b-y")

# stats_df %>%
#   group_by(Exact.Site, Month, Order) %>%
#   summarise(Y = sum(Abundnace)) -> tmp
# tmp <- as.data.frame(tmp)

sampling_days <- read.csv("./data/sampling_days.csv")
sampling_days$Month <- lubridate::parse_date_time(sampling_days$Month_Year, "b-y")

tmp <- full_join(stats_df, sampling_days, by = c("Exact.Site" = "Site", "Month" = "Month"))
tmp %>%
  select(-c(Month_Year.x, Month_Year.y)) -> tmp1
tmp1 <- as.data.frame(tmp1)

# tmp1 %>%
#   expand(Exact.Site, Month, BIN) -> tmp2
# data <- right_join(tmp1, tmp2)
# data$Days_on_Trap[is.na(data$Days_on_Trap)] <- 0
# 
# order <- cbind(stats_df$Order, stats_df$BIN)
# colnames(order) <- c("Order", "BIN")
# order <- as.data.frame(order)
# data <- left_join(data, order, by = "BIN")

tmp1 <- tmp1[order(
  tmp1[,"BIN"],
  tmp1[,"Month"],
  tmp1[,"Exact.Site"]),]

# Generate keys
site.key <- data.frame(Site = sort(unique(tmp1$Exact.Site)))
site.key$Site_ID <- as.numeric(factor(site.key$Site))

month.key <- data.frame(Month_Year = sort(unique(tmp1$Month)))
month.key$Month_ID <- as.numeric(factor(month.key$Month_Year))

order.key <- data.frame(Order = sort(unique(tmp1$Order)))
order.key$Order_ID <- as.numeric(factor(order.key$Order))

# then replace values w/ numbers
tmp1$Exact.Site <- as.numeric(as.factor(tmp1$Exact.Site))
tmp1$Month <- as.numeric(as.factor(tmp1$Month))
tmp1$Order <- as.numeric(as.factor(tmp1$Order))
tmp1$BIN <- as.numeric(as.factor(tmp1$BIN))

# array creation
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
J[is.na(J)] <- 0

orders <- unique(tmp1[,c('Order','BIN')])

#### Run Model ####
data_list <- list(
  nSite = max(tmp1$Exact.Site),
  nTaxa = max(tmp1$BIN),
  nOrder = max(tmp1$Order),
  nMonth = max(tmp1$Month),
  order = orders$Order,
  y = y
)

inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      mu.mu.beta0 = rnorm(1),
      mu.tau.beta0 = rgamma(1,1,1),
      mu.mu.beta1 = rnorm(1),
      mu.tau.beta1 = rgamma(1,1,1),
      mu.mu.alpha0 = rnorm(1),
      mu.tau.alpha0 = rgamma(1,1,1),
      mu.beta0 = rnorm(data_list$nOrder),
      mu.beta1 = rnorm(data_list$nOrder),
      mu.alpha0 = rnorm(data_list$nOrder),
      tau.beta0 = rgamma(1,1,1),
      tau.beta1 = rgamma(1,1,1),
      tau.alpha0 = rgamma(1,1,1),
      beta0 = rnorm(data_list$nTaxa),
      beta1 = matrix(1,
                     nrow = data_list$nTaxa,
                     ncol = data_list$nSite),
      alpha0 = rnorm(data_list$nTaxa),
      N = array(2000, dim = c(data_list$nTaxa,
                              data_list$nSite,
                              data_list$nMonth)),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(
    switch(
      chain,
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}

library(runjags)
runjags.options(jagspath = "C:/Users/rlarson/AppData/Local/Programs/JAGS/JAGS-4.3.1/x64/bin")
my_mod <- runjags::run.jags(
  model = "./model/JAGS_model.R",
  monitor = c("mu.mu.beta0","mu.tau.beta0","mu.mu.beta1","mu.tau.beta1","mu.beta0","mu.beta1","tau.beta0","tau.beta1",
              "beta0","beta1",
              "mu.mu.alpha0","mu.tau.alpha0","mu.alpha0","tau.alpha0","alpha0",
              "Nsite"),
  data = data_list,
  n.chains = 3,
  inits = inits,
  burnin = 5000,
  sample = 2000,
  adapt = 100,
  modules = "glm",
  thin = 2,
  method = "parallel",
  jags = runjags.getOption("jagspath")
)
summary(my_mod)
