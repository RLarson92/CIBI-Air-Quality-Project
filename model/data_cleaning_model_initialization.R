library(dplyr)
library(tidyr)

stats_df <- read.csv("~/GitHub/CIBI-Air-Quality-Project/data/stats_df.csv")
stats_df$Month <- lubridate::parse_date_time(stats_df$Month_Year, "b-y")

stats_df %>%
  group_by(Exact.Site, Month, Order) %>%
  summarise(Y = sum(Abundnace)) -> tmp
tmp <- as.data.frame(tmp)

sampling_days <- read_csv("data/sampling_days.csv")
sampling_days$Month <- lubridate::parse_date_time(sampling_days$Month_Year, "b-y")

tmp1 <- full_join(tmp, sampling_days, by = c("Exact.Site" = "Site", "Month" = "Month"))
tmp1 %>%
  select(-Month_Year) -> tmp1
tmp1 <- as.data.frame(tmp1)

tmp1 %>%
  expand(Exact.Site, Month, Order) -> tmp2
data <- right_join(tmp1, tmp2)
data$Days_on_Trap[is.na(data$Days_on_Trap)] <- 0
