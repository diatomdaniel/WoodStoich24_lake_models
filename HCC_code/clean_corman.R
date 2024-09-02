# Cody to load and clean Corman et al data set used to force model for grid search
# DG, August 2024


# packages 
library(tidyverse)

# load data sets
# GPP
corman.gpp <- read_csv("Corman2023 metabolism result.csv")
# stream nutrients 
corman.stream.nutrients  <- read_csv("Corman2023 stream nutrients.csv")
# volume and DOC created from S1 T1 and T1 in paper
# no DOC for Lillinonah! used median DOC of the data set = 11.8
corman.vol <- tibble(
  Lake = c("Acton", "Crampton", "EastLong", "Feeagh", "Harp", "Langtjern", "Lillinonah", "Lillsjoliden",
           "Mangstrettjarn", "Mendota", "Morris", "Nastjarn", "Ovre", "Struptjarn", "Trout", "Vortsjarv"),
  z = c(3.9, 5.03, 4.04, 14.50, 13.30, 2, 13, 3.8, 5.3, 12.80, 2.42, 4.20, 4, 3.80, 14.60, 2.80),
  SA = c(2.5, 0.26, 0.03, 3.92, 0.71, 0.23, 6.26, 0.01, 0.02, 39.60, 0.06, 0.01, 0.05, 0.03, 15.70, 270), 
  DOC_mgL = c(4.3, 4.6, 10.6, 5.6, 5.2, 11.8, 11.8, 15.5, 11.9, 4.7, 17.3, 8.2, 23.2, 19.9, 2.8, 11.9)) %>%
  mutate(volume_m3 = SA * 1e6 * z)

### aggregate data and conversions
# GPP
corman.gpp <- corman.gpp %>%
  mutate(DateTime = ymd(DateTime),
         Year = year(DateTime),
         Month = month(DateTime),
         GPP = GPP * 12/32) %>%
  group_by(Lake, Year, Month) %>%
  summarise(
    GPP_sd = sd(GPP, na.rm = T),
    GPP = median(GPP, na.rm = T))
table(is.na(corman.gpp$GPP)) # 12 missing data points, imputed median GPP for a given lake
corman.gpp.filled <- corman.gpp %>%
  group_by(Lake) %>%
  mutate(GPP = ifelse(is.na(GPP), median(GPP, na.rm = TRUE), GPP)) %>%
  ungroup()

# sd GPP > median GPP
table(corman.gpp.filled$GPP_sd > corman.gpp.filled$GPP)
# sd of monthly GPP is generally greater than median values!

# stream nutrients
corman.stream.nutrients <- corman.stream.nutrients %>%
  mutate(DateTime = parse_date_time(Datetime, orders = c("ymd", "ymd HMS")),Year = year(DateTime),
         Month = month(DateTime)) %>%
  group_by(Lake, Year, Month) %>%
  summarise(
    TP_in = median(TP, na.rm = T),
    TN_in = median(TN, na.rm = T))
# missing 4 observations; imputed the median value for the lake/stream
# row 22 Feeagh TN and TP in Nov 2013
# row 97 Trout TN in June 2013
# row 98 Trout TP in June August
# row 100 Trout TP in July 2013
corman.stream.nutrients.filled <- corman.stream.nutrients %>%
  group_by(Lake) %>%
  mutate(TN_in = ifelse(is.na(TN_in), median(TN_in, na.rm = TRUE), TN_in),
         TP_in = ifelse(is.na(TP_in), median(TP_in, na.rm = TRUE), TP_in)) %>%
  ungroup()

## combine data sets
corman2 <- merge(corman.gpp.filled,
                 corman.stream.nutrients.filled, by = c("Lake", "Year", "Month"))
corman2 <- merge(corman2, corman.vol, by = "Lake")

# add calculated zmix
corman2$zmix <- 10^(-0.515 + log10(corman2$DOC_mgL) + 0.115 * log10(2 * sqrt(corman2$SA/pi + 0.991)))
corman2$zmix <- ifelse(corman2$zmix > corman2$z, corman2$z, corman2$zmix)
# if zmix > z, set zmix to z; this is also included in models

# remove everything but corman2 data set
rm(list=setdiff(ls(), "corman2"))