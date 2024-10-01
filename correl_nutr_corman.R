#Code to compare lake and stream TN and TP concentrations from Corman et al 2023 data set
# DG, August 2024

# download Corman et al data set from here: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1268.3
# need stream and lake nutrient concentrations
# packages 

# packages 
library(tidyverse)
library(lubridate)

# load data sets
corman.lake.nutrients <- read_csv("AllLakes_rawinlakenuts.csv")
# stream nutrients 
corman.stream.nutrients  <- read_csv("raw_stream_nutrients.csv")

### Monthly medians

corman.lake.nutrients <- corman.lake.nutrients %>%
  mutate(DateTime = parse_date_time(datetime, orders = c("ymd", "ymd HMS")),Year = year(DateTime),
         Month = month(DateTime)) %>%
  group_by(Lake, Year, Month) %>% 
  summarise(
    TP_lake = median(TP, na.rm = T),
    TN_lake = median(TN, na.rm = T))

# stream nutrients
corman.stream.nutrients <- corman.stream.nutrients %>%
  mutate(DateTime = parse_date_time(Datetime, orders = c("ymd", "ymd HMS")),Year = year(DateTime),
         Month = month(DateTime)) %>%
  group_by(Lake, Year, Month) %>%
  summarise(
    TP_stream = median(TP, na.rm = T),
    TN_stream = median(TN, na.rm = T))

### merge
corman.merged.nutrients <- merge(corman.lake.nutrients, corman.stream.nutrients, 
                                 by = c("Lake", "Year", "Month"))

### quick and dirty linear model
# add teeny-tiny offset for log10
tp.model <- lm(log10(TP_stream + 0.001) ~ log10(TP_lake + 0.001), data = corman.merged.nutrients)
tn.model <- lm(log10(TN_stream + 0.001) ~ log10(TN_lake + 0.001), data = corman.merged.nutrients)

# summaries
summary(tp.model)
summary(tn.model)

# plot relationships
(tp.plt <- corman.merged.nutrients %>%
  ggplot(aes(TP_lake, TP_stream)) + 
  scale_x_log10() + scale_y_log10() + 
  geom_point() + 
  annotate(geom = "text", x = 5, y = 100,
           label = "log10(y) ~ 0.86 + 0.44log10(x) \nR^2 = 0.32, p-value < 0.001, \nF = 34.67 on 70 df") + 
  geom_smooth(method = "lm") + 
  labs(x = "log10(TP lake) ug L^-1", y = "log10(TP stream) ug L^-1"))

(tn.plt <- corman.merged.nutrients %>%
    ggplot(aes(TN_lake, TN_stream)) + 
    scale_x_log10() + scale_y_log10() + 
    geom_point() + 
    annotate(geom = "text", x = 500, y = 20000,
             label = "log10(y) ~ -0.87 + 1.37log10(x) \nR^2 = 0.70, p-value < 0.001, \nF = 158.9 on 66 df") +
    geom_smooth(method = "lm") + 
    labs(x = "log10(TN lake) ug L^-1", y = "log10(TN stream) ug L^-1"))

## arrange figs
correl.nutrients <- ggpubr::ggarrange(
  plotlist = list(tp.plt, tn.plt), ncol = 2, nrow = 1, align = "hv", labels = c("a", "b")
)
correl.nutrients
