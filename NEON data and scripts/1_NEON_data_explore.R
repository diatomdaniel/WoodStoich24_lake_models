NEON_long <- read.csv("D:\\Research\\Projects\\Woodstoich_5\\NEON_long.csv", header = TRUE, sep = ",")

#FIRST RUN NEON_data_load#

#this script will make graphs for TN, TP, TDN, TDP
#and also 

# do all the prep work --------------------------------------------------
#load necessary packages
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("stringr")

library(tidyr)
library(dplyr)
library("stringr")
library("ggplot2")
library("ggpubr")

#setting a theme
mytheme <- theme_bw()+
  theme(panel.grid.major.x = , #get rid of background etc
        #panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank()
  ) +
  
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  
  theme(legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 23),
        strip.text = element_text(size = 20)) 

# get data for model validation -----------------------------------------------------------
#grab data needed for TN, TP and dissolved N and P
Water_lab <- water$swc_externalLabDataByAnalyte
Water_lab$date_siteID <- paste(Water_lab$siteID, Water_lab$collectDate, sep="_")
Water_lab$date_namedLoc <- paste(Water_lab$namedLocation, Water_lab$collectDate, sep="_")

# delete littoral ones
Water_lab <- Water_lab[grepl("buoy", Water_lab$namedLocation), ]

#make sure units are included in sample name
Water_lab$analyte_unit <- paste(Water_lab$analyte, Water_lab$analyteUnits , sep="_")
#extract necessary data from full frame
Water_lab_2 <- Water_lab  %>% select(date_siteID, #to identify the unique date+location
                                     siteID, collectDate,date_namedLoc,
                                     analyte, analyteUnits,
                                     analyte_unit,analyteConcentration,
                                     remarks)

## pivot data -----------------------------------------------------------

#pivot the table so it's usable
Water_lab_final <- 
  pivot_wider(Water_lab_2,
              id_cols = c(siteID, collectDate),
              names_from = analyte_unit,
              values_from = analyteConcentration,
              values_fn = ~ mean(.x, na.rm = TRUE) #if there are duplicates, take the mean value 
  )


## Calculate the mean analyte concentration + water depth per site -----------------------------------------------------------
Water_lab_Summary <- Water_lab_final %>%
  group_by(siteID) %>%
  summarise(
    across(
      c(TP_milligramsPerLiter, TN_milligramsPerLiter, TDP_milligramsPerLiter,TDN_milligramsPerLiter,
        DOC_milligramsPerLiter), 
      list(mean = ~ mean(.x, na.rm = TRUE)))
  )

Water_depth_Summary <- Water_field_final %>%
  group_by(siteID) %>%
  summarise(
    across(maxDepth, 
           list(mean = ~ mean(.x, na.rm = TRUE)))
  )
)

# plotting for exploration -----------------------------------------------------------

## Nutrients -----------------------------------------------------------

###TN and TP -----------------------------------------------------------

Lakes_TP <- ggplot(Water_lab_final, aes(x=siteID, y=TP_milligramsPerLiter, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("TP (mg/L)")    +
  xlab("Lakes")+
  theme(legend.position = "none")

Lakes_TN <- ggplot(Water_lab_final, aes(x=siteID, y=TN_milligramsPerLiter, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("TN (mg/L)")    +
  xlab("Lakes")+
  theme(legend.position = "none")

### dissolved N and P -----------------------------------------------------------
Lakes_TDN <- ggplot(Water_lab_final, aes(x=siteID, y=TDN_milligramsPerLiter, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("Total Dissolved N (mg/L)")    +
  xlab("Lakes")+
  theme(legend.position = "none")

Lakes_TDP <- ggplot(Water_lab_final, aes(x=siteID, y=TDP_milligramsPerLiter, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("Total Dissolved P (mg/L)")    +
  xlab("Lakes")+
  theme(legend.position = "none")

ggarrange(Lakes_TN,Lakes_TP,
          Lakes_TDN,Lakes_TDP)

## field measurements -----------------------------------------------------------
#oxygen
Lakes_oxygen <- ggplot(Water_field_2, aes(x=siteID, y=dissolvedOxygen, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("Dissolved Oxygen")    +
  xlab("Lakes")+
  theme(legend.position = "none")

#saturation
Lakes_dissolvedOxygenSaturation <- ggplot(Water_field_2, aes(x=siteID, y=dissolvedOxygenSaturation, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  xlab("Lakes")    +
  ylab("Dissolved oxygen (%)")+
  theme(legend.position = "none")

#water temp
Lakes_waterTemp <- ggplot(Water_field_2, aes(x=siteID, y=waterTemp, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("Water temperature")    +
  xlab("Lakes")+
  theme(legend.position = "none")

#water depth
Lakes_MaxWaterDepth <- ggplot(Water_field_2, aes(x=siteID, y=maxDepth, fill=siteID)) +
  geom_boxplot(outlier.shape = NA) +
  #  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  mytheme  +
  scale_fill_brewer(palette='YlOrRd')+
  ylab("Max water depth (m)")    +
  xlab("Lakes")+
  theme(legend.position = "none")

ggarrange(Lakes_dissolvedOxygenSaturation, Lakes_oxygen, Lakes_waterTemp,Lakes_MaxWaterDepth)

# put multivariate data analysis below --------------------------------------------------
