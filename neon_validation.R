### phytoSTOICH
### Seston C:N:P validation
### DG, September 2024

# set wd
rm(list = ls())
setwd("C:/Users/DanielGschwentner/Documents/GitHub/WoodStoich24_lake_models")
###############################################################################
# Setup

#load packages
pck <- c("deSolve", "tidyverse", "cowplot", "ggthemes","ggpubr")
lapply(pck, require, character.only = T)
theme_set(theme_pubr() + theme(legend.position = "bottom"))

# Load algae parameters
source("algae_param_vctrs.R")

# Load models
# static model with fixed stoichiometry
source("models/static_liebig_zmix.R") # model 1 in Carly's framework
# Droop model
source("models/dynamic_liebig_zmix.R") 
# set timesteps
times <- 1:2000 # for troubleshooting, initial runs

# load data
NEON_Flat <- read_csv("data4input/NEON_Flathead_forcing.csv")
NEON_Flat <- as.data.frame(NEON_Flat)

################################################################################

### base predictions w. median values

# dynamic model
(start <- Sys.time())
NEON.dynamic <-  lapply(list(dynamic.algae, dynamic.diatoms, dynamic.greens, dynamic.cyanos), function(x) {
  params <- x
  lapply(1:nrow(NEON_Flat), function(i) {
    params["Pin"] = NEON_Flat[i, "Pin_mgm-3"]
    params["Nin"] = NEON_Flat[i, "Nin_mgm-3"]
    params["DOC"] = NEON_Flat[i, "DOC_mgL-1"]
    params["z"] = NEON_Flat[i, "mean depth_m"]
    params["SA"] = NEON_Flat[i, "SA_km2"]
    params["HRT"] = NEON_Flat[i, "HRT_years_mixedlayer"]
    y <- c("A1" = 100, "P" =NEON_Flat[i, "Pin_mgm-3"], "N" = NEON_Flat[i, "Nin_mgm-3"], 
           "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# extract from list and convert to data-frame
NEON.dynamic.average <- as_data_frame(do.call(rbind, NEON.dynamic[[1]]))
NEON.dynamic.diatoms <- as_data_frame(do.call(rbind, NEON.dynamic[[2]]))
NEON.dynamic.greens <- as_data_frame(do.call(rbind, NEON.dynamic[[3]]))
NEON.dynamic.cyanos <- as_data_frame(do.call(rbind, NEON.dynamic[[4]]))

################################################################################

### bind everything together and clean data
NEON.CNP <- bind_rows(NEON.dynamic.average, 
                      NEON.dynamic.diatoms,
                      NEON.dynamic.greens,
                      NEON.dynamic.cyanos)
# add-in information
NEON.CNP$Lake <- rep(NEON_Flat$lake, 4)
NEON.CNP$species <- rep(c("average", "diatoms", "greens", "cyanos"), each = nrow(NEON_Flat))

# save for Helena ternary plot
#write_csv(NEON.CNP, "figures/NEON_Flathead_model_seston_estimates.csv")

# some mutating and cleaning
NEON.CNP <- NEON.CNP %>%
  mutate(CN_mass_model = 1/QN1, 
         CP_mass_model = 1/QP1, 
         NP_mass_model = QN1/QP1) %>%
  rename(P_mod_ugL = P, N_mod_ugL = N) %>%
  mutate(species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")))

# bind with original data set
NEON.valid <- merge(NEON_Flat, NEON.CNP, by.x = "lake", by.y = "Lake")

#################################################################################

# plot predictions following Kelly et al 2018 F5 panels C and D

(obs.plt <- NEON_Flat %>%
   select(Pin = `Pin_mgm-3`, Nin = `Nin_mgm-3`, "C:N" = Seston_CN_mass, "C:P" = Seston_CP_mass, "N:P" = Seston_NP_mass, lake) %>%
   gather("Seston", "Seston.Value", -lake, -Pin, -Nin) %>%
   ggplot(aes(Pin, Seston.Value, col = Nin/Pin)) + 
   geom_point(size = 2) + 
   scale_x_log10() + scale_y_log10() + 
   scale_color_viridis_c() + 
   ggh4x::facet_grid2(.~Seston, scales = "free", independent = "y") + 
   labs(x = expression("Load TP ug L"^-1), 
        y = expression("C:N:P (mass)"),
        col = "Load N:P (mass)"))

(mod.plt <- NEON.valid %>%
    select(Pin = `Pin_mgm-3`, Nin = `Nin_mgm-3`, "C:N" = CN_mass_model, "C:P" = CP_mass_model, "N:P" = NP_mass_model, lake, species) %>%
    gather("Seston", "Seston.Value", -lake, -Pin, -Nin, -species) %>%
    ggplot(aes(Pin, Seston.Value, col = Nin/Pin, pch = species)) + 
    geom_point(size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    scale_color_viridis_c() + 
    ggh4x::facet_grid2(.~Seston, scales = "free", independent = "y") + 
    labs(x = expression("Load TP ug L"^-1), 
         y = expression("C:N:P (mass)"),
         col = "Load N:P (mass)"))

### arrange
subplot.b <- ggarrange(plotlist = list(obs.plt, mod.plt), 
                       nrow = 1, ncol = 2, align = "hv",
                       widths = c(1/2, 1/2),
                       labels = c("b", "c"), 
                       common.legend = T, legend = "bottom")
subplot.b

# 1:1 plot
NEON.valid.wide <- NEON.valid %>% select(lake, Pin = `Pin_mgm-3`, Nin = `Nin_mgm-3`, species,
                          Obs.CN = Seston_CN_mass, Obs.CP = Seston_CP_mass, Obs.NP = Seston_NP_mass,
                          Mod.CN = CN_mass_model, Mod.CP = CP_mass_model, Mod.NP = NP_mass_model) %>%
  gather("Key", "Value", -Pin, -Nin, -lake, -species) %>%
  mutate(source = grepl("Obs", Key), 
         source = ifelse(source == T, "Observed", "Modeled"),
         seston = ifelse(grepl("CN", Key), "C:N",
                         ifelse(grepl("CP", Key), "C:P", "N:P"))) %>%
  pivot_wider(id_cols = c(seston, lake, species), names_from = source, values_from = Value)

(predicted.plt <-  NEON.valid.wide %>%
  ggplot(aes(Observed, Modeled, pch = species)) + 
  geom_point(size = 2) + 
  scale_x_log10() + scale_y_log10() + 
  ggh4x::facet_grid2(.~seston, scales = "free", independent = "all") + 
  labs(x = "Observed C:N:P (mass)", y = "Modeled C:N:P (mass)", shape = NULL))

#################################################################################

### arrange everything into one figure
subplot.a <- predicted.plt
neon.validation.plt <- ggarrange(plotlist = list(subplot.a, subplot.b),
                                   nrow = 2, ncol = 1, align = "hv", 
                                   labels = c("a", ""))
neon.validation.plt
#save_plot("figures/seston_validation.png", neon.validation.plt, base_height = 8, base_width = 14)




                   