
### phytoSTOICH
### GPP and algae stoichiometry framework simulations
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
source("models/static_liebig_no_light.R") # model 1 in Carly's framework
# Droop model
source("models/dynamic_liebig_no_light.R") # model 3 in Carly's framework
# set timesteps
times <- 1:2000 # updated to 2k for publication


################################################################################
# Simulations

# Changed sims to correspond to the "linear" segment of phytoplankton stoichiometry from Klausmeier et al 2004 (also reported in Meunier et al 2014: https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0107737&type=printable)
# molar N:P = approx 20 to 40/50 = N:P mass of 9 to 20~ish
# set loads
# low, mid, high P supply --> trophic state
loads <- expand.grid(Pin = c(50, 100, 500),
                     NP_inflow = seq(5, 30, 1))
loads$Nin <- loads$Pin * loads$NP_inflow


# Static model
static.runs <-  lapply(list(static.algae, static.diatoms, static.greens, static.cyanos), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    y <- c("A1" = 100, "P" = loads[i, "Pin"],"N" = loads[i, "Nin"])
    run <- ode(y, times, parms = params, func = static.stoich)
    return(run[max(times),])
  })
})
# extract from list and convert to data-frame
static.runs.average <- as_data_frame(do.call(rbind, static.runs[[1]]))
static.runs.diatoms <- as_data_frame(do.call(rbind, static.runs[[2]]))
static.runs.greens <- as_data_frame(do.call(rbind, static.runs[[3]]))
static.runs.cyanos <- as_data_frame(do.call(rbind, static.runs[[4]]))

# Dynamic model

dynamic.runs <-  lapply(list(dynamic.algae, dynamic.diatoms, dynamic.greens, dynamic.cyanos), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    y <- c("A1" = 100, "P" = loads[i, "Pin"],"N" = loads[i, "Nin"], "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich)
    return(run[max(times),])
  })
})
# extract from list and convert to data-frame
dynamic.runs.average <- as_data_frame(do.call(rbind, dynamic.runs[[1]]))
dynamic.runs.diatoms <- as_data_frame(do.call(rbind, dynamic.runs[[2]]))
dynamic.runs.greens <- as_data_frame(do.call(rbind, dynamic.runs[[3]]))
dynamic.runs.cyanos <- as_data_frame(do.call(rbind, dynamic.runs[[4]]))

## bind all together
gpp.sims <- bind_rows(static.runs.average, static.runs.diatoms, static.runs.greens, static.runs.cyanos, 
                      dynamic.runs.average, dynamic.runs.diatoms, dynamic.runs.greens, dynamic.runs.cyanos) %>%
  mutate(model = rep(c("static", "dynamic"), each = nrow(loads) * 4),
         species = rep(c("average", "diatoms", "greens", "cyanos","average", "diatoms", "greens", "cyanos"),
                       each = nrow(loads)),
         Pin = rep(loads$Pin, 8), 
         Pin = paste0("Pin = ", Pin),
         Pin = factor(Pin, levels = paste0("Pin = ", c(50, 100, 500))),
         Nin = rep(loads$Nin, 8), 
         NP_inflow = rep(loads$NP_inflow, 8)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")),
         species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")),
         Pin = factor(Pin),
         NPin_molar = (NP_inflow/14.007)/(1/30.974))


################################################################################
# Plotting

# subset data for easy plotting
average <- gpp.sims[gpp.sims$species == "average",]
diatoms <- gpp.sims[gpp.sims$species == "diatoms",]
greens <- gpp.sims[gpp.sims$species == "greens",]
cyanos <- gpp.sims[gpp.sims$species == "cyanos",]

# data set for seston
seston <- gpp.sims  %>%
  # add in molar conversions
  mutate("C:N" = (1/12.01)/(QN1/14.007),"C:P" = (1/12.001)/(QP1/30.974), 
         "N:P" = (QN1/14.007)/(QP1/30.974), NP_inflow = NPin_molar) %>%
  select(`C:N`, `C:P`, `N:P`,Pin, NP_inflow, species, GPP) %>%
  gather("seston", "value", -Pin, -NP_inflow, -species, -GPP) %>%
  mutate(seston = factor(seston, levels = c("C:N", "C:P", "N:P"))) %>% 
  drop_na()

# subset for plotting
seston.average <- seston[seston$species == "average",]
seston.diatoms <- seston[seston$species == "diatoms",]
seston.greens <- seston[seston$species == "greens",]
seston.cyanos <- seston[seston$species == "cyanos",]

# consumption vctr
consump.vctr <- tibble(
  model = rep(c("static", "dynamic"), each = 4), 
  species = rep(c("average", "diatoms", "greens", "cyanos"), 2), 
  "minQN_minQP" = rep(c(
    ((0.09/14.007)/(1/12.001))/((0.0105/30.974)/(1/12.001)),
    ((0.155/14.007)/(1/12.001))/((0.02/30.974)/(1/12.001)),
    ((0.025/14.007)/(1/12.001))/((0.001/30.974)/(1/12.001)),
    ((0.01/14.007)/(1/12.001))/((0.001/30.974)/(1/12.001))), 2),
  half_sat_P = rep(c(0.005, 0.028, 0.026, 0.0165)/30.974 * 1000, 2),
  half_sat_N = rep(c(0.064, 0.036, 0.033, 0.05)/14.007 * 1000,2))
  

# create seston summary w. minQN and minQP
consump.vctr.seston <- seston %>%
  merge(consump.vctr, by = "species") %>%
  group_by(species, seston, Pin) %>%
  summarise(minQN_minQP = max(minQN_minQP), 
            CNP = max(value),
            GPP = max(GPP),
            half_sat_P = min(half_sat_P),
            half_sat_N = min(half_sat_N))

# manual legend
species.legend = c("average", "diatoms","greens", "cyanos")

# plot GPP across supply N:P for each group
(gpp.plt <- ggplot() + 
  # # add in the consumption vctrs 
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin), 
               lwd = 1, lty = "dashed") +
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species, group = Pin), 
               lwd = 0.75, lty = "dashed") +
  # plot GPP
  geom_line(data = gpp.sims, aes(x = NPin_molar, y = GPP,col = species, 
                                 group = interaction(species, Pin)), lwd = 0.75) +
  geom_point(data = gpp.sims, aes(x = NPin_molar, y = GPP, fill = species, pch = species, 
                                  group = Pin), size = 2) + 
  scale_x_log10() + scale_y_log10() + 
  ggh4x::facet_grid2(Pin~model) + 
  scale_shape_manual(values = c(21, 22, 24, 25)) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() +
  labs(x = "Load N:P (molar)", 
       y = expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
       fill = NULL, col = NULL, pch = NULL))


# plot CNP across supply N:P for each species
# Cedric wants this plot broken up by C:N, C:P and N:P

# C:N 
(cn.plt <- ggplot() +
    # add in the consumption vctrs
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "C:N",],
                 aes(y = CNP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "C:N",],
                 aes(y = CNP, x = minQN_minQP, xend = 0, yend = CNP, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_line(data = seston[seston$seston == "C:N",], aes(NP_inflow, value, col = species, group = interaction(species, Pin)),
              lwd = 0.75) +
    geom_point(data = seston[seston$seston == "C:N",], aes(NP_inflow, value, fill = species, pch = species, group = Pin), size = 2) +
    scale_x_log10() + scale_y_log10() +
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(x = "Load N:P (molar)", y = "Phytoplankton C:N (molar)", fill = NULL, col = NULL, pch = NULL) + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
    ))

(cp.plt <- ggplot() +
    # add in the consumption vctrs
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "C:P",],
                 aes(y = CNP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "C:P",],
                 aes(y = CNP, x = minQN_minQP, xend = 0, yend = CNP, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_line(data = seston[seston$seston == "C:P",], aes(NP_inflow, value, col = species, group = interaction(species, Pin)),
              lwd = 0.75) +
    geom_point(data = seston[seston$seston == "C:P",], aes(NP_inflow, value, fill = species, pch = species, group = Pin), size = 2) +
    scale_x_log10() + scale_y_log10() +
    #ggh4x::facet_grid2(.~seston, scales = "free", independent = "y") +
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") +
    labs(x = "Load N:P (molar)", y = "Phytoplankton C:P (molar)", fill = NULL, col = NULL, pch = NULL) + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
    ))

(np.plt <- ggplot() +
    # add in the consumption vctrs
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "N:P",],
                 aes(y = CNP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "N:P",],
                 aes(y = CNP, x = minQN_minQP, xend = 0, yend = CNP, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_line(data = seston[seston$seston == "N:P",], aes(NP_inflow, value, col = species, group = interaction(species, Pin)),
              lwd = 0.75) +
    geom_point(data = seston[seston$seston == "N:P",], aes(NP_inflow, value, fill = species, pch = species, group = Pin), size = 2) +
    scale_x_log10() + scale_y_log10() +
    #ggh4x::facet_grid2(.~seston, scales = "free", independent = "y") +
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") +
    labs(x = "Load N:P (molar)", y = "Phytoplankton N:P (molar)", fill = NULL, col = NULL, pch = NULL))

cnp.plt <- ggarrange(plotlist = list(cn.plt, cp.plt, np.plt), 
                     ncol = 3, labels = c("a", "b", "c"), 
                     common.legend = T, legend = "bottom", 
                     align = "hv")
cnp.plt

################################################################################

# save figures individually
if(!dir.exists("/figures")) {dir.create("/figures")}
save_plot("figures/gpp_framework.png", gpp.plt, base_width = 5, base_height = 8)
save_plot("figures/seston_framework.png", cnp.plt, base_height = 6, base_width = 8)


################################################################################
