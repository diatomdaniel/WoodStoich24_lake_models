
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
  
  # 
  # "minQN_minQP" = rep(c((0.09/14.007)/(0.0105/30.974), 
  #                       (0.155/14.007)/(0.02/30.974), 
  #                       (0.025/14.007)/(0.001/30.974), 
  #                       (0.01/14.007)/(0.001/30.974)), 2))

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

# plot GPP across supply N:P for each species
(gpp.plt <- ggplot() + 
  # # add in the consumption vctrs ala Tilmann
  # note that consumption vctrs won't match for the static model, only the dynamic as cell quotas change..
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin), 
               lwd = 1, lty = "dashed") +
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species, group = Pin), 
               lwd = 0.75, lty = "dashed") +
  # geom_line(data = average, aes(NP_inflow, GPP, col = Pin, group = Pin), lwd = 1) + 
  # geom_point(data = average, aes(NP_inflow, GPP, fill = Pin, group = Pin), pch = 21, size = 2) + 
  # geom_line(data = diatoms, aes(NP_inflow, GPP, col = Pin, group = Pin), lwd = 1) + 
  # geom_point(data = diatoms, aes(NP_inflow, GPP, fill = Pin, group = Pin), pch = 22, size = 2) + 
  # geom_line(data = greens, aes(NP_inflow, GPP, col = Pin, group = Pin), lwd = 1) + 
  # geom_point(data = greens, aes(NP_inflow, GPP, fill = Pin, group = Pin), pch = 23, size = 2) + 
  # geom_line(data = cyanos, aes(NP_inflow, GPP, col = Pin, group = Pin), lwd = 1) + 
  # geom_point(data = cyanos, aes(NP_inflow, GPP, fill = Pin, group = Pin), pch = 24, size = 2) + 
  #geom_vline(data = consump.vctr, aes(xintercept = minQN_minQP), lwd = 0.75, lty = "dashed") + 
  geom_line(data = gpp.sims, aes(x = NPin_molar, y = GPP,
                                 col = species, group = interaction(species, Pin)),
            lwd = 0.75) +
  geom_point(data = gpp.sims, aes(x = NPin_molar, y = GPP,
                                 fill = species, pch = species, group = Pin),
             size = 2) + 
  scale_x_log10() + scale_y_log10() + 
  #ggh4x::facet_grid2(.~model) + 
  ggh4x::facet_grid2(Pin~model) + 
  scale_shape_manual(values = c(21, 22, 24, 25)) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + 
  #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
  labs(x = "Load N:P (molar)", 
       y = expression("GPP (mg C L"^-1 ~ " day"^-1~")")))


# plot CNP across supply N:P for each species
# Cedric wants this plot broken up by C:N, C:P and N:P
# (cnp.plt <- ggplot() + 
#     # add in the consumption vctrs ala Tilmann
#     geom_segment(data = consump.vctr.seston, 
#                  aes(y = CNP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin), 
#                  lwd = 1, lty = "dashed") + 
#     geom_segment(data = consump.vctr.seston, 
#                  aes(y = CNP, x = minQN_minQP, xend = 0, yend = CNP, col = species, group = Pin), 
#                  lwd = 1, lty = "dashed") + 
#     # data
#     # geom_line(data = seston.average, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.average, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 21, size = 2) + 
#     # geom_line(data = seston.diatoms, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.diatoms, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 22, size = 2) + 
#     # geom_line(data = seston.greens, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.greens, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 23, size = 2) + 
#     # geom_line(data = seston.cyanos, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.cyanos, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 24, size = 2) +
#     geom_line(data = seston, aes(NP_inflow, value, col = species, group = interaction(species, Pin)), 
#               lwd = 0.75) + 
#     geom_point(data = seston, aes(NP_inflow, value, fill = species, pch = species, group = Pin), size = 2) + 
#     scale_x_log10() + scale_y_log10() + 
#     #ggh4x::facet_grid2(.~seston, scales = "free", independent = "y") + 
#     ggh4x::facet_grid2(Pin~seston, scales = "free", independent = "y") + 
#     scale_shape_manual(values = c(21, 22, 23, 24)) + 
#     scale_color_viridis_d() + 
#     scale_fill_viridis_d() + 
#     #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
#     labs(x = "Load N:P (mass)", y = "C:N:P (mass)", col = "species"))


(cn.plt <- ggplot() +
    # add in the consumption vctrs ala Tilmann
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
    #ggh4x::facet_grid2(.~seston, scales = "free", independent = "y") +
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") +
    labs(x = "Load N:P (molar)", y = "Phytoplankton C:N (molar)", col = "species") + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
    ))

(cp.plt <- ggplot() +
    # add in the consumption vctrs ala Tilmann
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
    labs(x = "Load N:P (molar)", y = "Phytoplankton C:P (molar)", col = "species") + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
    ))

(np.plt <- ggplot() +
    # add in the consumption vctrs ala Tilmann
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
    labs(x = "Load N:P (molar)", y = "Phytoplankton N:P (molar)", col = "species"))

cnp.plt <- ggarrange(plotlist = list(cn.plt, cp.plt, np.plt), 
                     ncol = 3, labels = c("a", "b", "c"), 
                     common.legend = T, legend = "bottom", 
                     align = "hv")
cnp.plt

# plot GPP across CNP for each species
# quick edit to consumption vctr.
consump.vctr.seston <- consump.vctr.seston %>% mutate(minQN_minQP = ifelse(seston == "N:P", minQN_minQP, NA))
# (gpp.cnp.plt <- ggplot() + 
#     # # add in the consumption vctrs ala Tilmann
#     geom_segment(data = consump.vctr.seston,
#                  aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin), 
#                  lwd = 0.75, lty = "dashed") +
#     geom_segment(data = consump.vctr.seston,
#                  aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species, group = Pin), 
#                  lwd = 0.75, lty = "dashed") +
#     # data
#     # geom_line(data = seston.average, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.average, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 21, size = 2) + 
#     # geom_line(data = seston.diatoms, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.diatoms, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 22, size = 2) + 
#     # geom_line(data = seston.greens, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.greens, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 23, size = 2) + 
#     # geom_line(data = seston.cyanos, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
#     # geom_point(data = seston.cyanos, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 24, size = 2) + 
#     geom_line(data = seston, aes(value, GPP, col = species, alpha = NP_inflow, group = species),
#               lwd = 0.75) + 
#     geom_point(data = seston, aes(value, GPP, fill = species, pch = species, alpha = NP_inflow, group = species), 
#                size = 2) + 
#     scale_x_log10() + scale_y_log10() + 
#     #ggh4x::facet_grid2(.~seston) + 
#     ggh4x::facet_grid2(Pin~seston, scales = "free", independent = "y") + 
#     scale_shape_manual(values = c(21, 22, 24, 25)) + 
#     scale_color_viridis_d() + 
#     scale_fill_viridis_d() + 
#     #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
#     scale_alpha(range=c(0.5,1), na.value = 0) + 
#     labs(x = "C:N:P (mass)", y =  expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
#          col = "species", alpha = "Load N:P (mass)"))

# Cedric wants figures plotted individually

(gpp.cn <- ggplot() +
    geom_line(data = seston[seston$seston == "C:N",], 
              aes(value, GPP, col = species, alpha = NP_inflow, group = species),
              lwd = 0.75) + 
    geom_point(data = seston[seston$seston == "C:N",], 
               aes(value, GPP, fill = species, pch = species, alpha = NP_inflow, group = species), 
               size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~seston) + 
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") + 
    scale_shape_manual(values = c(21, 22, 24, 25)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    scale_alpha(range=c(0.5,1), na.value = 0) + 
    labs(x = "Phytoplankton C:N (molar)", y =  expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
         col = "species", alpha = "Load N:P (molar)") + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
))

(gpp.cp <- ggplot() +
    geom_line(data = seston[seston$seston == "C:P",], 
              aes(value, GPP, col = species, alpha = NP_inflow, group = species),
              lwd = 0.75) + 
    geom_point(data = seston[seston$seston == "C:P",], 
               aes(value, GPP, fill = species, pch = species, alpha = NP_inflow, group = species), 
               size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~seston) + 
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") + 
    scale_shape_manual(values = c(21, 22, 24, 25)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    scale_alpha(range=c(0.5,1), na.value = 0) + 
    labs(x = "Phytoplankton C:P (molar)", y =  expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
         col = "species", alpha = "Load N:P (molar)") + 
    theme(strip.background = element_blank(),  # Remove the panel border
          strip.text = element_blank()
    ))

(gpp.np <- ggplot() +
        geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "N:P",],
                     aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                     lwd = 0.75, lty = "dashed") +
        geom_segment(data = consump.vctr.seston[consump.vctr.seston$seston == "N:P",],
                     aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species, group = Pin),
                     lwd = 0.75, lty = "dashed") +
    geom_line(data = seston[seston$seston == "N:P",], 
              aes(value, GPP, col = species, alpha = NP_inflow, group = species),
              lwd = 0.75) + 
    geom_point(data = seston[seston$seston == "N:P",], 
               aes(value, GPP, fill = species, pch = species, alpha = NP_inflow, group = species), 
               size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~seston) + 
    ggh4x::facet_grid2(Pin~., scales = "free", independent = "y") + 
    scale_shape_manual(values = c(21, 22, 24, 25)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    scale_alpha(range=c(0.5,1), na.value = 0) + 
    labs(x = "Phytoplankton N:P (molar)", y =  expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
         col = "species", alpha = "Load N:P (molar)"))


cnp.gpp.plt <- ggarrange(plotlist = list(gpp.cn, gpp.cp, gpp.np), 
                     ncol = 3, labels = c("d", "e", "f"), 
                     common.legend = T, legend = "bottom", 
                     align = "hv")
cnp.gpp.plt

################################################################################

# arrange figures
# seston plots go together
seston.figs <- ggarrange(plotlist = list(cnp.plt, cnp.gpp.plt), align = "hv", nrow = 2, 
                         common.legend = T, legend = "bottom")
seston.figs

# save figures individually
#save_plot("figures/gpp_framework.png", gpp.plt, base_width = 5, base_height = 8)
#save_plot("figures/seston_framework.png", seston.figs, base_height = 10, base_width = 7)


################################################################################

### Plot dissolved nutrients

# plot GPP across supply N:P for each species
(p.plt <- ggplot() + 
   # # add in the consumption vctrs ala Tilmann
   # note that consumption vctrs won't match for the static model, only the dynamic as cell quotas change..
   geom_segment(data = consump.vctr.seston,
                aes(y = half_sat_P, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                lwd = 1, lty = "dashed") +
   geom_segment(data = consump.vctr.seston,
                aes(y = half_sat_P, x = minQN_minQP, xend = 0, yend = half_sat_P, col = species, group = Pin),
                lwd = 0.75, lty = "dashed") +
   geom_line(data = gpp.sims, aes(x = NPin_molar, y = P,
                                  col = species, group = interaction(species, Pin)),
             lwd = 0.75) +
   geom_point(data = gpp.sims, aes(x = NPin_molar, y = P,
                                   fill = species, pch = species, group = Pin),
              size = 2) + 
   scale_x_log10() + scale_y_log10() + 
   #ggh4x::facet_grid2(.~model) + 
   ggh4x::facet_grid2(Pin~model) + 
   scale_shape_manual(values = c(21, 22, 24, 25)) + 
   scale_color_viridis_d() + 
   scale_fill_viridis_d() + 
   #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
   labs(x = "Load N:P (molar)", 
        y = expression("Residual P ug L"^-1)))

(n.plt <- ggplot() + 
    # # add in the consumption vctrs ala Tilmann
    # note that consumption vctrs won't match for the static model, only the dynamic as cell quotas change..
    geom_segment(data = consump.vctr.seston,
                 aes(y = half_sat_N, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species, group = Pin),
                 lwd = 1, lty = "dashed") +
    geom_segment(data = consump.vctr.seston,
                 aes(y = half_sat_N, x = minQN_minQP, xend = 0, yend = half_sat_N, col = species, group = Pin),
                 lwd = 0.75, lty = "dashed") +
    geom_line(data = gpp.sims, aes(x = NPin_molar, y = N,
                                   col = species, group = interaction(species, Pin)),
              lwd = 0.75) +
    geom_point(data = gpp.sims, aes(x = NPin_molar, y = N,
                                    fill = species, pch = species, group = Pin),
               size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~model) + 
    ggh4x::facet_grid2(Pin~model) + 
    scale_shape_manual(values = c(21, 22, 24, 25)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    labs(x = "Load N:P (molar)", 
         y = expression("Residual N ug L"^-1)))
