### phytoSTOICH
### minQP and minQN sensitivity; see Figure 4 in olson et al 2022
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

### Simulations for varying minQP
# set loads
# low, mid, high P supply --> trophic state
loads <- expand.grid(Pin = c(50, 100, 500),
                     NP_inflow = seq(5, 30, 1),
                     QP1 = c(0.001, 0.01, 0.05, 0.1))
loads$Nin <- loads$Pin * loads$NP_inflow


# Static model
static.runs <-  lapply(list(static.algae, static.diatoms, static.greens, static.cyanos), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    params["QP1"] = loads[i, "QP1"]
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
    params["minQP1"] = loads[i, "QP1"] # naming is a bit different in load data frame
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
qp.sims <- bind_rows(static.runs.average, static.runs.diatoms, static.runs.greens, static.runs.cyanos, 
                      dynamic.runs.average, dynamic.runs.diatoms, dynamic.runs.greens, dynamic.runs.cyanos) %>%
  mutate(model = rep(c("static", "dynamic"), each = nrow(loads) * 4),
         species = rep(c("average", "diatoms", "greens", "cyanos","average", "diatoms", "greens", "cyanos"),
                       each = nrow(loads)),
         Pin = rep(loads$Pin, 8), 
         Pin = paste0("Pin = ", Pin),
         Pin = factor(Pin, levels = paste0("Pin = ", c(50, 100, 500))),
         Nin = rep(loads$Nin, 8), 
         NP_inflow = rep(loads$NP_inflow, 8),
         minQP = rep(loads$QP1, 8)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")),
         species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")),
         Pin = factor(Pin),
         NPin_molar = (NP_inflow/14.007)/(1/30.974))

################################################################################

### Simulations for varying minQN
# set loads
# low, mid, high P supply --> trophic state
loads <- expand.grid(Pin = c(50, 100, 500),
                     NP_inflow = seq(5, 30, 1),
                     QN1 = c(0.01, 0.05, 0.1, 0.2))
loads$Nin <- loads$Pin * loads$NP_inflow


# Static model
static.runs <-  lapply(list(static.algae, static.diatoms, static.greens, static.cyanos), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    params["QN1"] = loads[i, "QN1"]
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
    params["minQN1"] = loads[i, "QN1"] # naming is a bit different in load data frame
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
qn.sims <- bind_rows(static.runs.average, static.runs.diatoms, static.runs.greens, static.runs.cyanos, 
                     dynamic.runs.average, dynamic.runs.diatoms, dynamic.runs.greens, dynamic.runs.cyanos) %>%
  mutate(model = rep(c("static", "dynamic"), each = nrow(loads) * 4),
         species = rep(c("average", "diatoms", "greens", "cyanos","average", "diatoms", "greens", "cyanos"),
                       each = nrow(loads)),
         Pin = rep(loads$Pin, 8), 
         Pin = paste0("Pin = ", Pin),
         Pin = factor(Pin, levels = paste0("Pin = ", c(50, 100, 500))),
         Nin = rep(loads$Nin, 8), 
         NP_inflow = rep(loads$NP_inflow, 8),
         minQN = rep(loads$QN1, 8)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")),
         species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")),
         Pin = factor(Pin),
         NPin_molar = (NP_inflow/14.007)/(1/30.974))

################################################################################


################################################################################

### Simulations for co-varying miNQN minQP
# set loads
# low, mid, high P supply --> trophic state
loads <- expand.grid(Pin = 100,
                     NP_inflow = c(5,10, 30),
                     QN1 = seq(0.01, 0.5, 0.05),
                     QP1 = seq(0.001, 0.1, 0.005))
loads$Nin <- loads$Pin * loads$NP_inflow


# Static model
static.runs <-  lapply(list(static.algae), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    params["QP1"] = loads[i, "QP1"]
    params["QN1"] = loads[i, "QN1"]
    y <- c("A1" = 100, "P" = loads[i, "Pin"],"N" = loads[i, "Nin"])
    run <- ode(y, times, parms = params, func = static.stoich)
    return(run[max(times),])
  })
})
# extract from list and convert to data-frame
static.runs.average <- as_data_frame(do.call(rbind, static.runs[[1]]))

# Dynamic model

dynamic.runs <-  lapply(list(dynamic.algae), function(x) {
  params <- x
  lapply(1:nrow(loads), function(i) {
    params["Pin"] = loads[i, "Pin"]
    params["Nin"] = loads[i, "Nin"]
    params["minQP1"] = loads[i, "QP1"] # naming is a bit different in load data frame
    params["minQN1"] = loads[i, "QN1"] # naming is a bit different in load data frame
    y <- c("A1" = 100, "P" = loads[i, "Pin"],"N" = loads[i, "Nin"], "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich)
    return(run[max(times),])
  })
})
# extract from list and convert to data-frame
dynamic.runs.average <- as_data_frame(do.call(rbind, dynamic.runs[[1]]))

## bind all together
qp.qn.sims <- bind_rows(static.runs.average,dynamic.runs.average) %>%
  mutate(model = rep(c("static", "dynamic"), each = nrow(loads)), 
         Pin = rep(loads$Pin, 2), 
         Pin = paste0("Pin = ", Pin),
         Nin = rep(loads$Nin, 2), 
         NP_inflow = rep(loads$NP_inflow, 2),
         minQP = rep(loads$QP1, 2),
         minQN = rep(loads$QN1, 2)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")),
         NPin_molar = (NP_inflow/14.007)/(1/30.974),
         NP_inflow2 = paste0("N:P supply = ", NP_inflow))

################################################################################

# plot QP sims
(qp.qn.plot <- qp.qn.sims %>% 
   ggplot(aes(x = minQP, y = minQN, fill = log(GPP))) +
   geom_tile() + 
   geom_abline(aes(slope = NP_inflow, intercept = 0), col = "white", lty = "dashed") + 
   ggh4x::facet_grid2(model~NP_inflow2) +
   scale_fill_viridis_c())


# plot QP sims
(qn.plot <- qn.sims %>%
    filter(species == "average", Pin == "Pin = 100") %>%
    ggplot() +
    geom_line(aes(x = NPin_molar, y = GPP, col = as.factor(minQN), group = as.factor(minQN)), lwd = 0.75) +
    geom_point(aes(x = NPin_molar, y = GPP, fill = as.factor(minQN), pch = as.factor(minQN), group = as.factor(minQN)), size = 2) +
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~model) + 
    ggh4x::facet_grid2(.~model) + 
    scale_shape_manual(values = c(21, 22, 24, 25)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    #guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    labs(x = "Load N:P (molar)", 
         y = expression("GPP (mg C L"^-1 ~ " day"^-1~")"),
         col = "minQN (mg N mg C^-1)",
         pch = "minQN (mg N mg C^-1)",
         fill = "minQN (mg N mg C^-1)"))

### combine figures
minQP.QN <- ggarrange(plotlist = list(qp.plot, qn.plot), nrow = 1, ncol = 2, 
                      align = "hv", labels = c("a", "b"))

minQP.QN
