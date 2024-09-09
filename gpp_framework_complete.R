
### phytoSTOICH
### GPP and algae stoichiometry framework simulations
### DG, September 2024


###############################################################################
# Setup
#load packages
pck <- c("deSolve", "tidyverse", "cowplot", "ggthemes","ggpubr")
lapply(pck, require, character.only = T)

# Load algae parameters
source("algae_param_vctrs.R")

# Load models
# static model with fixed stoichiometry
source("models/static_liebig_no_light.R") # model 1 in Carly's framework
# Droop model
source("models/dynamic_liebig_no_light.R") # model 3 in Carly's framework
# set timesteps
times <- 1:1000 # for troubleshooting, initial runs


################################################################################
# Simulations

# set loads
# low, mid, high P supply --> trophic state
loads <- expand.grid(Pin = c(50, 100, 500),
                     NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)))
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
         Nin = rep(loads$Nin, 8), 
         NP_inflow = rep(loads$NP_inflow, 8)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")),
         species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")),
         Pin = factor(Pin))


################################################################################
# Plotting

# subset data for easy plotting
average <- gpp.sims[gpp.sims$species == "average",]
diatoms <- gpp.sims[gpp.sims$species == "diatoms",]
greens <- gpp.sims[gpp.sims$species == "greens",]
cyanos <- gpp.sims[gpp.sims$species == "cyanos",]

# data set for seston
seston <- gpp.sims  %>%
  mutate("C:N" = 1/QN1,"C:P" = 1/QP1, "N:P" = QN1/QP1, NP_inflow = NP_inflow) %>%
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
  "minQN_minQP" = rep(c(0.09/0.0105, 0.155/0.02, 0.025/0.001, 0.01/0.001), 2))

# create seston summary w. minQN and minQP
consump.vctr.seston <- seston %>%
  merge(consump.vctr, by = "species") %>%
  group_by(species, seston) %>%
  summarise(minQN_minQP = max(minQN_minQP), 
            CNP = max(value),
            GPP = max(GPP))

# manual legend
species.legend = c("average", "diatoms","greens", "cyanos")

# plot GPP across supply N:P for each species
(gpp.plt <- ggplot() + 
  # # add in the consumption vctrs ala Tilmann
  # note that consumption vctrs won't match for the static model, only the dynamic as cell quotas change..
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species), 
               lwd = 1, lty = "dashed") +
  geom_segment(data = consump.vctr.seston,
               aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species), 
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
  geom_line(data = gpp.sims, aes(x = NP_inflow, y = GPP,
                                 col = species, group = interaction(species, Pin)),
            lwd = 0.75) +
  geom_point(data = gpp.sims, aes(x = NP_inflow, y = GPP,
                                 fill = species, pch = Pin, group = Pin),
             size = 2) + 
  scale_x_log10() + scale_y_log10() + 
  #ggh4x::facet_grid2(.~model) + 
  ggh4x::facet_grid2(Pin~model) + 
  scale_shape_manual(values = c(21, 22, 24)) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + 
  guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
  labs(x = "N:P inflow mass", y = "GPP mg C L^-1 day^-1", col = "P in ug L^-1"))


# plot CNP across supply N:P for each species
(cnp.plt <- ggplot() + 
    # add in the consumption vctrs ala Tilmann
    geom_segment(data = consump.vctr.seston, 
                 aes(y = CNP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species), 
                 lwd = 1, lty = "dashed") + 
    geom_segment(data = consump.vctr.seston, 
                 aes(y = CNP, x = minQN_minQP, xend = 0, yend = CNP, col = species), 
                 lwd = 1, lty = "dashed") + 
    # data
    # geom_line(data = seston.average, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.average, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 21, size = 2) + 
    # geom_line(data = seston.diatoms, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.diatoms, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 22, size = 2) + 
    # geom_line(data = seston.greens, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.greens, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 23, size = 2) + 
    # geom_line(data = seston.cyanos, aes(NP_inflow, value, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.cyanos, aes(NP_inflow, value, fill = Pin, group = Pin), pch = 24, size = 2) +
    geom_line(data = seston, aes(NP_inflow, value, col = species, group = interaction(species, Pin)), 
              lwd = 0.75) + 
    geom_point(data = seston, aes(NP_inflow, value, fill = species, pch = Pin, group = Pin), size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~seston, scales = "free", independent = "y") + 
    ggh4x::facet_grid2(Pin~seston, scales = "free", independent = "y") + 
    scale_shape_manual(values = c(21, 22, 23)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    labs(x = "N:P inflow mass", y = "C:N:P", col = "P in ug L^-1"))


# plot GPP across CNP for each species
# quick edit to consumption vctr.
consump.vctr.seston <- consump.vctr.seston %>% mutate(minQN_minQP = ifelse(seston == "N:P", minQN_minQP, NA))
(gpp.cnp.plt <- ggplot() + 
    # # add in the consumption vctrs ala Tilmann
    geom_segment(data = consump.vctr.seston,
                 aes(y = GPP, x = minQN_minQP, xend = minQN_minQP, yend = 0, col = species), 
                 lwd = 0.75, lty = "dashed") +
    geom_segment(data = consump.vctr.seston,
                 aes(y = GPP, x = minQN_minQP, xend = 0, yend = GPP, col = species), 
                 lwd = 0.75, lty = "dashed") +
    # data
    # geom_line(data = seston.average, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.average, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 21, size = 2) + 
    # geom_line(data = seston.diatoms, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.diatoms, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 22, size = 2) + 
    # geom_line(data = seston.greens, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.greens, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 23, size = 2) + 
    # geom_line(data = seston.cyanos, aes(value, GPP, col = Pin, group = Pin), lwd = 1) + 
    # geom_point(data = seston.cyanos, aes(value, GPP, fill = Pin, group = Pin, alpha = log(NP_inflow)), pch = 24, size = 2) + 
    geom_line(data = seston, aes(value, GPP, col = species, alpha = NP_inflow, group = interaction(species, Pin)),
              lwd = 0.75) + 
    geom_point(data = seston, aes(value, GPP, fill = species, pch = Pin, alpha = NP_inflow, group = interaction(species, Pin)), 
               size = 2) + 
    scale_x_log10() + scale_y_log10() + 
    #ggh4x::facet_grid2(.~seston) + 
    ggh4x::facet_grid2(Pin~seston) + 
    scale_shape_manual(values = c(21, 22, 24)) + 
    scale_color_viridis_d() + 
    scale_fill_viridis_d() + 
    guides(fill = "none", color = "none", pch = "none", alpha = "none") + 
    scale_alpha(range=c(0.5,1), na.value = 0) + 
    labs(x = "C:N:P molar", y = "GPP mg C L^-1 day_1", col = "P in ug L^-1", size = "Inflow log(N:P) mass"))

################################################################################

# arrange figures
# seston plots go together
seston.figs <- ggarrange(plotlist = list(cnp.plt, gpp.cnp.plt), align = "hv", ncol = 2, labels = c("b", "c"))
seston.figs

# add to gpp plot
#gpp.plt2 <- ggarrange(plotlist = list(gpp.plt, ggplot() + geom_blank()), nrow = 2, labels = c("a", ""))

framework.fig <- ggarrange(plotlist = list(gpp.plt, seston.figs),ncol = 2, labels = c("a", ""), 
                           align = "hv", widths = c(1/3, 2/3))
framework.fig
# manually add legend in powerpoint