
### phytoSTOICH
### GPP and C:N:P sensitivity analysis
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

### Static model
# KP
kp.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        KP1 = c(1, 5, 10, 15, 30, 50))
kp.loads$Nin <- kp.loads$Pin * kp.loads$NP_inflow

kp.static <- lapply(1:nrow(kp.loads), function(i) {
  #print(i)
  static.algae["Pin"] = kp.loads[i, "Pin"]
  static.algae["Nin"] = kp.loads[i, "Nin"]
  static.algae["KP1"] = kp.loads[i, "KP1"]
  y <- c("A1" = 100, "P" = kp.loads[i, "Pin"],"N" = kp.loads[i, "Nin"])
  run <- ode(y, times, parms = static.algae, func = static.stoich)
  return(run[max(times),])
})

kp.static <- as_data_frame(do.call(rbind, kp.static))
kp.static$Pin <- kp.loads$Pin
kp.static$Nin <- kp.loads$Nin
kp.static$NP_inflow <- kp.loads$NP_inflow
kp.static$trait = "KP"
kp.static$trait.value = kp.loads$KP1

# KN
kn.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        KN1 = c(10, 20, 30, 50, 100))
kn.loads$Nin <- kn.loads$Pin * kn.loads$NP_inflow

kn.static <- lapply(1:nrow(kn.loads), function(i) {
  #print(i)
  static.algae["Pin"] = kn.loads[i, "Pin"]
  static.algae["Nin"] = kn.loads[i, "Nin"]
  static.algae["KN1"] = kn.loads[i, "KN1"]
  y <- c("A1" = 100, "P" = kn.loads[i, "Pin"],"N" = kn.loads[i, "Nin"])
  run <- ode(y, times, parms = static.algae, func = static.stoich)
  return(run[max(times),])
})

kn.static <- as_data_frame(do.call(rbind, kn.static))
kn.static$Pin <- kn.loads$Pin
kn.static$Nin <- kn.loads$Nin
kn.static$NP_inflow <- kn.loads$NP_inflow
kn.static$trait = "KN"
kn.static$trait.value = kn.loads$KN1

# QP
qp.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        QP1 = c(0.001, 0.005, 0.01, 0.05, 0.1))
qp.loads$Nin <- qp.loads$Pin * qp.loads$NP_inflow

qp.static <- lapply(1:nrow(qp.loads), function(i) {
  #print(i)
  static.algae["Pin"] = qp.loads[i, "Pin"]
  static.algae["Nin"] = qp.loads[i, "Nin"]
  static.algae["QP1"] = qp.loads[i, "QP1"]
  y <- c("A1" = 100, "P" = qp.loads[i, "Pin"],"N" = qp.loads[i, "Nin"])
  run <- ode(y, times, parms = static.algae, func = static.stoich)
  return(run[max(times),])
})

qp.static <- as_data_frame(do.call(rbind, qp.static))
qp.static$Pin <- qp.loads$Pin
qp.static$Nin <- qp.loads$Nin
qp.static$NP_inflow <- qp.loads$NP_inflow
qp.static$trait = "minQP"
qp.static$trait.value = qp.loads$QP1


# QN
qn.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        QN1 = c(0.01, 0.05, 0.1, 0.5))
qn.loads$Nin <- qn.loads$Pin * qn.loads$NP_inflow

qn.static <- lapply(1:nrow(qn.loads), function(i) {
  #print(i)
  static.algae["Pin"] = qn.loads[i, "Pin"]
  static.algae["Nin"] = qn.loads[i, "Nin"]
  static.algae["QN1"] = qn.loads[i, "QN1"]
  y <- c("A1" = 100, "P" = qn.loads[i, "Pin"],"N" = qn.loads[i, "Nin"])
  run <- ode(y, times, parms = static.algae, func = static.stoich)
  return(run[max(times),])
})

qn.static <- as_data_frame(do.call(rbind, qn.static))
qn.static$Pin <- qn.loads$Pin
qn.static$Nin <- qn.loads$Nin
qn.static$NP_inflow <- qn.loads$NP_inflow
qn.static$trait = "minQN"
qn.static$trait.value = qn.loads$QN1

# rmax
umax.loads <- expand.grid(Pin = seq(20, 500, 100),
                          NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                          umax1 = c(0.1, 0.2, 0.5, 0.8, 1))
umax.loads$Nin <- umax.loads$Pin * umax.loads$NP_inflow

umax.static <- lapply(1:nrow(umax.loads), function(i) {
  #print(i)
  static.algae["Pin"] = umax.loads[i, "Pin"]
  static.algae["Nin"] = umax.loads[i, "Nin"]
  static.algae["umax1"] = umax.loads[i, "umax1"]
  y <- c("A1" = 100, "P" = umax.loads[i, "Pin"],"N" = umax.loads[i, "Nin"])
  run <- ode(y, times, parms = static.algae, func = static.stoich)
  return(run[max(times),])
})

umax.static <- as_data_frame(do.call(rbind, umax.static))
umax.static$Pin <- umax.loads$Pin
umax.static$Nin <- umax.loads$Nin
umax.static$NP_inflow <- umax.loads$NP_inflow
umax.static$trait = "umax"
umax.static$trait.value = umax.loads$umax1


# bind all together
static.trait.var <- bind_rows(kp.static, kn.static, qp.static, qn.static, umax.static)
static.trait.var$model <- "static"

### Dynamic model

# KP
kp.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        KP1 = c(1, 5, 10, 15, 30, 50))
kp.loads$Nin <- kp.loads$Pin * kp.loads$NP_inflow

kp.dynamic <- lapply(1:nrow(kp.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = kp.loads[i, "Pin"]
  dynamic.algae["Nin"] = kp.loads[i, "Nin"]
  dynamic.algae["KP1"] = kp.loads[i, "KP1"]
  y <- c("A1" = 100, "P" = kp.loads[i, "Pin"],"N" = kp.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

kp.dynamic <- as_data_frame(do.call(rbind, kp.dynamic))
kp.dynamic$Pin <- kp.loads$Pin
kp.dynamic$Nin <- kp.loads$Nin
kp.dynamic$NP_inflow <- kp.loads$NP_inflow
kp.dynamic$trait = "KP"
kp.dynamic$trait.value = kp.loads$KP1

# KN
kn.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        KN1 = c(10, 20, 30, 50, 100))
kn.loads$Nin <- kn.loads$Pin * kn.loads$NP_inflow

kn.dynamic <- lapply(1:nrow(kn.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = kn.loads[i, "Pin"]
  dynamic.algae["Nin"] = kn.loads[i, "Nin"]
  dynamic.algae["KN1"] = kn.loads[i, "KN1"]
  y <- c("A1" = 100, "P" = kn.loads[i, "Pin"],"N" = kn.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

kn.dynamic <- as_data_frame(do.call(rbind, kn.dynamic))
kn.dynamic$Pin <- kn.loads$Pin
kn.dynamic$Nin <- kn.loads$Nin
kn.dynamic$NP_inflow <- kn.loads$NP_inflow
kn.dynamic$trait = "KN"
kn.dynamic$trait.value = kn.loads$KN1

# QP
qp.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        minQP1 = c(0.001, 0.005, 0.01, 0.05, 0.1))
qp.loads$Nin <- qp.loads$Pin * qp.loads$NP_inflow

qp.dynamic <- lapply(1:nrow(qp.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = qp.loads[i, "Pin"]
  dynamic.algae["Nin"] = qp.loads[i, "Nin"]
  dynamic.algae["minQP1"] = qp.loads[i, "minQP1"]
  y <- c("A1" = 100, "P" = qp.loads[i, "Pin"],"N" = qp.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

qp.dynamic <- as_data_frame(do.call(rbind, qp.dynamic))
qp.dynamic$Pin <- qp.loads$Pin
qp.dynamic$Nin <- qp.loads$Nin
qp.dynamic$NP_inflow <- qp.loads$NP_inflow
qp.dynamic$trait = "minQP"
qp.dynamic$trait.value = qp.loads$minQP1


# QN
qn.loads <- expand.grid(Pin = seq(20, 500, 100),
                        NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                        minQN1 = c(0.01, 0.05, 0.1, 0.5))
qn.loads$Nin <- qn.loads$Pin * qn.loads$NP_inflow

qn.dynamic <- lapply(1:nrow(qn.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = qn.loads[i, "Pin"]
  dynamic.algae["Nin"] = qn.loads[i, "Nin"]
  dynamic.algae["minQN1"] = qn.loads[i, "minQN1"]
  y <- c("A1" = 100, "P" = qn.loads[i, "Pin"],"N" = qn.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

qn.dynamic <- as_data_frame(do.call(rbind, qn.dynamic))
qn.dynamic$Pin <- qn.loads$Pin
qn.dynamic$Nin <- qn.loads$Nin
qn.dynamic$NP_inflow <- qn.loads$NP_inflow
qn.dynamic$trait = "minQN"
qn.dynamic$trait.value = qn.loads$minQN1

# VmaxP
vmaxP.loads <- expand.grid(Pin = seq(20, 500, 100),
                           NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                           upP1 = c(0.01, 0.05, 0.1, 0.5, 1))
vmaxP.loads$Nin <- vmaxP.loads$Pin * vmaxP.loads$NP_inflow

vmaxP.dynamic <- lapply(1:nrow(vmaxP.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = vmaxP.loads[i, "Pin"]
  dynamic.algae["Nin"] = vmaxP.loads[i, "Nin"]
  dynamic.algae["upP1"] = vmaxP.loads[i, "upP1"]
  y <- c("A1" = 100, "P" = vmaxP.loads[i, "Pin"],"N" = vmaxP.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

vmaxP.dynamic <- as_data_frame(do.call(rbind, vmaxP.dynamic))
vmaxP.dynamic$Pin <- vmaxP.loads$Pin
vmaxP.dynamic$Nin <- vmaxP.loads$Nin
vmaxP.dynamic$NP_inflow <- vmaxP.loads$NP_inflow
vmaxP.dynamic$trait = "VmaxP"
vmaxP.dynamic$trait.value = vmaxP.loads$upP1

# VmaxN
vmaxN.loads <- expand.grid(Pin = seq(20, 500, 100),
                           NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                           upN1 = c(0.01, 0.05, 0.1, 0.5, 1, 5))
vmaxN.loads$Nin <- vmaxN.loads$Pin * vmaxN.loads$NP_inflow

vmaxN.dynamic <- lapply(1:nrow(vmaxN.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = vmaxN.loads[i, "Pin"]
  dynamic.algae["Nin"] = vmaxN.loads[i, "Nin"]
  dynamic.algae["upN1"] = vmaxN.loads[i, "upN1"]
  y <- c("A1" = 100, "P" = vmaxN.loads[i, "Pin"],"N" = vmaxN.loads[i, "Nin"],  "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

vmaxN.dynamic <- as_data_frame(do.call(rbind, vmaxN.dynamic))
vmaxN.dynamic$Pin <- vmaxN.loads$Pin
vmaxN.dynamic$Nin <- vmaxN.loads$Nin
vmaxN.dynamic$NP_inflow <- vmaxN.loads$NP_inflow
vmaxN.dynamic$trait <- "VmaxN"
vmaxN.dynamic$trait.value = vmaxN.loads$upN1

# umax
umax.loads <- expand.grid(Pin = seq(20, 500, 100),
                          NP_inflow = c(1, 2, 3, 7.23, seq(5, 100, 5)), 
                          umax1 = c(0.1, 0.2, 0.5, 0.8, 1))
umax.loads$Nin <- umax.loads$Pin * umax.loads$NP_inflow

umax.dynamic <- lapply(1:nrow(umax.loads), function(i) {
  #print(i)
  dynamic.algae["Pin"] = umax.loads[i, "Pin"]
  dynamic.algae["Nin"] = umax.loads[i, "Nin"]
  dynamic.algae["umax1"] = umax.loads[i, "umax1"]
  y <- c("A1" = 100, "P" = umax.loads[i, "Pin"],"N" = umax.loads[i, "Nin"], "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich)
  return(run[max(times),])
})

umax.dynamic <- as_data_frame(do.call(rbind, umax.dynamic))
umax.dynamic$Pin <- umax.loads$Pin
umax.dynamic$Nin <- umax.loads$Nin
umax.dynamic$NP_inflow <- umax.loads$NP_inflow
umax.dynamic$trait = "umax"
umax.dynamic$trait.value = umax.loads$umax1


# bind all together
dynamic.trait.var <- bind_rows(kp.dynamic, kn.dynamic, qp.dynamic, qn.dynamic, 
                               vmaxP.dynamic, vmaxN.dynamic, umax.dynamic)
dynamic.trait.var$model <- "dynamic"

# combine static and dynamic models
trait.var.combined <- bind_rows(static.trait.var, dynamic.trait.var)
trait.var.combined$model <- factor(trait.var.combined$model, levels = c("static", "dynamic"))
trait.var.combined$trait <- factor(trait.var.combined$trait, 
                                   levels = c("KP", "KN", "minQP", "minQN", "VmaxN", "VmaxP", "umax"))

################################################################################

# Plotting

# plt sensitivity of GPP to changes in traits

(gpp.sensitivity <- trait.var.combined %>%
  ggplot() + 
  geom_path(aes(trait.value, GPP, group = interaction(NP_inflow, Pin)), alpha = 0.3, lwd = .75) + 
  geom_point(aes(trait.value, GPP, col = NP_inflow, group = Pin), size = 2) + 
  ggh4x::facet_grid2(model~trait,  scales = "free", independent = "all") + 
  scale_color_viridis_c() +
  scale_x_log10() + scale_y_log10() + 
  labs(x = "Trait value",
       y =  "GPP mg C L^-1 day^-1",
       col = "N:P inflow mass"))

trait.var.plot


# plt sensitivity of C:N:P to changing traits
# add stoichiometric ratios
trait.var.combined$CN_mass <- 1/trait.var.combined$QN1
trait.var.combined$CP_mass <- 1/trait.var.combined$QP1
trait.var.combined$NP_mass <- trait.var.combined$QN1/trait.var.combined$QP1

# to molar
trait.var.combined$CN_molar <- (trait.var.combined$CN_mass * 14.007)/(1*12.011)
trait.var.combined$CP_molar <- (trait.var.combined$CP_mass * 30.974)/(1*12.011)
trait.var.combined$NP_molar <- (trait.var.combined$QN1/14.007)/(trait.var.combined$QP1/30.974)

# CN
(cn.sensitivity <- trait.var.combined %>%
  filter(model == "dynamic") %>%
  ggplot() + 
  geom_path(aes(trait.value, CN_molar, group = interaction(NP_inflow, Pin)), alpha = 0.3, lwd = .75) + 
  geom_point(aes(trait.value,CN_molar, col = NP_inflow), size = 2) + 
  ggh4x::facet_grid2(.~trait,  scales = "free", independent = "all") + 
  scale_color_viridis_c() + 
  scale_x_log10() + scale_y_log10() + 
  labs(x = NULL,
       y =  "cell C:N molar",
       col = "Inflow N:P mass") + 
  theme(legend.position = "right"))

# C:P
(cp.sensitivity <- trait.var.combined %>%
    filter(model == "dynamic") %>%
    ggplot() + 
    geom_path(aes(trait.value, CP_molar, group = interaction(NP_inflow, Pin)), alpha = 0.3, lwd = .75) + 
    geom_point(aes(trait.value,CP_molar, col = NP_inflow), size = 2) + 
    ggh4x::facet_grid2(.~trait,  scales = "free", independent = "all") + 
    scale_color_viridis_c() + 
    scale_x_log10() + scale_y_log10() + 
    labs(x = NULL,
         y =  "cell C:P molar",
         col = "Inflow N:P mass") + 
    theme(legend.position = "right"))

# N:P
(np.sensitivity <- trait.var.combined %>%
    filter(model == "dynamic") %>%
    ggplot() + 
    geom_path(aes(trait.value, NP_molar, group = interaction(NP_inflow, Pin)), alpha = 0.3, lwd = .75) + 
    geom_point(aes(trait.value,NP_molar, col = NP_inflow), size = 2) + 
    ggh4x::facet_grid2(.~trait,  scales = "free", independent = "all") + 
    scale_color_viridis_c() + 
    scale_x_log10() + scale_y_log10() + 
    labs(x = NULL,
         y =  "cell N:P molar",
         col = "Inflow N:P mass") + 
    theme(legend.position = "right"))

# combine
cnp.sensitivity <- ggarrange(plotlist = list(cn.sensitivity, cp.sensitivity, np.sensitivity), 
                           align = "hv", nrow = 3, labels = c("a", "b", "c"))
cnp.sensitivity
