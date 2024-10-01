### phytoSTOICH
### GPP validation
### DG, September 2024

# set wd
rm(list = ls())
setwd("C:/Users/DanielGschwentner/Documents/GitHub/WoodStoich24_lake_models")
###############################################################################
# Setup
# load Corman data (needs to be run first!
source("clean_corman.R")
# to df for data input
in.grid <- as.data.frame(corman2)
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
times <- 1:1000 # for troubleshooting, initial runs

################################################################################

### QN:QP gird search
cell.quota.grid <- expand.grid(
  minQP1 = c(0.001, 0.005, 0.01, 0.05, 0.1),
  cell.NP = c(1, 3, 5, 7, 10, 15, 30)
)
cell.quota.grid$minQN1 <- cell.quota.grid$minQP1 * cell.quota.grid$cell.NP

# combine grid search with corman et al data
cell.quota.grid.ext <- cell.quota.grid %>% slice(rep(1:n(), each = nrow(corman2)))

# extend corman data set
corman2.ext <- corman2 %>% slice(rep(1:n(), each = nrow(cell.quota.grid)))

# bind together
in.grid <- bind_cols(corman2.ext, cell.quota.grid.ext)

################################################################################

### base predictions w. median values

# static model
(start <- Sys.time())
corman.static <-  lapply(list(static.algae), function(x) {
  params <- x
  lapply(1:nrow(in.grid), function(i) {
    print(i)
    params["Pin"] = in.grid[i, "TP_in"]
    params["Nin"] = in.grid[i, "TN_in"]
    params["DOC"] = in.grid[i, "DOC_mgL"]
    params["z"] = in.grid[i, "z"]
    params["SA"] = in.grid[i, "SA"]
    params["HRT"] = in.grid[i, "HRT"]
    params["QP1"] = in.grid[i, "minQP1"]
    params["QN1"] = in.grid[i, "minQN1"]
    y <- c("A1" = 100, "P" =in.grid[i, "TP_in"], "N" = in.grid[i, "TN_in"])
    run <- ode(y, times, parms = params, func = static.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# extract from list and convert to data-frame
static <- as_data_frame(do.call(rbind, corman.static[[1]]))

# dynamic model
(start <- Sys.time())
corman.dynamic <-  lapply(list(dynamic.algae, dynamic.diatoms, dynamic.greens, dynamic.cyanos), function(x) {
  params <- x
  lapply(1:nrow(in.grid), function(i) {
    params["Pin"] = in.grid[i, "TP_in"]
    params["Nin"] = in.grid[i, "TN_in"]
    params["DOC"] = in.grid[i, "DOC_mgL"]
    params["z"] = in.grid[i, "z"]
    params["SA"] = in.grid[i, "SA"]
    params["HRT"] = in.grid[i, "HRT"]
    params["QP1"] = in.grid[i, "minQP1"]
    params["QN1"] = in.grid[i, "minQN1"]
    y <- c("A1" = 100, "P" =in.grid[i, "TP_in"], "N" = in.grid[i, "TN_in"], 
           "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

dynamic  <- as_data_frame(do.call(rbind, corman.dynamic[[1]]))

################################################################################

### Add back into original data frame

in.grid$static.GPP <- static.average$GPP
in.grid$dynamic.GPP <- dynamic.average$GPP

x<- in.grid %>% group_by(minQP1, cell.NP, Lake) %>%
  summarise(r2.static = cor(static.GPP, GPP)^2) #
            #r2.dynamic = cor(dynamic.GPP, obs.GPP))
  
x %>% ggplot(aes(cell.NP, r2.static, col = factor(minQP1), group = factor(minQP1))) + 
  geom_point() + geom_path() + facet_wrap(Lake~.)

x %>% ggplot(aes(cell.NP, r2.dynamic, col = factor(minQP1), group = factor(minQP1))) + 
  geom_point() + geom_path() + facet_wrap(Lake~.)
