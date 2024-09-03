################################################################################
#
# Grid-search for phytoplankton trait parameters to run on UNL HCC
# Michaelis-Menten model
# Corman et al 2023 data set grid search
# 2024-08-12
#
################################################################################

##### Set up

## Packages
library(deSolve)
library(tidyverse)

## Load data
corman2 <- read_csv("corman2023_model_input.csv")

##### Grid search parameters

## Create grid of scenarios for simulation
# droop grid
grid.droop <- expand.grid(
  KP1 = seq(1, 20, 0.5 ),
  KN1 = seq(1, 50, 1), 
  minQN1 = seq(0.1, 5, 0.1),
  minQP1 = seq(0.1, 5, 0.1),
  upN1 = c(0.1, 0.5, 1),
  upP1 = c(0.1, 0.5, 1)
)

# repeat grid nrow(corman2) times
grid.droop.new <- grid.droop %>% 
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
grid.droop.new$Lake <- rep(corman2$Lake, nrow(grid.droop))
grid.droop.new$GPP <- rep(corman2$GPP, nrow(grid.droop))
grid.droop.new$TN_load <- rep(corman2$TN_load, nrow(grid.droop))
grid.droop.new$TP_load <- rep(corman2$TP_load, nrow(grid.droop))
grid.droop.new$SA <- rep(corman2$SA, nrow(grid.droop))
grid.droop.new$DOC_mgL <- rep(corman2$DOC_mgL, nrow(grid.droop))
grid.droop.new$z <- rep(corman2$z, nrow(grid.droop))
grid.droop.new$z_mix <- rep(corman2$zmix, nrow(grid.droop))

#### Model parameterization

## Model structure
# easier to copy and paste model here than source file externally
droop.single <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  P <- y["P"]
  N <- y["N"]
  QP1 <- y["QP1"]
  QN1 <- y["QN1"]
  
  # lake parameters
  SA= params["SA"]
  z = params["z"]
  DOC = params["DOC"]
  Pin = params["Pin"]
  Nin = params["Nin"]
  
  # algae physiology parameters
  umax1 = params["umax1"]
  lA = params["lA"]
  v = params["v"]
  
  # half sat. constant P
  KP1 = params["KP1"]
  # min cell quota P
  minQP1 = params["minQP1"]
  # uptake rate P
  upP1 = params["upP1"]
  # half sat constant N
  KN1 = params["KN1"]
  # min cell quota N
  minQN1 = params["minQN1"]
  # uptake rate n
  upN1 = params["upN1"]
  
  # zmix
  zmix <- 10^(-0.515 + log10(DOC) + 0.115 * log10(2 * sqrt(SA/pi + 0.991)))
  # if zmix exceeds depth, set zmix to z.
  zmix <- ifelse(zmix > z, z, zmix)
  
  # In/output
  Qin=SA*1e6*zmix/365	# m^3 day^-1
  
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * z
  
  # biomass specific growth for entire mixed layer
  prod1 = (umax1 * min(1 - minQN1/QN1, 1 - minQP1/QP1 ))	# d-1
  GPP = prod1 * A1/1000 # this is the GPP rate! mg C L^-1 day^-1
  
  # model biomass
  dA1.dt=A1*prod1-lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# mg C m-3
  
  # cell quota P  
  dQP1.dt = upP1 * (P/(KP1 + P)) - prod1  * QP1
  # cell quota N  
  dQN1.dt = upN1 * (N/(KN1 + N)) - prod1  * QN1
  
  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P) + A1 * (-upP1 * (P/(KP1 + P))  + lA * QP1)  # mg P m-3 (in epi);
  
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N)+ A1 * (-upN1 * (N/(KN1 + N)) +  lA * QN1)  # mg N m-3 (in epi);
  
  # return objects 
  dY=c(d.GPP = GPP, dPdt=dP.dt, dNdt = dN.dt, dQP1dt = dQP1.dt, dQ1Ndt = dQN1.dt)
  gpp = c(GPP = GPP)
  names(gpp) = "GPP"
  #lim=c(Plim = Plim, Nlim = Nlim, Llim = Llim, NP_moles = NP_moles, growth.rate = prod)
  return(list(dY, gpp))
  
} 

# nr. of steps
times <- 1:10000

## Set up model parameters
## droop model for mean phytoplankton traits

params.droop <- c(
  # lake parameters
  SA= NA,		# lake surface area in km2
  #zmix = 2, # lake mixing depth in m
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  DOC = NA,
  z = NA,
  
  # algae physiology parameters
  umax1 = 0.665,
  lA=0.1,			# mortality rate day-1
  v= 0.1,			# m d-1; sinking loss of algae
  KP1 = NA, # phosphorus half sat constant in mg P m^-3 f
  minQP1 = NA, # algae cell P quota in mg P mg^-1 C^-1 f
  upP1 = NA, # max uptake rate P per day in mg P mg C^-1 day^-1 
  KN1 = NA, # nitrogen half sat constant in mg N m^-3 
  minQN1 = NA, # algae cell N quota in mg N mg^-1 C^-1 f
  upN1 = NA # max uptake rate N per day in mg N mg C^-1 day^-1 
)



names(params.droop) <- c("SA", "Pin", "Nin", "DOC", "z", "umax1", "lA", "v", "KP1",
                         "minQP1", "upP1", "KN1","minQN1", "upN1")
names(params.droop)
params.droop <- unlist(params.droop)

##### Perform grid search and return output

#takes over an hour to run!!
(start <- Sys.time())
droop.grid <- lapply(1:10, function(i){
  # indexing
  print(i)
  params.droop["SA"] = grid.droop.new[i, "SA"]
  params.droop["DOC"] = grid.droop.new[i, "DOC_mgL"]
  params.droop["z"] = grid.droop.new[i, "z"]
  params.droop["KP1"] = grid.droop.new[i, "KP1"]
  params.droop["KN1"] = grid.droop.new[i, "KN1"]
  params.droop["minQP1"] = grid.droop.new[i, "minQP1"]
  params.droop["minQN1"] = grid.droop.new[i, "minQN1"]
  params.droop["upP1"] = grid.droop.new[i, "upP1"]
  params.droop["upN1"] = grid.droop.new[i, "upN1"]
  params.droop["Pin"] = grid.droop.new[i, "TP"]
  params.droop["Nin"] = grid.droop.new[i, "TN"]
  # starting values
  y <- c("A1" = 100, "P" = grid.droop.new[i, "TP_load"],
         "N" = grid.droop.new[i, "TN_load"],
         "QP1" = 0.015,
         "QN1" = 0.1)
  run <- ode(y, times, parms = params.droop, func = droop.single)
  return(run[max(times),])
})
# convert to df
droop.grid <- do.call(rbind, droop.grid)
droop.grid <- as_data_frame(droop.grid)
(end <- Sys.time())
(time.elapsed <- (end - start))
## add to grid
grid.droop.new$est_GPP <- droop.grid$GPP
# save
save(droop.grid,file = "corman2023_Droop_gridsearch_raw.Rdata")
save(grid.droop.new, file = "corman2023_Droop_gridsearch_matched.Rdata")