
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
# michaelis menten grid
grid <- expand.grid(
  KP1 = seq(1, 20, 0.5 ),
  KN1 = seq(1, 50, 1), 
  minQN1 = seq(0.1, 5, 0.1),
  minQP1 = seq(0.1, 5, 0.1)
)

# repeat grid nrow(corman2) times
grid.new <- grid %>% 
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
grid.new$Lake <- rep(corman2$Lake, nrow(grid))
grid.new$GPP <- rep(corman2$GPP, nrow(grid))
grid.new$TN_load <- rep(corman2$TN_load, nrow(grid))
grid.new$TP_load <- rep(corman2$TP_load, nrow(grid))
grid.new$SA <- rep(corman2$SA, nrow(grid))
grid.new$DOC_mgL <- rep(corman2$DOC_mgL, nrow(grid))
grid.new$z_m <- rep(corman2$z_m, nrow(grid))

##### Model parameterization

## Model structure
# easier to copy and paste model here than source file externally
mich.single <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  P <- y["P"]
  N <- y["N"]
  
  # lake parameters
  SA= params["SA"]
  z <- params["z"]
  DOC = params["DOC"]
  Pin = params["Pin"]
  Nin = params["Nin"]
  
  # algae physiology parameters
  umax1 = params["umax1"]
  lA = params["lA"]
  v = params["v"]
  
  # species 1
  KP1 = params["KP1"]
  QP1 = params["QP1"]
  KN1 = params["KN1"]
  QN1 = params["QN1"]
  Klight = params["KLight"]
  
  # zmix
  zmix <- 10^(0.515 * log10(DOC) + 0.115 * log10(2 * sqrt(SA/pi + 0.991)))
  # if zmix exceeds depth, set zmix to z.
  zmix <- ifelse(zmix > z, z, zmix)
  # In/output
  Qin=SA*1e6*zmix/365	# m^3 day^-1
  
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * zmix
  # biomass specific growth for entire mixed layer
  # species 1
  prod1= (umax1) * (P/(P+KP1))*(N/(N + KN1))	# d-1
  GPP = prod1 * A1/1000 # this is the GPP rate! mg C L^-1 day^-1
  
  # model biomass
  # species 1
  dA1.dt=A1*prod1-lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# total biomass in mg C m^-3
  
  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P)+QP1*lA*A1-QP1*A1*prod1 #mg P m-3 (in epi);
  
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N) + QN1*lA*A1-QN1*A1*prod1  # mg N m-3 (in epi);
  
  # # indicators of limitation
  # Plim <- 1 - (P/(P + KP)) # unitless
  # Nlim <- 1 - (N/(N + KN))
  # Llim <- 1 - (1/(kD * zmix)) * log((KLight + I0)/(KLight + Izmix)) # unitless
  # NP_moles <- (N/14.007)/(P/30.974)
  
  # return objects
  dY=c(A1= dA1.dt, dPdt=dP.dt, dNdt = dN.dt)
  gpp = c(GPP = GPP)
  names(gpp) = "GPP"
  #lim=c(Plim = Plim, Nlim = Nlim, Llim = Llim, NP_moles = NP_moles, growth.rate = prod)
  return(list(dY, gpp))
  
}

# nr. of steps
times <- 1:10000

## Set up model parameters
## michaelis-menten model for mean phytoplankton traits
params.mich <- c(
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
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = NA, # phosphorus half sat constant in mg P m^-3 f
  QP1 = NA, # algae cell P quota in mg P mg^-1 C^-1
  KN1 = NA, # nitrogen half sat constant in mg N m^-3 
  QN1 = NA # algae cell N quota in mg N mg^-1 C^-1 
)
names(params.mich) <- c("SA", "Pin", "Nin", "DOC", "z", "umax1", "lA", "v", 
                        "KP1",  "QP1", "KN1",  "QN1")
names(params.mich)
params.mich <- unlist(params.mich)

##### Perform grid search and return output

#takes over an hour to run!!
(start <- Sys.time())
mm.grid <- lapply(1:nrow(grid.new), function(i){
  # indexing
  params.mich["SA"] = grid.new[i, "SA"]
  params.mich["DOC"] = grid.new[i, "DOC_mgL"]
  params.mich["z"] = grid.new[i, "z_m"]
  params.mich["KP1"] = grid.new[i, "KP1"]
  params.mich["KN1"] = grid.new[i, "KN1"]
  params.mich["QP1"] = grid.new[i, "minQP1"]
  params.mich["QN1"] = grid.new[i, "minQN1"]
  params.mich["Pin"] = grid.new[i, "TP_load"]
  params.mich["Nin"] = grid.new[i, "TN_load"]
  # starting values
  y <- c("A1" = 100, "P" = grid.new[i, "TP_load"], "N" = grid.new[i, "TN_load"])
  run <- ode(y, times, parms = params.mich, func = mich.single)
  return(run[max(times),])
})
# convert to df
mm.grid <- do.call(rbind, mm.grid)
mm.grid <- as_data_frame(mm.grid)
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))
# add to grid
grid.new$est_GPP <- mm.grid$GPP
# # grids take too long to run, gonna save to r data file and reload for convenience
save(mm.grid,file = "corman2023_MM_gridsearch_raw.Rdata")
save(grid.new, file = "corman2023_MM_gridsearch_matched.Rdata")
