
### HCC static model grid search

## Packages
library(tidyverse)
library(deSolve)


# call data/functions
source("clean_corman.R")
source("static_liebig_zmix.R")

# paramterize model
times <- 1:1000
static.algae <- c(
  # lake parameters
  SA= NA,		# lake surface area in km2
  # varies by simulation
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  DOC = NA, 
  z = NA,
  
  # algae physiology parameters
  umax1 = NA,
  lA=0.1,			# mortality rate day-1
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = 0.0165 * 1000, # phosphorus half sat constant in mg P m^-3 f
  QP1 =0.0105, # algae cell P quota in mg P mg^-1 C^-1
  KN1 = 0.05 * 1000, # nitrogen half sat constant in mg N m^-3 
  QN1 = 0.09 # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.algae) <- c("SA", "Pin", "Nin", "DOC", "z", "umax1", "lA", "v", 
                         "KP1",  "QP1", "KN1",  "QN1")
names(static.algae)
static.algae <- unlist(static.algae)

# create input grid
static.grid.search <- expand.grid(
  minQP1 =  seq(0.01, 0.2, 0.01),
  minQN1 = seq(0.01, 0.2, 0.01),
  umax1 = seq(0.1, 1, 0.1 )
)
static.grid.search$KP1 = 16.5
static.grid.search$KN1 = 50

# repeat grid nrow(corman2) times
static.grid.search2 <- static.grid.search %>% 
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
static.grid.search2$Lake <- rep(corman2$Lake, nrow(static.grid.search))
static.grid.search2$GPP <- rep(corman2$GPP, nrow(static.grid.search))
static.grid.search2$TN_in <- rep(corman2$TN_in, nrow(static.grid.search))
static.grid.search2$TP_in <- rep(corman2$TP_in, nrow(static.grid.search))
static.grid.search2$SA <- rep(corman2$SA, nrow(static.grid.search))
static.grid.search2$DOC_mgL <- rep(corman2$DOC_mgL, nrow(static.grid.search))
static.grid.search2$z <- rep(corman2$z, nrow(static.grid.search))

#this would take over 10 hrs to run....maybe best to outsource to HCC??? but can do this later
(start <- Sys.time())
static.grid.search.out <- lapply(1:nrow(static.grid.search2), function(i){
  print(i)
  # indexing
  static.algae["SA"] = static.grid.search2[i, "SA"]
  static.algae["DOC"] = static.grid.search2[i, "DOC_mgL"]
  static.algae["z"] = static.grid.search2[i, "z"]
  static.algae["KP1"] = static.grid.search2[i, "KP1"]
  static.algae["KN1"] = static.grid.search2[i, "KN1"]
  static.algae["QP1"] = static.grid.search2[i, "minQP1"]
  static.algae["QN1"] = static.grid.search2[i, "minQN1"]
  static.algae["umax1"] = static.grid.search2[i, "umax1"]
  static.algae["Pin"] = static.grid.search2[i, "TP_in"]
  static.algae["Nin"] = static.grid.search2[i, "TN_in"]
  # starting values
  y <- c("A1" = 100, "P" = static.grid.search2[i, "TP_in"], "N" = static.grid.search2[i, "TN_in"])
  run <- ode(y, times, parms = static.algae, func = static.stoich.zmix)
  return(run[max(times),])
})
# convert to df
static.grid.search.out <- do.call(rbind, static.grid.search.out)
static.grid.search.out <- as_data_frame(static.grid.search.out)
(end <- Sys.time())
time.elapsed <- (end - start)
print(paste0("Time elapsed = ", time.elapsed, " hours!"))
# add to grid
static.grid.search2$est_GPP <- static.grid.search.out$GPP

# save
saveRDS(static.grid.search.out, "static_liebig_grid_search_out.RDS", static.grid.search.out)