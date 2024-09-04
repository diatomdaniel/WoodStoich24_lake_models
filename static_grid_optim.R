
### R-script for static model grid optimization
# DG, September, 2024
# phytoStoich

# load input and validation data set 
source("clean_corman.R")

# load packages
pck <- c("deSolve", "tidyverse", "cowplot", "ggthemes", "wesandersen", "ggpubr")
lapply(pck, require, character.only = T)

# load model and parameterize
source("models/static_liebig_zmix.R") 
times <- 1:1000

#parameterize models
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
  KP1 = NA, # phosphorus half sat constant in mg P m^-3 f
  QP1 = NA, # algae cell P quota in mg P mg^-1 C^-1
  KN1 = NA, # nitrogen half sat constant in mg N m^-3 
  QN1 = NA # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.algae) <- c("SA", "Pin", "Nin", "DOC", "z", "umax1", "lA", "v", 
                         "KP1",  "QP1", "KN1",  "QN1")
names(static.algae)
static.algae <- unlist(static.algae)

# create in-put grid across which model is run
static.grid.search <- expand.grid(
  KP1 = c(1, 5, 10, 15, 20),
  KN1  = c(20, 30, 40, 50, 60),
  minQP1 =  c(0.001, 0.005, 0.01, 0.05, 0.1),
  minQN1 = c(0.01, 0.05, 0.1, 0.15, 0.2),
  umax1 = c(0.2, 0.5, 1, 1.5)
)

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

 # run model
(start <- Sys.time())
static.grid.search.out <- lapply(1:nrow(static.grid.search2), function(i){
  #print(i)
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

# output data frame
static.grid.opt <- static.grid.search2
static.grid.opt$est.GPP <- static.grid.search.out$GPP
static.grid.opt$est.A <- static.grid.search.out$A1
static.grid.opt$est.P <- static.grid.search.out$P
static.grid.opt$est.N <- static.grid.search.out$N

# calculate error metrics to find the "best" model runs
static.grid.opt.metrics <- static.grid.opt %>%
  group_by(minQP1, minQN1, umax1, KN1, KP1) %>%
  summarise(MAE = mean(abs(GPP - est.GPP)), # lower is better
            RMSE =  sqrt((mean(GPP - est.GPP)^2)), # lower is better
            NSE = 1 - sum((GPP - est.GPP)^2) / sum((GPP - mean(GPP))^2)  # closer to one is better; below 0 is not good
  )  %>%
  mutate(NSE = 1/(2 - NSE)) # normalize NSE to 0 to 1

# save files
save(static.grid.opt, static.grid.opt.metrics,  file = "corman_static_grid_optim.Rdata")

