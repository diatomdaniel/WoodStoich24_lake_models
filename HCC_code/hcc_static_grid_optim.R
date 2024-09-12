
### R-script for static model grid optimization
# DG, September, 2024
# phytoStoich

# load input and validation data set 
source("hcc_clean_corman.R")
source("hcc_algae_param_vctrs.R")

# load packages
pck <- c("deSolve", "tidyverse")
lapply(pck, require, character.only = T)

# load model and parameterize
source("hcc_static_liebig_zmix.R") 
times <- 1:1000

# create in-put grid across which model is run
static.grid.search <- expand.grid(
  KP1 =  seq(1, 20, 1),
  KN1  = seq(5, 50, 5),
  minQP1 =  seq(0.01, 1, 0.01),
  minQN1 = seq(0.01, 1, 0.01),
  umax1 = seq(0.1, 1.5, 0.1)
)

# repeat grid nrow(corman2) times
static.grid.search2 <- static.grid.search %>%
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
static.grid.search2$Lake <- rep(corman2$Lake, nrow(static.grid.search))
static.grid.search2$GPP <- rep(corman2$GPP, nrow(static.grid.search))
static.grid.search2$TN_in <- rep(corman2$TN_in, nrow(static.grid.search))
static.grid.search2$TP_in <- rep(corman2$TP_in, nrow(static.grid.search))
static.grid.search2$HRT <- rep(corman2$HRT, nrow(static.grid.search))
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
  static.algae["HRT"] = static.grid.search2[i, "HRT"]
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
save(static.grid.opt, static.grid.opt.metrics,  file = "hcc_corman_static_grid_optim.Rdata")

# quick overview of MAE, RMSE, NSE
# hist(static.grid.opt.metrics$MAE)
# hist(static.grid.opt.metrics$RMSE)
# hist(static.grid.opt.metrics$NSE)
