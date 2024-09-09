
### R-script for dynamic model grid optimization
# DG, September, 2024
# phytoStoich

# load input and validation data set 
source("clean_corman.R")
source("algae_param_vctrs.R")

# load packages
pck <- c("deSolve", "tidyverse")
lapply(pck, require, character.only = T)

# load model and parameterize
source("models/dynamic_liebig_zmix.R") 
times <- 1:1000


# create in-put grid across which model is run
dynamic.grid.search <- expand.grid(
  KP1 = c(1, 5, 10, 15, 20),
  KN1  = c(20, 30, 40, 50, 60),
  minQP1 =  c(0.001, 0.005, 0.01, 0.05, 0.1),
  minQN1 = c(0.01, 0.05, 0.1, 0.15, 0.2),
  upP1 = 0.4255, # c(0.1, 0.3, 0.6), 
  upN1 =  2.453, # c(0.1, 0.5, 1),
  umax1 = c(0.2, 0.5, 1, 1.5)
)

# repeat grid nrow(corman2) times
dynamic.grid.search2 <- dynamic.grid.search %>%
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
dynamic.grid.search2$Lake <- rep(corman2$Lake, nrow(dynamic.grid.search))
dynamic.grid.search2$GPP <- rep(corman2$GPP, nrow(dynamic.grid.search))
dynamic.grid.search2$TN_in <- rep(corman2$TN_in, nrow(dynamic.grid.search))
dynamic.grid.search2$TP_in <- rep(corman2$TP_in, nrow(dynamic.grid.search))
dynamic.grid.search2$HRT <- rep(corman2$HRT, nrow(dynamic.grid.search))
dynamic.grid.search2$SA <- rep(corman2$SA, nrow(dynamic.grid.search))
dynamic.grid.search2$DOC_mgL <- rep(corman2$DOC_mgL, nrow(dynamic.grid.search))
dynamic.grid.search2$z <- rep(corman2$z, nrow(dynamic.grid.search))

# run model
(start <- Sys.time())
dynamic.grid.search.out <- lapply(1:nrow(dynamic.grid.search2), function(i){
  #print(i)
  # indexing
  dynamic.algae["SA"] = dynamic.grid.search2[i, "SA"]
  dynamic.algae["DOC"] = dynamic.grid.search2[i, "DOC_mgL"]
  dynamic.algae["z"] = dynamic.grid.search2[i, "z"]
  dynamic.algae["KP1"] = dynamic.grid.search2[i, "KP1"]
  dynamic.algae["KN1"] = dynamic.grid.search2[i, "KN1"]
  dynamic.algae["minQP1"] = dynamic.grid.search2[i, "minQP1"]
  dynamic.algae["minQN1"] = dynamic.grid.search2[i, "minQN1"]
  dynamic.algae["upP1"] = dynamic.grid.search2[i, "upP1"]
  dynamic.algae["upN1"] = dynamic.grid.search2[i, "upN1"]
  dynamic.algae["umax1"] = dynamic.grid.search2[i, "umax1"]
  dynamic.algae["Pin"] = dynamic.grid.search2[i, "TP_in"]
  dynamic.algae["Nin"] = dynamic.grid.search2[i, "TN_in"]
  dynamic.algae["HRT"] = dynamic.grid.search2[i, "HRT"]
  # starting values
  y <- c("A1" = 100, "P" = dynamic.grid.search2[i, "TP_in"], "N" = dynamic.grid.search2[i, "TN_in"],
         "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich.zmix)
  return(run[max(times),])
})
# convert to df
dynamic.grid.search.out <- do.call(rbind, dynamic.grid.search.out)
dynamic.grid.search.out <- as_data_frame(dynamic.grid.search.out)
(end <- Sys.time())
time.elapsed <- (end - start)
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# output data frame
dynamic.grid.opt <- dynamic.grid.search2
dynamic.grid.opt$est.GPP <- dynamic.grid.search.out$GPP
dynamic.grid.opt$est.A <- dynamic.grid.search.out$A1
dynamic.grid.opt$est.P <- dynamic.grid.search.out$P
dynamic.grid.opt$est.N <- dynamic.grid.search.out$N

# calculate error metrics to find the "best" model runs
dynamic.grid.opt.metrics <- dynamic.grid.opt %>%
  group_by(minQP1, minQN1, umax1, KN1, KP1, upP1, upN1) %>%
  summarise(MAE = mean(abs(GPP - est.GPP)), # lower is better
            RMSE =  sqrt((mean(GPP - est.GPP)^2)), # lower is better
            NSE = 1 - sum((GPP - est.GPP)^2) / sum((GPP - mean(GPP))^2)  # closer to one is better; below 0 is not good
  )  %>%
  mutate(NSE = 1/(2 - NSE)) # normalize NSE to 0 to 1

# save files
#save(dynamic.grid.opt, dynamic.grid.opt.metrics, file =  "corman_dynamic_grid_optim.Rdata")

# quick overview of MAE, RMSE, NSE
hist(dynamic.grid.opt.metrics$MAE)
hist(dynamic.grid.opt.metrics$RMSE)
hist(dynamic.grid.opt.metrics$NSE)
