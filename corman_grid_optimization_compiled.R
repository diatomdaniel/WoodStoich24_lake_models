
### phytoSTOICH
### Corman et al GPP grid search
### DG, September 2024


###############################################################################
# Setup
# load Corman data (needs to be run first!
source("clean_corman.R")
# to df for data input
in.grid <- as.data.frame(corman2)
#load packages
pck <- c("deSolve", "tidyverse", "cowplot", "ggthemes","ggpubr")
lapply(pck, require, character.only = T)

# Load algae parameters
source("algae_param_vctrs.R")

# Load models
# static model with fixed stoichiometry
source("models/static_liebig_zmix.R") # model 1 in Carly's framework
# Droop model
source("models/dynamic_liebig_zmix.R") # model 3 in Carly's framework
# set timesteps
times <- 1:1000 # for troubleshooting, initial runs

################################################################################

### base predictions w. median values

# static model
(start <- Sys.time())
corman.static <-  lapply(list(static.algae, static.diatoms, static.greens, static.cyanos), function(x) {
  params <- x
  lapply(1:nrow(in.grid), function(i) {
    params["Pin"] = in.grid[i, "TP_in"]
    params["Nin"] = in.grid[i, "TN_in"]
    params["DOC"] = in.grid[i, "DOC_mgL"]
    params["z"] = in.grid[i, "z"]
    params["SA"] = in.grid[i, "SA"]
    params["HRT"] = in.grid[i, "HRT"]
    y <- c("A1" = 100, "P" =in.grid[i, "TP_in"], "N" = in.grid[i, "TN_in"])
    run <- ode(y, times, parms = params, func = static.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# extract from list and convert to data-frame
corman.static.average <- as_data_frame(do.call(rbind, corman.static[[1]]))
corman.static.diatoms <- as_data_frame(do.call(rbind, corman.static[[2]]))
corman.static.greens <- as_data_frame(do.call(rbind, corman.static[[3]]))
corman.static.cyanos <- as_data_frame(do.call(rbind, corman.static[[4]]))

# create summary data frame
static.out <- corman2
static.out$average <- corman.static.average$GPP
static.out$diatoms <- corman.static.diatoms$GPP
static.out$greens <- corman.static.greens$GPP
static.out$cyanos <- corman.static.cyanos$GPP
static.out$model <- "static"

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
    y <- c("A1" = 100, "P" =in.grid[i, "TP_in"], "N" = in.grid[i, "TN_in"], 
           "QP1" = 0.015, "QN1" = 0.1)
    run <- ode(y, times, parms = params, func = dynamic.stoich.zmix)
    return(run[max(times),])
  })
})
(end <- Sys.time())
time.elapsed <- (end - start)/60/60
print(paste0("Time elapsed = ", time.elapsed, " hours!"))

# extract from list and convert to data-frame
corman.dynamic.average <- as_data_frame(do.call(rbind, corman.dynamic[[1]]))
corman.dynamic.diatoms <- as_data_frame(do.call(rbind, corman.dynamic[[2]]))
corman.dynamic.greens <- as_data_frame(do.call(rbind, corman.dynamic[[3]]))
corman.dynamic.cyanos <- as_data_frame(do.call(rbind, corman.dynamic[[4]]))

# create summary data frame
dynamic.out <- corman2
dynamic.out$average <- corman.dynamic.average$GPP
dynamic.out$diatoms <- corman.dynamic.diatoms$GPP
dynamic.out$greens <- corman.dynamic.greens$GPP
dynamic.out$cyanos <- corman.dynamic.cyanos$GPP
dynamic.out$model <- "dynamic"

################################################################################

# Plot base predictions
base.predictions.corman <- bind_rows(static.out, dynamic.out)
base.predictions.corman$model <- factor(base.predictions.corman$model, levels = c("static", "dynamic"))
# wide to long
base.predictions.corman <- base.predictions.corman %>%
  gather("species", "est_GPP", -Lake, -c(1:16), -model) %>%
  mutate(species = factor(species, levels = c("average", "diatoms", "greens", "cyanos")))

# create summary/RMSE across species and models
base.rmse <- base.predictions.corman %>%
  mutate(sqe = (est_GPP - GPP)^2) %>%
  group_by(model, species) %>%
  summarize(rmse = sqrt(mean(sqe))) %>%
  pivot_wider(id_cols = model, names_from = species, values_from = rmse)
# print  rmse
base.rmse
base.rmse[,-1] <- apply(base.rmse[,-1], 2, round, 2)
# calculate rsq values
rsq <- expand.grid(model = c("static", "dynamic"), 
                   species = c("average", "diatoms", "greens", "cyanos"), 
                   Lake = unique(corman2$Lake))
rsq$rsq <- NA
r <- lapply(1:nrow(rsq), function(i) {
  subset <- base.predictions.corman %>% filter(species == rsq[i, "species"]  & model == rsq[i, "model"] & Lake ==rsq[i, "Lake"])
  subset <- subset %>%
    mutate(log10GPP = log10(GPP), 
           log10estGPP = log10(est_GPP)) %>%
    mutate(log10GPP = ifelse(is.infinite(log10GPP), NA, log10GPP))
  r <- summary(lm(log10estGPP ~ log10GPP, data = subset))$adj.r.squared
  return(r)
})
rsq$rsq <- round(unlist(r), 2)
rsq %>% ggplot(aes(Lake, rsq, col = species)) + geom_point() + facet_wrap(model~.)
rsq.long <- rsq %>% pivot_wider(id_cols = c(Lake, model), names_from = species, values_from = rsq)

### Plot base predictions
(predicted.plt <- base.predictions.corman %>%
    #filter(Lake != "Feeagh" & Lake != "Acton" & Lake != "Langtjern" & Lake != "Trout") %>%
    ggplot() + 
    geom_smooth(aes(GPP, est_GPP), method = "lm", alpha = 0.3) + 
    geom_point(aes(GPP, est_GPP, pch = Lake), size = 2) + 
    ggh4x::facet_grid2(model ~ species,  scales = "free", independent = "y") + 
    scale_x_log10() + scale_y_log10() + 
    labs(x = "Measured GPP mg C L^-1 day^-1", 
         y = "Modelled GPP mg C L^-1 day^-1",
         shape = NULL) + 
    scale_shape_manual(
      values = c(
        "Acton" = 0,
        "EastLong" = 1,
        "Feeagh" = 2,
        "Harp" = 3,
        "Langtjern" = 4,
        "Lillinoah" = 5,
        "Lillsjoliden" = 6,
        "Mangstrettjarn" = 7,
        "Mendota" = 8,
        "Morris" = 9,
        "Struptjarn" = 10,
        "Trout" = 11,
        "Vortsjarv" = 12
      )
    ) + 
    theme(legend.position = "right")) 

################################################################################

# Perform grid optimization

### Static mdoel

#source("static_grid_optim.R") # takes almost 4hrs to run
load("corman_static_grid_optim.Rdata")

# select top 10 runs from static grid optim
# select top 2.5 percentile of runs based on Rsq
static.grid.optim <- arrange(static.grid.opt.metrics,RMSE )
#static.top10 <- static.grid.optim[1:10,]
static.top10 <- static.grid.optim[static.grid.optim$RMSE < quantile(static.grid.optim$RMSE, 0.025), ]

# re-run ten top scenarios to generate ensemble 
static.ensemble.in <- static.top10 %>%
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
static.ensemble.in$Lake <- rep(corman2$Lake, nrow(static.top10))
static.ensemble.in$Year <- rep(corman2$Year, nrow(static.top10))
static.ensemble.in$Month <- rep(corman2$Month, nrow(static.top10))
static.ensemble.in$GPP <- rep(corman2$GPP, nrow(static.top10))
static.ensemble.in$sd.GPP <- rep(corman2$GPP_sd, nrow(static.top10))
static.ensemble.in$TN_in <- rep(corman2$TN_in, nrow(static.top10))
static.ensemble.in$TP_in <- rep(corman2$TP_in, nrow(static.top10))
static.ensemble.in$HRT <- rep(corman2$HRT, nrow(static.top10))
static.ensemble.in$SA <- rep(corman2$SA, nrow(static.top10))
static.ensemble.in$DOC_mgL <- rep(corman2$DOC_mgL, nrow(static.top10))
static.ensemble.in$z <- rep(corman2$z, nrow(static.top10))
static.ensemble.in <- as.data.frame(static.ensemble.in)

# run models
static.ensemble <- lapply(1:nrow(static.ensemble.in), function(i){
  #print(i)
  # indexing
  static.algae["SA"] = static.ensemble.in[i, "SA"]
  static.algae["DOC"] = static.ensemble.in[i, "DOC_mgL"]
  static.algae["z"] = static.ensemble.in[i, "z"]
  static.algae["KP1"] = static.ensemble.in[i, "KP1"]
  static.algae["KN1"] = static.ensemble.in[i, "KN1"]
  static.algae["QP1"] = static.ensemble.in[i, "minQP1"]
  static.algae["QN1"] = static.ensemble.in[i, "minQN1"]
  static.algae["umax1"] = static.ensemble.in[i, "umax1"]
  static.algae["Pin"] = static.ensemble.in[i, "TP_in"]
  static.algae["Nin"] = static.ensemble.in[i, "TN_in"]
  static.algae["HRT"] = static.ensemble.in[i, "HRT"]
  # starting values
  y <- c("A1" = 100, "P" = static.ensemble.in[i, "TP_in"], "N" = static.ensemble.in[i, "TN_in"])
  run <- ode(y, times, parms = static.algae, func = static.stoich.zmix)
  return(run[max(times),])
})

# convert to df
static.ensemble <- do.call(rbind, static.ensemble)
static.ensemble <- as.data.frame(static.ensemble)

# rename columns
colnames(static.ensemble) <- paste0("est.", colnames(static.ensemble))
# add-in variables
static.ensemble <- bind_cols(static.ensemble, static.ensemble.in)

# calculate mean GPP and se per lake
static.ensemble.aggregate <- static.ensemble %>%
  group_by(Lake, Year, Month) %>%
  mutate(ave.est.GPP = mean(est.GPP), 
         sd.est.GPP = sd(est.GPP)) %>%
  ungroup() %>%
  mutate(model = "static", species = "optimal")


### dynamic mdoel

#source("dynamic_grid_optim.R") # takes almost 4hrs to run
load("corman_dynamic_grid_optim.Rdata")

# select top 10 runs from dynamic grid optim
dynamic.grid.optim <- arrange(dynamic.grid.opt.metrics,RMSE )
# select top 2.5 percentile of runs based on Rsq
dynamic.grid.optim <- arrange(dynamic.grid.optim,RMSE )
#dynamic.top10 <- dynamic.grid.optim[1:10,]
dynamic.top10 <- dynamic.grid.optim[dynamic.grid.optim$RMSE < quantile(dynamic.grid.optim$RMSE, 0.025), ]


# re-run ten top scenarios to generate ensemble 
dynamic.ensemble.in <- dynamic.top10 %>%
  slice(rep(1:n(), each = nrow(corman2)))

# add in corman2 data
dynamic.ensemble.in$Lake <- rep(corman2$Lake, nrow(dynamic.top10))
dynamic.ensemble.in$Year <- rep(corman2$Year, nrow(dynamic.top10))
dynamic.ensemble.in$Month <- rep(corman2$Month, nrow(dynamic.top10))
dynamic.ensemble.in$GPP <- rep(corman2$GPP, nrow(dynamic.top10))
dynamic.ensemble.in$sd.GPP <- rep(corman2$GPP_sd, nrow(dynamic.top10))
dynamic.ensemble.in$TN_in <- rep(corman2$TN_in, nrow(dynamic.top10))
dynamic.ensemble.in$TP_in <- rep(corman2$TP_in, nrow(dynamic.top10))
dynamic.ensemble.in$HRT <- rep(corman2$HRT, nrow(dynamic.top10))
dynamic.ensemble.in$SA <- rep(corman2$SA, nrow(dynamic.top10))
dynamic.ensemble.in$DOC_mgL <- rep(corman2$DOC_mgL, nrow(dynamic.top10))
dynamic.ensemble.in$z <- rep(corman2$z, nrow(dynamic.top10))
dynamic.ensemble.in <- as.data.frame(dynamic.ensemble.in)

# run models
dynamic.ensemble <- lapply(1:nrow(dynamic.ensemble.in), function(i){
  #print(i)
  # indexing
  dynamic.algae["SA"] = dynamic.ensemble.in[i, "SA"]
  dynamic.algae["DOC"] = dynamic.ensemble.in[i, "DOC_mgL"]
  dynamic.algae["z"] = dynamic.ensemble.in[i, "z"]
  dynamic.algae["KP1"] = dynamic.ensemble.in[i, "KP1"]
  dynamic.algae["KN1"] = dynamic.ensemble.in[i, "KN1"]
  dynamic.algae["QP1"] = dynamic.ensemble.in[i, "minQP1"]
  dynamic.algae["QN1"] = dynamic.ensemble.in[i, "minQN1"]
  dynamic.algae["upP1"] = dynamic.ensemble.in[i, "upP1"]
  dynamic.algae["upN1"] = dynamic.ensemble.in[i, "upN1"]
  dynamic.algae["umax1"] = dynamic.ensemble.in[i, "umax1"]
  dynamic.algae["Pin"] = dynamic.ensemble.in[i, "TP_in"]
  dynamic.algae["Nin"] = dynamic.ensemble.in[i, "TN_in"]
  dynamic.algae["HRT"] = dynamic.ensemble.in[i, "HRT"]
  # starting values
  y <- c("A1" = 100, "P" = dynamic.ensemble.in[i, "TP_in"], "N" = dynamic.ensemble.in[i, "TN_in"],
         "QP1" = 0.015, "QN1" = 0.1)
  run <- ode(y, times, parms = dynamic.algae, func = dynamic.stoich.zmix)
  return(run[max(times),])
})

# convert to df
dynamic.ensemble <- do.call(rbind, dynamic.ensemble)
dynamic.ensemble <- as.data.frame(dynamic.ensemble)

# rename columns
colnames(dynamic.ensemble) <- paste0("est.", colnames(dynamic.ensemble))
# add-in variables
dynamic.ensemble <- bind_cols(dynamic.ensemble, dynamic.ensemble.in)

# calculate mean GPP and se per lake
dynamic.ensemble.aggregate <- dynamic.ensemble %>%
  group_by(Lake, Year, Month) %>%
  mutate(ave.est.GPP = mean(est.GPP), 
         sd.est.GPP = sd(est.GPP)) %>%
  ungroup() %>%
  mutate(model = "dynamic", species = "optimal")


### bind ensembles together
ensembles.all <- bind_rows(static.ensemble.aggregate, dynamic.ensemble.aggregate)

################################################################################
# plotting

# combine data set with baseline predictions for pretty plotting
all.predictions <- bind_rows(
  base.predictions.corman,
  ensembles.all %>%
    # only selects needed to match, discarding a lot of information here
    select(Lake, Month, Year, TP_in, TN_in, HRT, 
           SA, DOC_mgL, z, species, model, 
           GPP, est_GPP = ave.est.GPP, sd.est.GPP)) %>%
  mutate(model = factor(model, levels = c("static", "dynamic")), 
         species = factor(species, levels = c("average", "diatoms", "greens", "cyanos", "optimal")))

### Plot base predictions
(predicted.plt <- all.predictions %>%
    ggplot() + 
    geom_smooth(aes(GPP, est_GPP), method = "lm", alpha = 0.3) + 
    geom_point(aes(GPP, est_GPP, pch = Lake), size = 2) + 
    ggh4x::facet_grid2(model ~ species,  scales = "free", independent = "y") + 
    scale_x_log10() + scale_y_log10() + 
    labs(x = "Measured GPP mg C L^-1 day^-1", 
         y = "Modelled GPP mg C L^-1 day^-1",
         shape = NULL) + 
    scale_shape_manual(
      values = c(
        "Acton" = 0,
        "EastLong" = 1,
        "Feeagh" = 2,
        "Harp" = 3,
        "Langtjern" = 4,
        "Lillinoah" = 5,
        "Lillsjoliden" = 6,
        "Mangstrettjarn" = 7,
        "Mendota" = 8,
        "Morris" = 9,
        "Struptjarn" = 10,
        "Trout" = 11,
        "Vortsjarv" = 12
      )
    ) + 
    theme(legend.position = "right")) 

# generate rsq values
# calculate rsq values
rsq <- expand.grid(model = c("static", "dynamic"), 
                   species = c("average", "diatoms", "greens", "cyanos", "optimal"))
rsq$rsq <- NA
r <- lapply(1:nrow(rsq), function(i) {
  subset <- all.predictions %>% filter(species == rsq[i, "species"]  & model == rsq[i, "model"])
  subset <- subset %>%
    mutate(log10GPP = log10(GPP), 
           log10estGPP = log10(est_GPP)) %>%
    mutate(log10GPP = ifelse(is.infinite(log10GPP), NA, log10GPP))
  r <- summary(lm(log10estGPP ~ log10GPP, data = subset))$adj.r.squared
  return(r)
})
rsq$rsq <- round(unlist(r), 2)


# plot range of trait valyes
# lack of range in miNQP and minQN means best runs always have min possible cell quotas
# further evidence that minQN:miNQP drives GPP estimates
best.traits <- bind_rows(
  static.top10 %>% mutate(model = "static"),
  dynamic.top10 %>% mutate(model = "dynamic")
) %>% mutate(model = factor(model, levels = c("static", "dynamic")))

best.traits %>% 
  gather("Trait", "Trait.value", 
         -model, -RMSE, -MAE, -NSE) %>%
  ggplot() + 
  geom_point(aes(model, Trait.value)) + 
  facet_wrap(.~Trait, scales = "free_y")
