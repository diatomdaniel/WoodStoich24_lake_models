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
times <- 1:2000 # for troubleshooting, initial runs

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

# merge and unit conversions, GPP standardization
base.predictions.corman <- bind_rows(static.out, dynamic.out)
base.predictions.corman$model <- factor(base.predictions.corman$model, levels = c("static", "dynamic"))
# wide to long
base.predictions.corman <- base.predictions.corman %>%
  gather("species", "est_GPP", -Lake, -c(1:16), -model) %>%
  mutate(species = factor(species, levels = c("average", "diatoms", "greens", "cyanos"))) %>%
  group_by(model, Lake) %>%
  mutate(zscore_gpp = (est_GPP - mean(est_GPP))/sd(est_GPP)) %>%
  ungroup() %>%
  mutate(zscore_obs_gpp = (GPP - mean(GPP))/sd(GPP))

################################################################################

### Decided to visualize density distributions of z-scores for models and measured GPP; better for pattern comparison? 
(gpp.density <- base.predictions.corman %>%
  ggplot() + 
  geom_density(aes(log(GPP + 0.1)),  col = "black", alpha = 0.3) + 
  geom_density(aes(log(est_GPP + 0.1)), fill = "grey20", col = "black", alpha = 0.3) + 
  ggh4x::facet_grid2(model~species, independent = "y", scales = "free_y") + 
  labs(x = expression("log(GPP mg C L"^-1 ~ "day"^-1 ~"+ 0.1)"), y = "Density"))

save_plot("figures/gpp_validation_densities.png", gpp.density, base_height = 6, base_width = 9)

# add regression coefs (manually)
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "static" & base.predictions.corman$species == "average", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "static" & base.predictions.corman$species == "diatoms", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "static" & base.predictions.corman$species == "greens", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "static" & base.predictions.corman$species == "cyanos", ]))

summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "dynamic" & base.predictions.corman$species == "average", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "dynamic" & base.predictions.corman$species == "diatoms", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "dynamic" & base.predictions.corman$species == "greens", ]))
summary(lm(log(GPP + 0.1) ~ log(est_GPP + 0.1), data = base.predictions.corman[base.predictions.corman$model == "dynamic" & base.predictions.corman$species == "cyanos", ]))


################################################################################

# calculate error metrics and relate to inflows and lake-specific ,factors
err.metrics <- base.predictions.corman %>%
  group_by(Lake, species, model, Month) %>%
  summarise("N:P" = (TN_in/14.007)/(TP_in/30.974), 
            "TP umoles/L" = (TP_in/30.974), "TN umoles/L" = (TN_in),
            "Lake area km2" = (SA), "DOC mg/L" = (DOC_mgL),
            "Epi. HRT days" = (HRT)/365,
            MAE = sum(abs(est_GPP - GPP))/n())

err.metrics %>%
  gather("Key", "Value", -model, -Lake, -species, -MAE) %>%
  ggplot(aes(log(Value), log(MAE), pch = species)) + 
  geom_point() + 
  ggh4x::facet_grid2(model~Key, independent = "all", scales = "free")


summary(lm(log(MAE) ~ log(`TP umoles/L`), data = err.metrics[err.metrics$model == "static",]))
summary(lm(log(MAE) ~ log(`TP umoles/L`), data = err.metrics[err.metrics$model == "dynamic",]))

summary(lm(log(MAE) ~ log(`TN umoles/L`), data = err.metrics[err.metrics$model == "static",]))
summary(lm(log(MAE) ~ log(`TN umoles/L`), data = err.metrics[err.metrics$model == "dynamic",]))


summary(lm(log(MAE) ~ log(`N:P`), data = err.metrics[err.metrics$model == "static",]))
summary(lm(log(MAE) ~ log(`N:P`), data = err.metrics[err.metrics$model == "dynamic",]))


################################################################################

# plot predictions following Kelly et al 2018 F5 panels C and D

(obs.plt <- corman2 %>%
  ggplot(aes(log(TP_in/30.974), log(GPP + 0.1), col = (TN_in/14.007)/(TP_in/30.974))) + 
  geom_point(size = 2) + 
  #scale_x_log10() +  #scale_y_log10() + 
  scale_color_viridis_c() + 
  labs(x = expression("log(load TP umoles L"^-1 ~ ")"), 
       y = expression("log(GPP + 0.1) mg C L"^-1 ~ " day"^-1),
       col = "Load N:P (molar)"))

(mod.plt <- base.predictions.corman %>%
    ggplot(aes(log(TP_in/30.974), log(est_GPP + 0.1),col = (TN_in/14.007)/(TP_in/30.974))) + 
    geom_point(size = 2) + 
    ggh4x::facet_grid2(model~species, scales = "free", independent = "y") + 
    #() +  #scale_y_log10() +
    scale_color_viridis_c() + 
    labs(x = expression("log(load TP umoles L"^-1 ~ ")"), 
         y = expression("log(modeled GPP + 0.1) mg C L"^-1 ~ " day"^-1),
         col = "Load N:P (molar)") + 
    guides(shape = "none"))

### arrange
subplot.b <- ggarrange(plotlist = list(obs.plt, mod.plt), 
                       nrow = 1, ncol = 2, align = "hv",
                       widths = c(1/3, 2/3),
                       labels = c("b", "c"), 
                       common.legend = T, legend = "bottom")
subplot.b


### Plot base predictions
(predicted.plt <- base.predictions.corman %>%
    #filter(Lake != "Feeagh" & Lake != "Acton" & Lake != "Langtjern" & Lake != "Trout") %>%
    ggplot() + 
    #geom_smooth(aes(GPP, est_GPP), method = "lm", alpha = 0.3) + 
    geom_point(aes(log(est_GPP + 0.1), log(GPP + 0.1), pch = Lake), size = 2) + 
    ggh4x::facet_grid2(model ~ species,  scales = "free", independent = "x") + 
    #geom_abline(slope = 1, intercept = 0) + 
    #scale_x_log10() + scale_y_log10() + 
    labs(y = expression("Obs. GPP (mg C L"^-1 ~ " day"^-1~")"), 
         x = expression("log(mod. GPP + 0.1) mg C L"^-1 ~ " day"^-1),
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

subplot.a <- predicted.plt

################################################################################

### arrange everything into one figure
corman.validation.plt <- ggarrange(plotlist = list(subplot.a, subplot.b),
                                   nrow = 2, ncol = 1, align = "hv", 
                                   labels = c("a", ""))
corman.validation.plt
#save_plot("figures/gpp_validation.png", corman.validation.plt, base_height = 8, base_width = 14)
