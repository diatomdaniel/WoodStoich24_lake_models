

### Code to store algae parameter vectors

# traits are from Olson et al 2024 in prep.
# trait values aggregated from Edwards et al 2015 https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3562857?file=5635515

# create table with trait values
traits <- tibble::tibble(
  Trait = c("K_N_mgL", "K_P_mgL", "Qmax_mg_N_mg_C", "Qmax_mg_P_mg_C", "Qmin_mg_N_mg_C", 
            "Qmin_mg_P_mg_C", "Vmax_mg_N_mg_C", "Vmax_mg_P_mg_C", "umax_N", "umax_P", "K_light", "umax"),
  diatoms = c(0.064, 0.005, 0.201, 0.146, 0.155, 0.02, 4.49, 0.608, 0.86, 0.61, 167.8, 0.455),
  greens = c(0.036, 0.028, 0.068, 0.015, 0.025, 0.001, 0.416, 0.243, 1.32, 0.63, 165.25, 0.875),
  cyanos = c(0.033, 0.026, 0.023, 0.008, 0.01, 0.001, 0.034, 0.113, 0.86, 0.655, 88.6, 0.64),
  average = c(0.05, 0.0165, 0.1345, 0.0805, 0.09, 0.0105, 2.453, 0.4255, 1.09, 0.62, 166.525, 0.665)
)

# average algae static model

static.algae <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  # varies by simulation
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "average"],
  lA=0.1,			# mortality rate day-1
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "average"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  QP1 = traits[6, "average"], # algae cell P quota in mg P mg^-1 C^-1
  KN1 = traits[1, "average"] * 1000, # nitrogen half sat constant in mg N m^-3 
  QN1 = traits[5, "average"] # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.algae) <- c("SA", "zmix", "Pin", "Nin", "HRT", "umax1", "lA", "v", 
                         "KP1",  "QP1", "KN1",  "QN1")
names(static.algae)
static.algae <- unlist(static.algae)


# diatoms static model

static.diatoms <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  # varies by simulation
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "diatoms"],
  lA=0.1,			# mortality rate day-1
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "diatoms"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  QP1 = traits[6, "diatoms"], # algae cell P quota in mg P mg^-1 C^-1
  KN1 = traits[1, "diatoms"] * 1000, # nitrogen half sat constant in mg N m^-3 
  QN1 = traits[5, "diatoms"] # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.diatoms) <- c("SA", "zmix", "Pin", "Nin", "HRT", "umax1", "lA", "v", 
                           "KP1",  "QP1", "KN1",  "QN1")
names(static.diatoms)
static.diatoms <- unlist(static.diatoms)


# green algae static model

static.greens <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  # varies by simulation
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days,
  
  # algae physiology parameters
  umax1 = traits[12, "greens"],
  lA=0.1,			# mortality rate day-1
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "greens"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  QP1 = traits[6, "greens"], # algae cell P quota in mg P mg^-1 C^-1
  KN1 = traits[1, "greens"] * 1000, # nitrogen half sat constant in mg N m^-3 
  QN1 = traits[5, "greens"] # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.greens) <- c("SA", "zmix", "Pin", "Nin", "HRT", "umax1", "lA", "v", 
                          "KP1",  "QP1", "KN1",  "QN1")
names(static.greens)
static.greens <- unlist(static.greens)


# cyanobacteria static model

static.cyanos <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  # varies by simulation
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "cyanos"],
  lA=0.1,			# mortality rate day-1
  v=0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "cyanos"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  QP1 = traits[6, "cyanos"], # algae cell P quota in mg P mg^-1 C^-1
  KN1 = traits[1, "cyanos"] * 1000, # nitrogen half sat constant in mg N m^-3 
  QN1 = traits[5, "cyanos"] # algae cell N quota in mg N mg^-1 C^-1 
)
names(static.cyanos) <- c("SA", "zmix", "Pin", "Nin", "HRT", "umax1", "lA", "v",
                          "KP1",  "QP1", "KN1",  "QN1")
names(static.cyanos)
static.cyanos <- unlist(static.cyanos)


# average algae dynamic model

dynamic.algae <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "average"],
  lA=0.1,			# mortality rate day-1
  v= 0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "average"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  minQP1 = traits[6, "average"], # algae cell P quota in mg P mg^-1 C^-1 f
  upP1 = traits[8, "average"], # max uptake rate P per day in mg P mg C^-1 day^-1 
  KN1 = traits[1, "average"] * 1000, # nitrogen half sat constant in mg N m^-3 
  minQN1 = traits[5, "average"], # algae cell N quota in mg N mg^-1 C^-1 f
  upN1 = traits[7, "average"] # max uptake rate N per day in mg N mg C^-1 day^-1 
)

names(dynamic.algae) <- c("SA", "zmix", "Pin", "Nin", "HRT",
                          "umax1", "lA", "v",
                          "KP1", "minQP1", "upP1", "KN1",
                          "minQN1", "upN1")
names(dynamic.algae)
dynamic.algae <- unlist(dynamic.algae)


# diatoms dynamic model

dynamic.diatoms <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "diatoms"],
  lA=0.1,			# mortality rate day-1
  v= 0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "diatoms"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  minQP1 = traits[6, "diatoms"], # algae cell P quota in mg P mg^-1 C^-1 f
  upP1 = traits[8, "diatoms"], # max uptake rate P per day in mg P mg C^-1 day^-1 
  KN1 = traits[1, "diatoms"] * 1000, # nitrogen half sat constant in mg N m^-3 
  minQN1 = traits[5, "diatoms"], # algae cell N quota in mg N mg^-1 C^-1 f
  upN1 = traits[7, "diatoms"] # max uptake rate N per day in mg N mg C^-1 day^-1 
)

names(dynamic.diatoms) <- c("SA", "zmix", "Pin", "Nin", "HRT",
                            "umax1", "lA", "v","KP1", "minQP1", "upP1", "KN1",
                            "minQN1", "upN1")
names(dynamic.diatoms)
dynamic.diatoms <- unlist(dynamic.diatoms)


# greens dynamic model

dynamic.greens <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "greens"],
  lA=0.1,			# mortality rate day-1
  v= 0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "greens"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  minQP1 = traits[6, "greens"], # algae cell P quota in mg P mg^-1 C^-1 f
  upP1 = traits[8, "greens"], # max uptake rate P per day in mg P mg C^-1 day^-1 
  KN1 = traits[1, "greens"] * 1000, # nitrogen half sat constant in mg N m^-3 
  minQN1 = traits[5, "greens"], # algae cell N quota in mg N mg^-1 C^-1 f
  upN1 = traits[7, "greens"] # max uptake rate N per day in mg N mg C^-1 day^-1 
)

names(dynamic.greens) <- c("SA", "zmix", "Pin", "Nin", "HRT",
                           "umax1", "lA", "v", "KP1", "minQP1", "upP1", "KN1",
                           "minQN1", "upN1")
names(dynamic.greens)
dynamic.greens <- unlist(dynamic.greens)


# cyanos dynamic model

dynamic.cyanos <- c(
  # lake parameters
  SA= 1,		# lake surface area in km2
  zmix = 2, # lake mixing depth in m
  Pin = NA, # P inflow concentration in mg P m^-3
  Nin = NA, # N inflow concentration in mg N m^-3
  HRT = 365, # water residence time in days
  
  # algae physiology parameters
  umax1 = traits[12, "cyanos"],
  lA=0.1,			# mortality rate day-1
  v= 0.1,			# m d-1; sinking loss of algae
  KP1 = traits[2, "cyanos"] * 1000, # phosphorus half sat constant in mg P m^-3 f
  minQP1 = traits[6, "cyanos"], # algae cell P quota in mg P mg^-1 C^-1 f
  upP1 = traits[8, "cyanos"], # max uptake rate P per day in mg P mg C^-1 day^-1 
  KN1 = traits[1, "cyanos"] * 1000, # nitrogen half sat constant in mg N m^-3 
  minQN1 = traits[5, "cyanos"], # algae cell N quota in mg N mg^-1 C^-1 f
  upN1 = traits[7, "cyanos"] # max uptake rate N per day in mg N mg C^-1 day^-1 
)

names(dynamic.cyanos) <- c("SA", "zmix", "Pin", "Nin", "HRT",
                           "umax1", "lA", "v", "KP1", 
                           "minQP1", "upP1", "KN1", "minQN1", "upN1")
names(dynamic.cyanos)
dynamic.cyanos <- unlist(dynamic.cyanos)

