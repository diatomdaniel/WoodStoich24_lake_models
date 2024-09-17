######################################################################

# R script that contains mechanistic models

######################################################################

# light attenuation model
lightAtten<-function(z,I0,kD){
  Iz=I0*exp(-kD*z)
  return(Iz)
}

### Stoichiometry model using Droop equation
# Huisman and Weissing 1995, Kelly et al 2014, Hall et al 2007
# build model
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
  zmix = params["zmix"]
  Pin = params["Pin"]
  Nin = params["Nin"]
  
  # light parameters
  I0 <- params["I0"]
  kBg = params["kBg"]
  kA = params["kA"]
  
  # algae physiology parameters
  umax1 = params["umax1"]
  lA = params["lA"]
  v = params["v"]
  KLight = params["KLight"]
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
  
  # In/output
  Qin=SA*1e6*zmix/365	# m^3 day^-1
  # kd = function of algae, background light absorption, and light absorption due to ash
  kD <- kA * A1 + kBg 
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * zmix
  # light attenuation
  Izmix=lightAtten(z=zmix,I0=I0,kD=kD)
  
  # biomass specific growth for entire mixed layer
  prod1 =(umax1/(kD*zmix))*log((KLight + I0)/(KLight+Izmix))* (umax1 * min(1 - minQN1/QN1, 1 - minQP1/QP1 ))	# d-1
  
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
  dY=c(dA1dt=dA1.dt, dPdt=dP.dt, dNdt = dN.dt, dQP1dt = dQP1.dt, dQ1Ndt = dQN1.dt)
  return(list(dY))
  
} 

