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
droop.trip <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  A2 <- y["A2"]
  A3 <- y["A3"]
  P <- y["P"]
  N <- y["N"]
  QP1 <- y["QP1"]
  QP2 <- y["QP2"]
  QP3 <- y["QP3"]
  QN1 <- y["QN1"]
  QN2 <- y["QN2"]
  QN3 <- y["QN3"]
  
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
  umax2 = params["umax2"]
  umax3 = params["umax3"]
  
  lA = params["lA"]
  v = params["v"]
  KLight1 = params["KLight1"]
  KLight2 = params["KLight2"]
  KLight3 = params["KLight3"]
  # half-sat. constants for P
  KP1 = params["KP1"]
  KP2 = params["KP2"]
  KP3 = params["KP3"]
  # min. cell quota for P
  minQP1 = params["minQP1"]
  minQP2 = params["minQP2"]
  minQP3 = params["minQP3"]
  # uptake rate P
  upP1 = params["upP1"]
  upP2 = params["upP2"]
  upP3 = params["upP3"]
  # half-sat. constants for N
  KN1 = params["KN1"]
  KN2 = params["KN2"]
  KN3 = params["KN3"]
  # min cell quota for N
  minQN1 = params["minQN1"]
  minQN2 = params["minQN2"]
  minQN3 = params["minQN3"]
  # uptake rate N
  upN1 = params["upN1"]
  upN2 = params["upN2"]
  upN3 = params["upN3"]
  
  # In/output
  Qin=SA*1e6*zmix/365	# m^3 day^-1
  # kd = function of algae, background light absorption, and light absorption due to ash
  kD <- kA * (A1 +A2 + A3) + kBg
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * zmix
  # light attenuation
  Izmix=lightAtten(z=zmix,I0=I0,kD=kD)
  
  # biomass specific growth for entire mixed layer
  prod1 =(umax1/(kD*zmix))*log((KLight1 + I0)/(KLight1+Izmix))* (umax1 * min(1 - minQN1/QN1, 1 - minQP1/QP1 ))	# d-1
  prod2 =(umax2/(kD*zmix))*log((KLight2 + I0)/(KLight2+Izmix))* (umax2 * min(1 - minQN2/QN2, 1 - minQP2/QP2 ))	# d-1
  prod3 =(umax3/(kD*zmix))*log((KLight3 + I0)/(KLight3+Izmix))* (umax3 * min(1 - minQN3/QN3, 1 - minQP3/QP3 ))	# d-1
  
  # model biomass
  dA1.dt=A1*prod1 - lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# mg C m-3
  dA2.dt=A2*prod2 - lA*A2-v/zmix*A2-Qin/(zmix*SA*1e6)*A2	# mg C m-3
  dA3.dt=A3*prod3 - lA*A3-v/zmix*A3-Qin/(zmix*SA*1e6)*A3	# mg C m-3
  
  # cell quota P  
  dQP1.dt = upP1 * (P/(KP1 + P)) - prod1  * QP1
  dQP2.dt = upP2 * (P/(KP2 + P)) - prod2  * QP2
  dQP3.dt = upP3 * (P/(KP3 + P)) - prod3  * QP3
  
  # cell quota N  
  dQN1.dt = upN1 * (N/(KN1 + N)) - prod1  * QN1
  dQN2.dt = upN2 * (N/(KN2 + N)) - prod2  * QN2
  dQN3.dt = upN3 * (N/(KN3 + N)) - prod3  * QN3

  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P) + A1 * (-upP1 * (P/(KP1 + P))  + lA * QP1) + A2 * (-upP2 * (P/(KP2 + P))  + lA * QP2) + A3 * (-upP3 * (P/(KP3 + P))  + lA * QP3)  # mg P m-3 (in epi);
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N) + A1 * (-upN1 * (N/(KN1 + N)) +  lA * QN1) + A2 * (-upN2 * (N/(KN2 + N)) +  lA * QN2) + A3 * (-upN3 * (N/(KN3 + N)) +  lA * QN3)  # mg N m-3 (in epi);

  # return objects 
  dY=c(dA1dt=dA1.dt, dA2dt=dA2.dt, dA3dt=dA3.dt, 
       dPdt=dP.dt, dNdt = dN.dt, 
       dQP1dt = dQP1.dt, dQP2dt = dQP2.dt, dQP3dt = dQP3.dt,
       dQ1Ndt = dQN1.dt, dQ2Ndt = dQN2.dt, dQ3Ndt = dQN3.dt)
  return(list(dY))
  
} 

