######################################################################

# R script that contains mechanistic models

# Models are based on from Huissman and Weissing 1995 Am. Nat., Kelly et al 2018 Ecosyst, Jäger and Diehl 2014 Ecology, Hall et al 2007 Ecology, and Olson et al 2022.
#Daniel Gschwentner January 2024

# run this file as source to call models

# light attenuation model
lightAtten<-function(z,I0,kD){
  Iz=I0*exp(-kD*z)
  return(Iz)
}

### nitrogen-phosphorus model
# Huisman and Weissing 1995, Kelly et al 2014, Jäger and Diehl 2014
# build model
mich.trip <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  A2 <- y["A2"]
  A3 <- y["A3"]
  P <- y["P"]
  N <- y["N"]
  
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
  
  # species 1
  KP1 = params["KP1"]
  QP1 = params["QP1"]
  KN1 = params["KN1"]
  QN1 = params["QN1"]
  # species 2
  KP2 = params["KP2"]
  QP2 = params["QP2"]
  KN2 = params["KN2"]
  QN2 = params["QN2"]
  # species 3
  KP3 = params["KP3"]
  QP3 = params["QP3"]
  KN3 = params["KN3"]
  QN3 = params["QN3"]
  
  # In/output
  Qin=SA*1e6*zmix/365	# m^3 day^-1
  # kd = function of algae, background light absorption, and light absorption due to ash
  kD <- kA * (A1 +A2 + A3) + kBg
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * zmix
  # light attenuation
  Izmix=lightAtten(z=zmix,I0=I0,kD=kD)
  
  # biomass specific growth for entire mixed layer
  # species 1
  prod1=(umax1/(kD*zmix))*log((KLight1 + I0)/(KLight1+Izmix))*(P/(P+KP1))*(N/(N + KN1))	# d-1
  #prod1=(P/(P+KP1))*(N/(N + KN1))	# d-1
  
  # species 2
  prod2=(umax2/(kD*zmix))*log((KLight2 + I0)/(KLight2+Izmix))*(P/(P+KP2))*(N/(N + KN2))	# d-1
  #prod2=(P/(P+KP2))*(N/(N + KN2))
  # species 3
  prod3=(umax3/(kD*zmix))*log((KLight3 + I0)/(KLight3+Izmix))*(P/(P+KP3))*(N/(N + KN3))	# d-
  #prod3=(P/(P+KP3))*(N/(N + KN3))
  
  # model biomass
  # species 1
  dA1.dt=A1*prod1-lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# mg C m-3
  # species 2
  dA2.dt=A2*prod2-lA*A2-v/zmix*A2-Qin/(zmix*SA*1e6)*A2	# mg C m-3
  # species 3
  dA3.dt=A3*prod3-lA*A3-v/zmix*A3-Qin/(zmix*SA*1e6)*A3	# mg C m-3
  
  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P)+QP1*lA*A1-QP1*A1*prod1 + QP2*lA*A2-QP2*A2*prod2 +QP3*lA*A3-QP3*A3*prod3 #mg P m-3 (in epi);
  
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N) + QN1*lA*A1-QN1*A1*prod1 + QN2*lA*A2-QN2*A2*prod2 + QN3*lA*A3-QN3*A3*prod3  # mg N m-3 (in epi);
  
  # # indicators of limitation
  # Plim <- 1 - (P/(P + KP)) # unitless
  # Nlim <- 1 - (N/(N + KN))
  # Llim <- 1 - (1/(kD * zmix)) * log((KLight + I0)/(KLight + Izmix)) # unitless
  # NP_moles <- (N/14.007)/(P/30.974)
  
  # return objects
  dY=c(dA1dt=dA1.dt,dA2dt=dA2.dt,dA3dt=dA3.dt, dPdt=dP.dt, dNdt = dN.dt)
  #lim=c(Plim = Plim, Nlim = Nlim, Llim = Llim, NP_moles = NP_moles, growth.rate = prod)
  return(list(dY))
  
}
