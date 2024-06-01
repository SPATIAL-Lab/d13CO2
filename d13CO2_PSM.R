model {
  
# Likelihood block 
####################################################################################################
  
  # Benthic forams
  for (i in 1:length(d13Cbf.data)){
    d13Cbf.data[i] ~ dnorm(d13Cbf[ai.d13Cbf[i], si.d13Cbf[i]], d13Cbf.p[i])
    d13Cbf.p[i] = 1/d13Cbf.sd[i]^2
  }

  # Planktic forams
  for (i in 1:length(d13Cpf.data)){
    d13Cpf.data[i] ~ dnorm(d13Cpf[ai.d13Cpf[i], si.d13Cpf[i]], d13Cpf.p[i])
    d13Cpf.p[i] = 1/d13Cpf.sd[i]^2
  }
  
  
# Constants
####################################################################################################  
  # Equation 8 of Tipple et al. (2010)
  eps.dic_cc <- -1   
  
  # Benthic disequilibrium effects; totals to 2.8 (mean) from Tipple et al. (2010)
  asd ~ dnorm(1, 1/0.2^2)T(0,)  # air sea disequilibrium 
  bpump ~ dnorm(1.2, 1/0.4^2)T(0,) # biological pump
  remin ~ dnorm(0.6, 1/0.3^2)T(0,) # remineralization and oxidation 
  A = 2.8 # asd + bpump + remin
  
  
# Proxy system model 
####################################################################################################
  for (i in 1:n.steps){
    for (j in 1:n.sites){
      # Calculate various carbonate archive d13C values from d13CO2 and temperature
      
      # Equations 4, 5, and 3 (respectively) of Tipple et al. (2010)
      eps.bicarb_co2[i,j] <- -0.1141*tempC[i,j] + 10.78                   
      eps.ci_co2[i,j] <- 0.0049*tempC[i,j] - 1.31
      eps.dic_co2[i,j] <- (0.91*eps.bicarb_co2[i,j]) + (0.08*eps.ci_co2[i,j])
      
      # Benthic forams; equation 7 of Tipple et al. (2010)
      d13Cbf[i,j] <- ((d13CO2[i]+1000)*((eps.dic_co2[i,j]/1000)+1)) - eps.dic_cc - A - 1000
      
      # Planktic forams; equation 9 of Tipple et al. (2010)
      d13Cpf[i,j] <- ((d13CO2[i]+1000)*((eps.dic_co2[i,j]/1000)+1)) - eps.dic_cc - 1000
    }
  }
  
  
# Time evolution model  
####################################################################################################
  # d13C of atmospheric CO2
  d13CO2[1] ~ dnorm(d13CO2.m, 1/d13CO2.sd^2)T(-10,-6)    
  d13CO2.phi ~ dbeta(5,2) 
  d13CO2.eps[1] = 0 
  d13CO2.tau ~ dgamma(100, 10)
  
  # Global mean surface temperature 
  GMST[1] ~ dnorm(GMST.m[1], 1/GMST.sd[1]^2)T(GMST.m[1]-GMST.sd[1]*2,GMST.m[1]+GMST.sd[1]*2) 
  
  # Temp in C
  for (j in 1:n.sites){
    toff[1,j] ~ dnorm(toff.m[1,j], toff.sd[1,j])
    tempC[1,j] = GMST[1] + toff[1,j]
  }
  
  for (i in 2:n.steps){
    # d13C of atmospheric CO2
    d13CO2.pc[i] <- d13CO2.tau*((1-d13CO2.phi^2)/(1-d13CO2.phi^(2*dt[i-1])))
    d13CO2.eps[i] ~ dnorm(d13CO2.eps[i-1]*(d13CO2.phi^dt[i-1]), d13CO2.pc[i])T(-2, 2)
    d13CO2[i] <- d13CO2[1] + d13CO2.eps[i]
    
    # Global mean surface temperature 
    GMST[i] ~ dnorm(GMST.m[i], 1/GMST.sd[i]^2)T(GMST.m[i]-GMST.sd[i]*2,GMST.m[i]+GMST.sd[i]*2)  
    
    for (j in 1:n.sites){
      # Temp in C
      toff[i,j] ~ dnorm(toff.m[i,j], toff.sd[i,j])
      tempC[i,j] = GMST[i] + toff[i,j]
    }
  }
}






