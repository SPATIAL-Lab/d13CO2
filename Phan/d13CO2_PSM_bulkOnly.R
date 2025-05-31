model {
  
  # Likelihood block 
  ####################################################################################################
  
  # Bulk carbonate
  d13Cbulk.analyt.sd <- 0.05
  d13Cbulk.sd <- sqrt(d13Cbulk.analyt.sd^2 + bulk.alteration^2)
  d13Cbulk.p <- 1/(d13Cbulk.sd)^2
  for (i in 1:n.d13Cbulk){
    d13Cbulk.data[i] ~ dnorm(d13Cbulk[ri.d13Cbulk[i]], d13Cbulk.p)
  }
  # 
  # # Benthic forams
  # d13Cbf.analyt.sd <- 0.05
  # d13Cbf.sd <- sqrt(d13Cbf.analyt.sd^2 + bf.alteration^2)
  # d13Cbf.p <- 1/d13Cbf.sd^2
  # for (i in 1:n.d13Cbf){
  #   d13Cbf.data[i] ~ dnorm(d13Cbf[ri.d13Cbf[i]], d13Cbf.p)
  # }
  # 
  # # Brachiopods
  # d13Cbrach.analyt.sd <- 0.05
  # d13Cbrach.sd <- sqrt(d13Cbrach.analyt.sd^2 + brach.alteration^2)
  # d13Cbrach.p <- 1/0.05^2
  # for (i in 1:n.d13Cbrach){
  #   d13Cbrach.data[i] ~ dnorm(d13Cbrach[ri.d13Cbrach[i]], d13Cbrach.p)
  # }
  
  
  # Constants
  ####################################################################################################  
  # Equation 8 of Tipple et al. (2010)
  eps.dic_cc <- -1   
  
  # Benthic disequilibrium effects; totals to 2.8 (mean) from Tipple et al. (2010)
  asd ~ dnorm(1, 1/0.2^2)T(0,)  # air sea disequilibrium 
  bpump ~ dnorm(1.2, 1/0.4^2)T(0,) # biological pump
  remin ~ dnorm(0.6, 1/0.3^2)T(0,) # remineralization and oxidation 
  A = asd + bpump + remin
  
  # Bulk carbonate uncertainty terms from Barral et al. (2017)
  coeff1 ~ dnorm(0.0144, 1/0.01^2)
  coeff2 ~ dnorm(0.107, 1/0.002^2)
  eps.int ~ dnorm(10.53, 1/0.05^2)
  f_co2 ~ dnorm(0.12, 1/0.04^2)T(0.05,0.20)
  
  
  # Proxy system model 
  ####################################################################################################
  for (i in 1:length(ai.flat)){
   
      # Calculate various carbonate archive d13C values from d13CO2 and temperature
      # Equations 4, 5, and 3 (respectively) of Tipple et al. (2010)
      eps.bicarb_co2[i] <- -0.1141*tempC[i] + 10.78                   
      eps.ci_co2[i] <- 0.0049*tempC[i] - 1.31
      eps.dic_co2[i] <- (0.91*eps.bicarb_co2[i]) + (0.08*eps.ci_co2[i])
      
      # bottom water carb chem
      eps.bicarb_co2_bot[i] <- -0.1141*tempC_bot[i] + 10.78                   
      eps.ci_co2_bot[i] <- 0.0049*tempC_bot[i] - 1.31
      eps.dic_co2_bot[i] <- (0.91*eps.bicarb_co2_bot[i]) + (0.08*eps.ci_co2_bot[i])
      
      # Benthic forams; equation 7 of Tipple et al. (2010)
      d13Cbf[i] <- ((d13CO2[ai.flat[i]]+1000)*((eps.dic_co2_bot[i]/1000)+1)) - eps.dic_cc - A - 1000
      
      # Bulk carbonate; equation 1 of Barral et al. (2017)
      eps.dic_co2_surf[i] <- coeff1*tempC[i]*f_co2 - coeff2*tempC[i] + eps.int
      d13Cbulk[i] <- d13CO2[ai.flat[i]] + eps.dic_co2_surf[i]
      
      # Brachiopods
      d13Cbrach[i] <- d13CO2[ai.flat[i]] + eps.dic_co2_surf[i]
  }
  
  
  # Time evolution model  
  ####################################################################################################
  # d13C of atmospheric CO2
  d13CO2[1] ~ dnorm(d13CO2.m, 1/d13CO2.sd^2)   
  d13CO2.phi ~ dbeta(5,2) 
  d13CO2.eps[1] = 0 
  d13CO2.tau ~ dgamma(100, 10)
  
  # Global mean surface temperature 
  GMST[1] ~ dnorm(GMST.m[1], 1/GMST.sd[1]^2) 
  
  # Bottom water temperature 
  BWT[1] ~ dnorm(BWT.m[1], 1/BWT.sd[1]^2) 
  
  for (i in 2:n.steps){
    # d13C of atmospheric CO2
    d13CO2.pc[i] <- d13CO2.tau*((1-d13CO2.phi^2)/(1-d13CO2.phi^(2*dt[i-1])))
    d13CO2.eps[i] ~ dnorm(d13CO2.eps[i-1]*(d13CO2.phi^dt[i-1]), d13CO2.pc[i])
    d13CO2[i] <- d13CO2[i-1] + d13CO2.eps[i]
    
    # Global mean surface temperature 
    GMST[i] ~ dnorm(GMST.m[i], 1/GMST.sd[i]^2) 
    
    # Bottom water temperature 
    BWT[i] ~ dnorm(BWT.m[i], 1/BWT.sd[i]^2) 
  }
  
  for (i in 1:length(ai.flat)){
    toff[i] ~ dnorm(toff.m[i], 1/toff.sd[i]^2)
    toff_bot[i] ~ dnorm(0, 1/toff_sd_uniform_bot^2)
    
    tempC[i] <- GMST[ai.flat[i]] + toff[i]
    tempC_bot[i] <- BWT[ai.flat[i]] + toff_bot[i]
  }
}
