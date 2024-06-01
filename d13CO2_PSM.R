model {
  
# Likelihood block 
####################################################################################################
  
  for (i in 1:length(d13Cc.data)){
    d13Cc.data[i] ~ dnorm(d13Cc[ai.d13Cc[i], si.d13Cc[i]], d13Cc.p[ai.d13Cc[i]])
    d13Cc.p[i] = 1/d13Cc.sd[i]^2
  }
  
  
# Proxy system model 
####################################################################################################
  for (i in 1:n.steps){
    for (j in 1:n.sites){
      # Calculate carbonate d13C from env parm d13CO2; rearranged eqn of Romanek et al. (1992)
      # UPDATE THIS RELATIONSHIP
      d13Cc[i,j] <- ((11.98 - 0.12*tempC[i,j])/1000 + 1)*(d13CO2[i] + 1000) - 1000 
    }
  }
  
  
# Time evolution model  
####################################################################################################
  # d13C of atmospheric CO2
  d13CO2[1] ~ dnorm(d13CO2.m, 1/d13CO2.sd^2)T(-10,-6)    
  d13CO2.phi ~ dbeta(5,2) 
  d13CO2.eps[1] = 0 
  d13CO2.tau ~ dgamma(1e3, 1e-3)
  
  # Temp in C
  for (j in 1:n.sites){
    toff[1,j] ~ dnorm(toff.m[1,j], toff.sd[1,j])
    GMST[1,j] ~ dnorm(GMST.m[1], 1/GMST.sd[1]^2) 
    tempC[1,j] = GMST[1,j] + toff[1,j]
  }
  
  for (i in 2:n.steps){
    # d13C of atmospheric CO2
    d13CO2.pc[i] <- d13CO2.tau*((1-d13CO2.phi^2)/(1-d13CO2.phi^(2*dt[i-1])))
    d13CO2.eps[i] ~ dnorm(d13CO2.eps[i-1]*(d13CO2.phi^dt[i-1]), d13CO2.pc[i])T(-2, 2)
    d13CO2[i] <- d13CO2[1] + d13CO2.eps[i]
    
    for (j in 1:n.sites){
      # Temp in C
      toff[i,j] ~ dnorm(toff.m[i,j], toff.sd[i,j])
      GMST[i,j] ~ dnorm(GMST.m[i], 1/GMST.sd[i]^2) 
      tempC[i,j] = GMST[i,j] + toff[i,j]
    }
  }
}






