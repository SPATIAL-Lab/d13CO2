model {
  
# Likelihood block 
####################################################################################################
  
  for (i in 1:length(d13Cc.data)){
    d13Cc.data[i] ~ dnorm(d13Cc[ai.d13Cc[i], si.d13Cc[i]], d13Cc.p[ai.d13Cc[i]])
    d13Cc.p[i] = 1/d13Cc.sd[i]^2
  }
  
  
# Proxy System Model 
####################################################################################################
  for (i in 1:n.steps){
    for (j in 1:n.sites){
      # Calculate foram d13C from env parm d13C.co2; rearranged eqns of Romanek et al. (1992) and Mook et al. (1974)
      d13Cc[i,j] <- ((11.98 - 0.12*tempC[i,j])/1000 + 1)*(d13C.co2[i] + 1000) - 1000 
    }
  }
  
  
# Time evolution model  
####################################################################################################
  # d13C of atmospheric CO2
  d13C.co2[1] ~ dnorm(d13C.co2.m, 1/d13C.co2.sd^2)T(-10,-6)    
  d13C.co2.phi ~ dbeta(5,2) 
  d13C.co2.eps[1] = 0 
  d13C.co2.tau ~ dgamma(1e3, 1e-3)
  
  # Temp in C
  for (j in 1:n.sites){
    t.off[1,j] ~ dnorm(t.off.m[1,j], t.off.sd[1,j])
    GMST[1,j] ~ dnorm(GMST.m[1], 1/GMST.sd[1]^2) 
    tempC[1,j] = GMST[1,j] + t.off[1,j]
  }
  
  for (i in 2:n.steps){
    # d13C of atmospheric CO2
    d13C.co2.pc[i] <- d13C.co2.tau*((1-d13C.co2.phi^2)/(1-d13C.co2.phi^(2*dt[i-1])))
    d13C.co2.eps[i] ~ dnorm(d13C.co2.eps[i-1]*(d13C.co2.phi^dt[i-1]), d13C.co2.pc[i])T(-2, 2)
    d13C.co2[i] <- d13C.co2[1] + d13C.co2.eps[i]
    
    for (j in 1:n.sites){
      # Temp in C
      t.off[i,j] ~ dnorm(t.off.m[i,j], t.off.sd[i,j])
      GMST[i,j] ~ dnorm(GMST.m[i], 1/GMST.sd[i]^2) 
      tempC[i,j] = GMST[i,j] + t.off[i,j]
    }
  }
}






