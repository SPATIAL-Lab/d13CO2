model {
  
  # Likelihood block 
  ####################################################################################################
  

  # Benthic forams
  d13Cbf.analyt.sd <- 0.05
  d13Cbf.p <- 1/d13Cbf.analyt.sd^2
  for (i in 1:n.d13Cbf){
    d13Cbf.data[i] ~ dnorm(d13Cbf[ri.d13Cbf[i]], d13Cbf.p)
  }

  # micrite open ocean
  d13Cmicrite.analyt.sd <- 0.05
  d13Cmicrite.p <- 1/(d13Cmicrite.analyt.sd)^2
  for (i in 1:n.d13Cmicrite){
    d13Cmicrite.data[i] ~ dnorm(d13Cmicrite[ri.d13Cmicrite[i]], d13Cmicrite.p)
  }

  # Bulk carbonate open ocean 
  d13Cbulk.analyt.sd <- 0.05
  d13Cbulk.p <- 1/d13Cbulk.analyt.sd^2
  for (i in 1:n.d13Cbulk){
    d13Cbulk.data[i] ~ dnorm(d13Cbulk[ri.d13Cbulk[i]], d13Cbulk.p)
  }

  # Bulk carbonate semi restricted
  d13Cbulk_sr.analyt.sd <- 0.05
  d13Cbulk_sr.p <- 1/(d13Cbulk_sr.analyt.sd)^2
  for (i in 1:n.d13Cbulk_sr){
    d13Cbulk_sr.data[i] ~ dnorm(d13Cbulk_sr[ri.d13Cbulk_sr[i]], d13Cbulk_sr.p)
  }

  # Bulk carbonate marginal sea
  d13Cbulk_marg.analyt.sd <- 0.05
  d13Cbulk_marg.p <- 1/(d13Cbulk_marg.analyt.sd)^2
  for (i in 1:n.d13Cbulk_marg){
    d13Cbulk_marg.data[i] ~ dnorm(d13Cbulk_marg[ri.d13Cbulk_marg[i]], d13Cbulk_marg.p)
  }


  # Constants
  ####################################################################################################  
  # Equation 8 of Tipple et al. (2010)
  eps.dic_cc <- -1   
  
  # Bulk carbonate uncertainty terms from Romanek et al. (1992) 
  cc_co2_constant1 ~ dnorm(11.98, 1/0.13^2)T(11.72,12.24)
  cc_co2_coeff1 ~ dnorm(0.12, 1/0.01^2)T(0.1,0.14)
  
  # Benthic disequilibrium effects; totals to 2.8 (mean) from Tipple et al. (2010)
  asd ~ dnorm(1, 1/0.2^2)T(0,)  # air sea disequilibrium 
  bpump ~ dnorm(1.2, 1/0.4^2)T(0,) # biological pump
  remin ~ dnorm(0.6, 1/0.3^2)T(0,) # remineralization and oxidation 
  A <- asd + bpump + remin
  f_co3 ~ dnorm(0.12, 1/0.04^2)T(0.04,0.20)
  f_carbacid ~ dnorm(0.01, 1/0.005^2)T(0,0.02)
  
  # Non-secular bias uncertainty - hyperpriors set the across-site scatter (partial pooling)
  bf.nsb_sigma ~ dunif(0.3*bf.nsb.sd, bf.nsb.sd)         
  bf.nsb_tau <- 1 / (bf.nsb_sigma^2)
  bulk.nsb_sigma ~ dunif(0.3*bulk.nsb.sd, bulk.nsb.sd)         
  bulk.nsb_tau <- 1 / (bulk.nsb_sigma^2)
  micrite.nsb_sigma ~ dunif(0.3*micrite.nsb.sd, micrite.nsb.sd)         
  micrite.nsb_tau <- 1 / (micrite.nsb_sigma^2)
  bulk_sr.nsb_sigma ~ dunif(0.3*bulk_sr.nsb.sd, bulk_sr.nsb.sd)         
  bulk_sr.nsb_tau <- 1 / (bulk_sr.nsb_sigma^2)
  bulk_marg.nsb_sigma ~ dunif(0.3*bulk_marg.nsb.sd, bulk_marg.nsb.sd)         
  bulk_marg.nsb_tau <- 1 / (bulk_marg.nsb_sigma^2)

  # Site-level NSBs
  for (i in 1:n.sites){
    bf.nsb_site[i] ~ dt(bf.nsb.m, bf.nsb_tau, 3)
    bulk.nsb_site[i] ~ dt(bulk.nsb.m, bulk.nsb_tau, 3)
    micrite.nsb_site[i] ~ dt(micrite.nsb.m, micrite.nsb_tau, 3)
    bulk_sr.nsb_site[i] ~ dt(bulk_sr.nsb.m, bulk_sr.nsb_tau, 3)
    bulk_marg.nsb_site[i] ~ dt(bulk_marg.nsb.m, bulk_marg.nsb_tau, 3)
  }
  
  # Proxy system model 
  ####################################################################################################
  for (i in 1:length(ai.flat)){
    
      # Equations 4, 5, and 3 (respectively) of Tipple et al. (2010); bottom water carb chem fractionation
      eps.bicarb_co2_bot[i] <- -0.1141*tempC_bot[i] + 10.78                   
      eps.ci_co2_bot[i] <- 0.0049*tempC_bot[i] - 1.31
      eps.dic_co2_bot[i] <- (1-f_co3-f_carbacid)*(eps.bicarb_co2_bot[i]) + (f_co3*eps.ci_co2_bot[i])
      
      
      # Romanek et al. (1992) bulk calcite archives
      eps.cc_co2_surf[i] <- cc_co2_constant1 - cc_co2_coeff1*tempC[i]

      # Calculate various carbonate archive d13C values
      d13Cbulk_secular[i] <- d13CO2[ai.flat[i]] + eps.cc_co2_surf[i]
      d13Cbulk[i] <- d13Cbulk_secular[i] + bulk.nsb_site[si.flat[i]]
      d13Cmicrite[i] <- d13Cbulk_secular[i] + micrite.nsb_site[si.flat[i]]
      d13Cbulk_sr[i] <- d13Cbulk_secular[i] + bulk_sr.nsb_site[si.flat[i]]
      d13Cbulk_marg[i] <- d13Cbulk_secular[i] + bulk_marg.nsb_site[si.flat[i]]
      
      # Benthic forams; equation 7 of Tipple et al. (2010)
      d13Cbf[i] <- ((d13CO2[ai.flat[i]]+1000)*((eps.dic_co2_bot[i]/1000)+1)) - eps.dic_cc 
                    - A - 1000 + bf.nsb_site[si.flat[i]]
  }
  
  
  # Time evolution model  
  ####################################################################################################
  d13CO2[1] ~ dnorm(d13CO2.m, 1/d13CO2.sd^2)
  d13CO2_sigma ~ dunif(0, 0.5)         
  d13CO2_tau <- 1 / (d13CO2_sigma^2)
  
  GMST[1] ~ dnorm(GMST.m[1], 1/GMST.sd[1]^2) 
  BWT[1] ~ dnorm(BWT.m[1], 1/BWT.sd[1]^2) 
  
  for (i in 2:n.steps){
    d13CO2[i] ~ dnorm(d13CO2[i-1], d13CO2_tau)
    GMST[i] ~ dnorm(GMST.m[i], 1/GMST.sd[i]^2) 
    BWT[i] ~ dnorm(BWT.m[i], 1/BWT.sd[i]^2) 
  }
  
  for (i in 1:length(ai.flat)){
    toff[i] ~ dnorm(toff.m[i], 1/toff.sd[i]^2)
    toff_bot[i] ~ dnorm(0, 1/toff_sd_uniform_bot^2)
    
    tempC[i] <- GMST[ai.flat[i]] + toff[i]
    tempC_bot[i] <- BWT[ai.flat[i]] + toff_bot[i]
  }
}
