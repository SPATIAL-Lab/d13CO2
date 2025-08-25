
# Plots time series, or individual sample posteriors from MCMC output object


### INPUT ###
############################################################################################
### MCMC object from MCMC inversion output 
mcmcout <- inv.out


### List of parameters to plot. Parameter names must match saved parameters in MCMC object
parms2plot <- c("d13CO2", "GMST")
format.parm.labels <- TRUE

if (format.parm.labels == TRUE){
  parm.labels <- c(expression(paste(delta^13, "C of ", CO[2], " (‰)")), 
                   (expression(paste("GMST (", degree,"C)"))))

} else if (format.parm.labels == FALSE){
  parm.labels <- parms2plot
}


## Output the table of sites, site index #s
print(sites)

############################################################################################


### PLOT THE RAW d13C DATA, COLORED BY SITE ###
############################################################################################
plot(prox.in$age, prox.in$d13C, col=site.index, xlim = c(age.max, age.min), cex = 0.3)

############################################################################################


### PULL DATA FROM MCMC OBJECT AND PLOT SELECTED PARAMETERS ###
############################################################################################
for (i in 1:length(parms2plot)){
  parm <- parms2plot[i]
  parm.label <- parm.labels[i]
  parm.sl <- mcmcout[["BUGSoutput"]][["sims.list"]][[parm]]
  parm.med <- mcmcout[["BUGSoutput"]][["median"]][[parm]]
  
  if (length(dim(parm.sl)) > 2){
    col.pal <- viridis(9)
    plot(ages/1e3, parm.med[,1], type="l", xlab = "age (Ma)", col =adjustcolor("dodgerblue", alpha.f = 0),
         ylab = parm.label, xlim = rev(range(ages/1e3)), ylim = c(min(parm.med),max(parm.med)))
    grid(nx = NA, ny = NULL) 
    for (i in 1:length(parm.med[1,])){
      lines(ages/1e3, parm.med[,i], col=col.pal[i], lwd=1)
    }
  }
  
  if (length(parm.med) < 2){
    plot(density(parm.sl), xlab = parm.label, col = rgb(0,0,0, alpha = 0))
    grid(nx = NA, ny = NULL) 
    polygon(density(parm.sl), col = adjustcolor("dodgerblue", alpha.f = 0.40))
    abline(v=parm.med, col="black")
    
  } else if (length(dim(parm.sl)) < 3){
      parm.q025 <- vector(mode = "numeric", length = length(parm.sl[1,]))
      parm.q250 <- vector(mode = "numeric", length = length(parm.sl[1,]))
      parm.q750 <- vector(mode = "numeric", length = length(parm.sl[1,]))
      parm.q975 <- vector(mode = "numeric", length = length(parm.sl[1,]))
      
      for (i in 1:length(parm.sl[1,])){
        parm.q025[i] <- quantile(parm.sl[,i], 0.025)
        parm.q250[i] <- quantile(parm.sl[,i], 0.250)
        parm.q750[i] <- quantile(parm.sl[,i], 0.750)
        parm.q975[i] <- quantile(parm.sl[,i], 0.975)
      }
    }
  if (length(dim(parm.sl)) < 3){
    plot(ages/1e3, parm.med, type="l", xlab = "age (Ma)", 
         ylab = parm.label, xlim = rev(range(ages/1e3)), ylim = c(min(parm.q025),max(parm.q975)))
    grid(nx = NA, ny = NULL) 
    polygon(x = c(ages/1e3, rev(ages/1e3)), y = c(parm.q025, rev(parm.q975)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
    polygon(x = c(ages/1e3, rev(ages/1e3)), y = c(parm.q250, rev(parm.q750)), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.80), border = NA)
    lines(ages/1e3, parm.med, col=rgb(red=0, green=0, blue=0, alpha=1), lwd=1)
  }  
}  


############################################################################################
# Density plots of individual timesteps of d13CO2

############################################################################################

timestep = 17

plot(x = c(-13,1), y =c(0,1.2), col= NULL, xlab = expression(paste(delta^{13}, "C (\u2030 VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, -6, 2)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$d13CO2[,timestep]),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$d13CO2[timestep],inv.out$BUGSoutput$median$d13CO2[timestep]), y = c(0,1.2),
      col = "black", lty = "dashed")
lines(x = c(-6,-6), y = c(0,0.2), col = "black", lty = "dashed")





