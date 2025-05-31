
# load("mcmcout_Li.RData")
# load("mcmcout_PhanDA_normal.RData")

mcmcout_Li <- mcmc_out

parm <- "d13CO2"
parm.label <- expression(paste(delta^13, "C of ", CO[2], " (‰)"))
parm.sl_saved <- mcmcout_Li[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_saved <- mcmcout_Li[["BUGSoutput"]][["median"]][["d13CO2"]]

parm.sl <- inv.out[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med <- inv.out[["BUGSoutput"]][["median"]][["d13CO2"]]

parm.q025_saved <- vector(mode = "numeric", length = length(parm.sl_saved[1,]))
parm.q975_saved <- vector(mode = "numeric", length = length(parm.sl_saved[1,]))

parm.q025 <- vector(mode = "numeric", length = length(parm.sl[1,]))
parm.q975 <- vector(mode = "numeric", length = length(parm.sl[1,]))

for (i in 1:length(parm.sl[1,])){
  parm.q025[i] <- quantile(parm.sl[,i], 0.025)
  parm.q975[i] <- quantile(parm.sl[,i], 0.975)
}

for (i in 1:length(parm.sl_saved[1,])){
  parm.q025_saved[i] <- quantile(parm.sl_saved[,i], 0.025)
  parm.q975_saved[i] <- quantile(parm.sl_saved[,i], 0.975)
}


plot(ages[10:549]/1e3, parm.med[10:549], type="l", xlab = "age (Ma)", 
     ylab = parm.label, xlim = rev(range(ages[10:549]/1e3)), ylim = c(min(parm.q025[10:549]),max(parm.q975[10:549])))
grid(nx = NA, ny = NULL) 
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025[10:549], rev(parm.q975[10:549])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
lines(ages[10:549]/1e3, parm.med[10:549], col="dodgerblue", lwd=1)
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025_saved[10:549], rev(parm.q975_saved[10:549])), 
        col = adjustcolor("black", alpha.f = 0.30), border = NA)
lines(ages[10:549]/1e3, parm.med_saved[10:549], col="black", lwd=1)
legend("topright", legend = c("PhanDA GMST", "Li22 GMST"), fill = c(adjustcolor("dodgerblue", alpha.f=0.40), 
       adjustcolor("black", alpha.f=0.30)), border = NA, bty="n")



###################################################################################################

mcmcout_norm <- mcmcout

parm.sl_saved <- mcmcout_norm[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_saved <- mcmcout_norm[["BUGSoutput"]][["median"]][["d13CO2"]]

parm.q025_saved <- vector(mode = "numeric", length = length(parm.sl_saved[1,]))
parm.q975_saved <- vector(mode = "numeric", length = length(parm.sl_saved[1,]))

for (i in 1:length(parm.sl_saved[1,])){
  parm.q025_saved[i] <- quantile(parm.sl_saved[,i], 0.025)
  parm.q975_saved[i] <- quantile(parm.sl_saved[,i], 0.975)
}


plot(ages[10:549]/1e3, parm.med[10:549], type="l", xlab = "age (Ma)", 
     ylab = parm.label, xlim = rev(range(ages[10:549]/1e3)), ylim = c(min(parm.q025[10:549]),max(parm.q975[10:549])))
grid(nx = NA, ny = NULL) 
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025[10:549], rev(parm.q975[10:549])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
lines(ages[10:549]/1e3, parm.med[10:549], col="dodgerblue", lwd=1)
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025_saved[10:549], rev(parm.q975_saved[10:549])), 
        col = adjustcolor("black", alpha.f = 0.30), border = NA)
lines(ages[10:549]/1e3, parm.med_saved[10:549], col="black", lwd=1)
legend("topright", legend = c("uniform d13Ca prior", "normal d13Ca prior"), fill = c(adjustcolor("dodgerblue", alpha.f=0.40), 
                                                                    adjustcolor("black", alpha.f=0.30)), border = NA, bty="n")





