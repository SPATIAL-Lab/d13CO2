

# load mcmc objects

load("Phan/chpc_output/ages_PhanDA.rda")

load("Phan/chpc_output/inv.out_LiGMST.rda")
mcmcout_Li <- mcmcout
rm(inv.out)

load("Phan/chpc_output/inv.out_PhanDA.rda")

###################################################################################################
# d13CO2 posterior


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
# GMST posterior


parm <- "GMST"
parm.label <- expression(paste("GMST (",degree,"C)"))
parm.sl_saved <- mcmcout_Li[["BUGSoutput"]][["sims.list"]][["GMST"]]
parm.med_saved <- mcmcout_Li[["BUGSoutput"]][["median"]][["GMST"]]

parm.sl <- inv.out[["BUGSoutput"]][["sims.list"]][["GMST"]]
parm.med <- inv.out[["BUGSoutput"]][["median"]][["GMST"]]

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


plot(ages[10:483]/1e3, parm.med[10:483], type="l", xlab = "age (Ma)", 
     ylab = parm.label, xlim = rev(range(ages[10:549]/1e3)), ylim = c(min(parm.q025[10:483]),max(parm.q975[10:483])))
grid(nx = NA, ny = NULL) 
polygon(x = c(ages[10:483]/1e3, rev(ages[10:483]/1e3)), y = c(parm.q025[10:483], rev(parm.q975[10:483])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
lines(ages[10:483]/1e3, parm.med[10:483], col="dodgerblue", lwd=1)
polygon(x = c(ages[10:483]/1e3, rev(ages[10:483]/1e3)), y = c(parm.q025_saved[10:483], rev(parm.q975_saved[10:483])), 
        col = adjustcolor("black", alpha.f = 0.30), border = NA)
lines(ages[10:483]/1e3, parm.med_saved[10:483], col="black", lwd=1)
legend("topright", legend = c("PhanDA GMST", "Li22 GMST"), fill = c(adjustcolor("dodgerblue", alpha.f=0.40), 
                                                                    adjustcolor("black", alpha.f=0.30)), border = NA, bty="n")



###################################################################################################
# sensitivity test normal vs. uniform distribution for d13CO2 prior 

load("Phan/chpc_output/ages_PhanDA.rda")

load("Phan/chpc_output/inv.out_PhanDA.rda")
mcmcout_norm <- inv.out
rm(inv.out)

load("Phan/chpc_output/inv.out_UniPrior.rda")


parm.label <- expression(paste(delta^13, "C of ", CO[2], " (‰)"))
parm.sl_saved <- mcmcout_norm[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_saved <- mcmcout_norm[["BUGSoutput"]][["median"]][["d13CO2"]]

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
legend("topright", legend = c("uniform d13Ca prior", "normal d13Ca prior"), fill = c(adjustcolor("dodgerblue", alpha.f=0.40),
                                                                    adjustcolor("black", alpha.f=0.30)), border = NA, bty="n")



###################################################################################################
# sensitivity test influence of paleogeographic model

load("Phan/chpc_output/ages_PhanDA.rda")

load("Phan/chpc_output/inv.out_PhanDA.rda")
mcmcout_Sc16 <- inv.out
rm(inv.out)

load("Phan/chpc_output/inv.out_TC17.rda")
mcmcout_TC17 <- inv.out
rm(inv.out)

load("Phan/chpc_output/inv.out_MU22.rda")
mcmcout_MU22 <- inv.out
rm(inv.out)

parm.label <- expression(paste(delta^13, "C of ", CO[2], " (‰)"))
parm.sl_Sc16 <- mcmcout_Sc16[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_Sc16 <- mcmcout_Sc16[["BUGSoutput"]][["median"]][["d13CO2"]]
parm.sl_TC17 <- mcmcout_TC17[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_TC17 <- mcmcout_TC17[["BUGSoutput"]][["median"]][["d13CO2"]]
parm.sl_MU22 <- mcmcout_MU22[["BUGSoutput"]][["sims.list"]][["d13CO2"]]
parm.med_MU22 <- mcmcout_MU22[["BUGSoutput"]][["median"]][["d13CO2"]]

parm.q025_Sc16 <- vector(mode = "numeric", length = length(parm.sl_Sc16[1,]))
parm.q975_Sc16 <- vector(mode = "numeric", length = length(parm.sl_Sc16[1,]))
parm.q025_TC17 <- vector(mode = "numeric", length = length(parm.sl_TC17[1,]))
parm.q975_TC17 <- vector(mode = "numeric", length = length(parm.sl_TC17[1,]))
parm.q025_MU22 <- vector(mode = "numeric", length = length(parm.sl_MU22[1,]))
parm.q975_MU22 <- vector(mode = "numeric", length = length(parm.sl_MU22[1,]))


for (i in 1:length(parm.sl_Sc16[1,])){
  parm.q025_Sc16[i] <- quantile(parm.sl_Sc16[,i], 0.025)
  parm.q975_Sc16[i] <- quantile(parm.sl_Sc16[,i], 0.975)
}

for (i in 1:length(parm.sl_TC17[1,])){
  parm.q025_TC17[i] <- quantile(parm.sl_TC17[,i], 0.025)
  parm.q975_TC17[i] <- quantile(parm.sl_TC17[,i], 0.975)
}

for (i in 1:length(parm.sl_MU22[1,])){
  parm.q025_MU22[i] <- quantile(parm.sl_MU22[,i], 0.025)
  parm.q975_MU22[i] <- quantile(parm.sl_MU22[,i], 0.975)
}


plot(ages[10:549]/1e3, parm.med_Sc16[10:549], type="l", xlab = "age (Ma)",
     ylab = parm.label, xlim = rev(range(ages[10:549]/1e3)), ylim = c(min(parm.q025_Sc16[10:549]),max(parm.q975_Sc16[10:549])))
grid(nx = NA, ny = NULL)
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025_Sc16[10:549], rev(parm.q975_Sc16[10:549])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
lines(ages[10:549]/1e3, parm.med_Sc16[10:549], col="black", lwd=1)
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025_TC17[10:549], rev(parm.q975_TC17[10:549])),
        col = adjustcolor("coral", alpha.f = 0.40), border = NA)
lines(ages[10:549]/1e3, parm.med_TC17[10:549], col="coral", lwd=1)
polygon(x = c(ages[10:549]/1e3, rev(ages[10:549]/1e3)), y = c(parm.q025_MU22[10:549], rev(parm.q975_MU22[10:549])),
        col = adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
lines(ages[10:549]/1e3, parm.med_MU22[10:549], col="dodgerblue", lwd=1)
legend("topright", legend = c("Scotese (2016)", "Torvik & Cocks (2017)", "Müller et al. (2022)"), fill = c(adjustcolor("black", alpha.f=0.30),
                                                                                     adjustcolor("coral", alpha.f=0.40),
                                                                                     adjustcolor("dodgerblue", alpha.f=0.40)), border = NA, bty="n")



