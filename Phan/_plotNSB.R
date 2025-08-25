
# non-secular bias of archival d13C from high-fidelity (true) value (per mille), standard deviation w/ fixed mean = 0 

# Benthic forams
bf.nsb.m <- 0
bf.nsb.sd <- 0.5

plot(x = c(-2,2), y =c(0,3.2), col= NULL, xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, bf.nsb.m, bf.nsb.sd)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$bf.nsb),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$bf.nsb,inv.out$BUGSoutput$median$bf.nsb), y = c(0,1),
      col = "darkred", lty = "dashed")
lines(x = c(bf.nsb.m,bf.nsb.m), y = c(0,0.78), col = "navyblue", lty = "dashed")
legend("topright", legend = c("b. foram NSB prior", "b. foram NSB post"), 
       fill = c(adjustcolor("slateblue", alpha.f=0.50), adjustcolor("coral", alpha.f=0.50)), border = NA, bty="n")



# Bulk carbonate - open ocean
bulk.nsb.m <- 0
bulk.nsb.sd <- 0.5

plot(x = c(-2,2), y =c(0,3.2), col= NULL, xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, bulk.nsb.m, bulk.nsb.sd)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$bulk.nsb),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$bulk.nsb,inv.out$BUGSoutput$median$bulk.nsb), y = c(0,1.85),
      col = "darkred", lty = "dashed")
lines(x = c(bulk.nsb.m,bulk.nsb.m), y = c(0,0.8), col = "navyblue", lty = "dashed")
legend("topright", legend = c("bulk open ocean NSB prior", "bulk open ocean NSB post"), 
       fill = c(adjustcolor("slateblue", alpha.f=0.50), adjustcolor("coral", alpha.f=0.50)), border = NA, bty="n")



# micrite - open ocean
micrite.nsb.m <- 0
micrite.nsb.sd <- 0.25

plot(x = c(-1,1), y =c(0,3), col= NULL, xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, micrite.nsb.m, micrite.nsb.sd)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$micrite.nsb),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$micrite.nsb,inv.out$BUGSoutput$median$micrite.nsb), y = c(0,1.9),
      col = "darkred", lty = "dashed")
lines(x = c(micrite.nsb.m,micrite.nsb.m), y = c(0,1.56), col = "navyblue", lty = "dashed")
legend("topright", legend = c("micrite open ocean NSB prior", "micrite open ocean NSB post"), 
       fill = c(adjustcolor("slateblue", alpha.f=0.50), adjustcolor("coral", alpha.f=0.50)), border = NA, bty="n")



# Bulk carbonate - semi-restricted
bulk_sr.nsb.m <- -1
bulk_sr.nsb.sd <- 1

plot(x = c(-4,2), y =c(0,1.8), col= NULL, xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, bulk_sr.nsb.m, bulk_sr.nsb.sd)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$bulk_sr.nsb),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$bulk_sr.nsb,inv.out$BUGSoutput$median$bulk_sr.nsb), y = c(0,1.45),
      col = "darkred", lty = "dashed")
lines(x = c(bulk_sr.nsb.m,bulk_sr.nsb.m), y = c(0,0.38), col = "navyblue", lty = "dashed")
legend("topright", legend = c("bulk semi-restricted NSB prior", "bulk semi-restricted NSB post"), 
       fill = c(adjustcolor("slateblue", alpha.f=0.50), adjustcolor("coral", alpha.f=0.50)), border = NA, bty="n")



# Bulk carbonate - marginal sea
bulk_marg.nsb.m <- -0.5
bulk_marg.nsb.sd <- 1

plot(x = c(-3,2), y =c(0,2), col= NULL, xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")
polygon(density(rnorm(100000, bulk_marg.nsb.m, bulk_marg.nsb.sd)), col =adjustcolor("slateblue", alpha.f = 0.50), border = NA)
polygon(density(inv.out$BUGSoutput$sims.list$bulk_marg.nsb),
        col = adjustcolor("coral1", alpha.f = 0.50), border = NA, title = NULL)
lines(x = c(inv.out$BUGSoutput$median$bulk_marg.nsb,inv.out$BUGSoutput$median$bulk_marg.nsb), y = c(0,1.5),
      col = "darkred", lty = "dashed")
lines(x = c(bulk_marg.nsb.m ,bulk_marg.nsb.m ), y = c(0,0.37), col = "navyblue", lty = "dashed")
legend("topright", legend = c("bulk marginal NSB prior", "bulk marginal NSB post"), 
       fill = c(adjustcolor("slateblue", alpha.f=0.50), adjustcolor("coral", alpha.f=0.50)), border = NA, bty="n")





