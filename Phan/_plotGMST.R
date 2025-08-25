
age.min <- 0
age.max <- 540000
step.int <- 1000
n.spinup <- 10

# uniform standard deviation applied to GMST from Li et al (2022); they do not provide uncertainty estimates
GMST_sd_Li22 <- 5
# uniform standard deviation applied to all surface temperature site offset values - no uncertainty in Li22 estimates
toff_sd_uniform <- 2 
# uniform standard deviation applied to all bottom temperature site offset values 
toff_sd_uniform_bot <- 1 


# load proxy data
prox.in <- as.data.frame(read.csv(file = "Phan/PhanData/PhanCompWithTemp.csv"))
prox.in <- cbind(prox.in[,3:18], rep(x = toff_sd_uniform, times = nrow(prox.in)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category", 
                    "paleolon","paleolat", "MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                    "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA", "temp_offset_sd")


# age index proxy data (in kyrs)
prox.in$age <- prox.in$age*1e3
prox.in <- prox.in[prox.in$age >= (age.min) & prox.in$age <= (age.max),]
prox.in <- transform(prox.in, ai = n.spinup + as.numeric(1 + floor((max(prox.in$age) - prox.in$age) / step.int)))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE),]
ages.short <- seq(from = max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
ages <- seq(from = n.spinup*step.int + max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int

# calculate symmetric uncertainty (GMST_PhanDA_sd)
PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                (prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo)) / 2


# interpolate GMST for ages 
GMST.m <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = ages, rule = 2)$y
GMST.sd <- approx(prox.in$age, PhanDA_sd, xout = ages, rule = 2)$y
GMST.low <- GMST.m - GMST.sd*2
GMST.hi <- GMST.m + GMST.sd*2

GMST.m.Li <- approx(prox.in$age, prox.in$GMST_Li22, xout = ages, rule = 2)$y
GMST.sd.Li <- rep(x = GMST_sd_Li22, times = length(ages)) 
GMST.low.Li <- GMST.m.Li - GMST.sd.Li*2
GMST.hi.Li <- GMST.m.Li + GMST.sd.Li*2

# plot 
parm.label <- expression(paste("GMST (",degree,"C)"))
plot(ages[10:483]/1e3, GMST.m[10:483], type="l", xlab = "age (Ma)", 
     ylab = parm.label, xlim = rev(range(ages[10:549]/1e3)), ylim = c(-3,56))
grid(nx = NA, ny = NULL) 
polygon(x = c(ages[10:483]/1e3, rev(ages[10:483]/1e3)), y = c(GMST.low[10:483], rev(GMST.hi[10:483])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.40), border = NA)
lines(ages[10:483]/1e3, GMST.m[10:483], col="dodgerblue", lwd=1)
polygon(x = c(ages[10:483]/1e3, rev(ages[10:483]/1e3)), y = c(GMST.low.Li[10:483], rev(GMST.hi.Li[10:483])), 
        col = adjustcolor("black", alpha.f = 0.30), border = NA)
lines(ages[10:483]/1e3, GMST.m.Li[10:483], col="black", lwd=1)
legend("topright", legend = c("PhanDA GMST prior", "Li22 GMST prior"), fill = c(adjustcolor("dodgerblue", alpha.f=0.40), 
                                                                    adjustcolor("black", alpha.f=0.30)), border = NA, bty="n")



