
# Driver script for PSM used to reconstruct d13C of atmospheric CO2 
# Dustin T. Harper


# Load libraries
############################################################################################
library(R2jags)

# Load and prepare proxy data 
############################################################################################
# Load proxy data
prox.in <- read.csv(file = "data/d13Ccarb.csv")
names(prox.in) <- c("site", "age", "d13C", "d13C.sd")

# Age and site index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE),]
prox.in <- transform(prox.in,site.index=as.numeric(factor(site)))
site.index <- c(prox.in$site.index)
n.sites <- as.numeric(length(unique(site.index)))
ages.prox <- sort(unique(prox.in$age), decreasing = TRUE)
dt <- abs(diff(ages.prox, lag=1))

# Clean and prepare proxy data 
clean.d13Cc <- prox.in[complete.cases(prox.in$d13C), ]
ai.d13Cc <- sort(c(as.integer(clean.d13Cc$ai)), decreasing = FALSE)    
n.steps <- as.numeric(length((ai.d13Cc)))
si.d13Cc <- clean.d13Cc$site.index
d13Cc.data <- clean.d13Cc$d13C
d13Cc.sd <- clean.d13Cc$d13C.sd
############################################################################################


# Load GMST and temperature offset matrix
############################################################################################
# GMST interpolated for time steps 
GMST.in <- as.data.frame(read.csv(file = "data/GMST.csv"))
GMST.m <- approx(GMST.in$age, GMST.in$GMST.m, xout=ages.prox, method="linear")
GMST.sd <- approx(GMST.in$age, GMST.in$GMST.sd, xout=ages.prox, method="linear")
GMST.m <- GMST.m[["y"]]
GMST.sd <- GMST.sd[["y"]]

# Matrix of spatio-temporally dependent temperature offsets from GMST
t.off.m.in <- as.data.frame(read.csv(file = "data/TempOffsetMean.csv")) # matrix [i,j] of mean values for temperature offset, for ages.prox (i) and sites (j) 
t.off.m = c(1:(length(ages.prox)*(length(t.off.m.in[,1])-1)))
dim(t.off.m) = c((length(ages.prox)), (length(t.off.m.in[,1])-1))
for (i in 2:length(t.off.m.in)){
  t.off.m[,i-1] <- approx(x = t.off.m.in$age, y = as.numeric(unlist(t.off.m.in[i])), xout=ages.prox, method="linear")[["y"]]
}

t.off.sd.in <- as.data.frame(read.csv(file = "data/TempOffsetSD.csv"))  # matrix [i,j] of sds for temperature offset, for ages.prox (i) and sites (j)
t.off.sd = c(1:(length(ages.prox)*(length(t.off.sd.in[,1])-1)))
dim(t.off.sd) = c((length(ages.prox)), (length(t.off.sd.in[,1])-1))
for (i in 2:length(t.off.sd.in)){
  t.off.sd[,i-1] <- approx(x = t.off.sd.in$age, y = as.numeric(unlist(t.off.sd.in[i])), xout=ages.prox, method="linear")[["y"]]
}
############################################################################################


# Environmental prior (d13C of atm CO2) 
############################################################################################
d13C.co2.m <- -8
d13C.co2.sd <- 2
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass = list("d13Cc.data" = d13Cc.data,   
                 "d13Cc.sd" = d13Cc.sd,   
                 "n.steps" = n.steps,
                 "dt" = dt,
                 "n.sites" = n.sites,
                 "ai.d13Cc" = ai.d13Cc,
                 "si.d13Cc" = si.d13Cc,
                 "GMST.m" = GMST.m,
                 "GMST.sd" = GMST.sd,
                 "t.off.m" = t.off.m,
                 "t.off.sd" = t.off.sd,
                 "d13C.co2.m" = d13C.co2.m,
                 "d13C.co2.sd" = d13C.co2.sd) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("tempC", "d13C.co2", "d13Cc")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "d13CO2_PSM.R", parameters.to.save = parms,
                        inits = NULL, n.chains = 3, n.iter = 1e4,
                        n.burnin = 5e3, n.thin = 1)

############################################################################################

  
  
  
  
  