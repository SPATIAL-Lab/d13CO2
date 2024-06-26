
# Driver script for PSM used to reconstruct d13C of atmospheric CO2 
# Dustin T. Harper


# Load libraries
############################################################################################
library(R2jags)

# Load and prepare proxy data 
############################################################################################
# Load proxy data
prox.in <- read.csv(file = "data/d13Ccarb.csv")
names(prox.in) <- c("site", "age", "d13C", "d13C.sd", "archive")

# Age and site index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1))))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE),]
prox.in <- transform(prox.in,site.index=as.numeric(factor(site, ordered = is.ordered(site))))
site.index <- c(prox.in$site.index)
n.sites <- as.numeric(length(unique(site.index)))
ages.prox <- sort(unique(prox.in$age), decreasing = TRUE)
dt <- abs(diff(ages.prox, lag=1))

# List the input sites and their associated site index
sites <- data.frame((sort(unique(prox.in$site), decreasing = FALSE)),seq(1:length(sort(unique(prox.in$site), decreasing = FALSE))))
names(sites) <- c("site", "site.index")
print(sites)
print("Make sure the site indexes in the temp offset .csv files are ordered as they are here! (Numerical increasing)")

# Clean and prepare proxy data 
clean.d13C <- prox.in[complete.cases(prox.in$d13C), ]

# planktic foraminifera
clean.d13Cpf <- clean.d13C[clean.d13C$archive == "pf",]
ai.d13Cpf <- sort(c(as.integer(clean.d13Cpf$ai)), decreasing = FALSE)    
si.d13Cpf <- clean.d13Cpf$site.index
d13Cpf.data <- clean.d13Cpf$d13C
d13Cpf.sd <- clean.d13Cpf$d13C.sd

# benthic foraminifera 
clean.d13Cbf <- clean.d13C[clean.d13C$archive == "bf",]
ai.d13Cbf <- sort(c(as.integer(clean.d13Cbf$ai)), decreasing = FALSE)    
si.d13Cbf <- clean.d13Cbf$site.index
d13Cbf.data <- clean.d13Cbf$d13C
d13Cbf.sd <- clean.d13Cbf$d13C.sd

# bulk carbonate 
clean.d13Cbulk <- clean.d13C[clean.d13C$archive == "bulk",]
ai.d13Cbulk <- sort(c(as.integer(clean.d13Cbulk$ai)), decreasing = FALSE)    
si.d13Cbulk <- clean.d13Cbulk$site.index
d13Cbulk.data <- clean.d13Cbulk$d13C
d13Cbulk.sd <- clean.d13Cbulk$d13C.sd

# parameters for all proxies
ai.all <- sort((c(ai.d13Cpf, ai.d13Cbf, ai.d13Cbulk)) , decreasing = FALSE)
n.steps <- as.numeric(length((ai.all)))
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
INPUTtoff.m <- as.data.frame(read.csv(file = "data/TempOffsetMean.csv")) # matrix [i,j] of mean values for temperature offset, for ages.prox (i) and sites (j) 
toff.m = c(1:(length(ages.prox)*(length(INPUTtoff.m[,1])-1)))
dim(toff.m) = c((length(ages.prox)), (length(INPUTtoff.m[,1])-1))
for (i in 2:length(INPUTtoff.m)){
  toff.m[,i-1] <- approx(x = INPUTtoff.m$age, y = as.numeric(unlist(INPUTtoff.m[i])), xout=ages.prox, method="linear")[["y"]]
}

INPUTtoff.sd <- as.data.frame(read.csv(file = "data/TempOffsetSD.csv"))  # matrix [i,j] of sds for temperature offset, for ages.prox (i) and sites (j)
toff.sd = c(1:(length(ages.prox)*(length(INPUTtoff.sd[,1])-1)))
dim(toff.sd) = c((length(ages.prox)), (length(INPUTtoff.sd[,1])-1))
for (i in 2:length(INPUTtoff.sd)){
  toff.sd[,i-1] <- approx(x = INPUTtoff.sd$age, y = as.numeric(unlist(INPUTtoff.sd[i])), xout=ages.prox, method="linear")[["y"]]
}
############################################################################################


# Environmental prior (d13C of atm CO2) 
############################################################################################
d13CO2.m <- -6
d13CO2.sd <- 2
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass = list("d13Cpf.data" = d13Cpf.data,   
                 "d13Cpf.sd" = d13Cpf.sd,  
                 "ai.d13Cpf" = ai.d13Cpf,
                 "si.d13Cpf" = si.d13Cpf,
                 "d13Cbf.data" = d13Cbf.data,   
                 "d13Cbf.sd" = d13Cbf.sd,
                 "ai.d13Cbf" = ai.d13Cbf,
                 "si.d13Cbf" = si.d13Cbf,
                 "d13Cbulk.data" = d13Cbulk.data,   
                 "d13Cbulk.sd" = d13Cbulk.sd,  
                 "ai.d13Cbulk" = ai.d13Cbulk,
                 "si.d13Cbulk" = si.d13Cbulk,
                 "ai.all" = ai.all,
                 "n.steps" = n.steps,
                 "dt" = dt,
                 "n.sites" = n.sites,
                 "GMST.m" = GMST.m,
                 "GMST.sd" = GMST.sd,
                 "toff.m" = toff.m,
                 "toff.sd" = toff.sd,
                 "d13CO2.m" = d13CO2.m,
                 "d13CO2.sd" = d13CO2.sd) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("tempC", "d13CO2", "d13Cpf", "d13Cbf", "GMST", "toff")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "d13CO2_PSM.R", parameters.to.save = parms,
                        inits = NULL, n.chains = 3, n.iter = 1e4,
                        n.burnin = 5e3, n.thin = 1)

############################################################################################

  
  
  
  
  