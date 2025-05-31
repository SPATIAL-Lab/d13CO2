

# Driver script for PSM used to reconstruct d13C of atmospheric CO2 
# Dustin T. Harper


# libraries
############################################################################################

library(R2jags)


# set age range, time step interval (resolution) for reconstruction, and number of spinup time steps
############################################################################################

age.min <- 0
age.max <- 540000
step.int <- 1000
n.spinup <- 10

# select from 'PhanDA' or 'Li22' 
GMST_model <- "PhanDA" 

# uniform standard deviation applied to GMST from Li et al (2022); they do not provide uncertainty estimates
GMST_sd_Li22 <- 5

# select 'PhanDA' for Li22 MAT minus PhanDA GMST; select 'Li22' for Li22 MAT minus Li22 GMST; 
temp_offset_model <- "Li22" 

# uniform standard deviation applied to all surface temperature site offset values - no uncertainty in Li22 estimates
toff_sd_uniform <- 2 

# uniform standard deviation applied to all bottom temperature site offset values 
toff_sd_uniform_bot <- 1 


# ----- Alteration -----
# alteration of archive d13C from 'pristine' value (per mille), standard deviation w/ fixed mean = 0 

# Benthic forams
bf.alteration <- 0.5

# Bulk carbonate
bulk.alteration <- 0.5

# Brachiopods
brach.alteration <- 0.5
############################################################################################


# load and prepare proxy and climate data, indexing vectors and matrices 
############################################################################################

# load proxy data
prox.in <- as.data.frame(read.csv(file = "PhanData/PhanCompWithTemp.csv"))
prox.in <- cbind(prox.in[,2:10], prox.in[,12:18], rep(x = toff_sd_uniform, times = nrow(prox.in)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "material", 
                    "paleolon","paleolat", "MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                    "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA", "temp_offset_sd")

# age index proxy data (in kyrs)
prox.in$age <- prox.in$age*1e3
prox.in <- prox.in[prox.in$age >= (age.min) & prox.in$age <= (age.max),]
prox.in <- transform(prox.in, ai = n.spinup + as.numeric(1 + floor((max(prox.in$age) - prox.in$age) / step.int)))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE),]

ages.short <- seq(from = max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
ages <- seq(from = n.spinup*step.int + max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
ai.all <- c(c(1:n.spinup), sort(unique(prox.in$ai), decreasing = FALSE))
age.indices <- as.data.frame(cbind(ai.all, ages))
names(age.indices) <- c("ai", "age")

dt <- abs(diff(unique(ages), lag=1))
n.steps <- as.numeric(length(dt)+1)
age.max.spinup <- age.max + step.int*n.steps 


# calculate symmetric uncertainty (GMST_PhanDA_sd)
PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                (prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo)) / 2


# interpolate GMST for ages 
if (GMST_model == "PhanDA"){
  GMST.m <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = ages, rule = 2)$y
  GMST.sd <- approx(prox.in$age, PhanDA_sd, xout = ages, rule = 2)$y
} else if (GMST_model == "Li22"){
  GMST.m <- approx(prox.in$age, prox.in$GMST_Li22, xout = ages, rule = 2)$y
  GMST.sd <- rep(x = GMST_sd_Li22, times = length(ages)) 
}


# site index proxy data
prox.in <- transform(prox.in,site.index=as.numeric(factor(site, ordered = is.ordered(site))))
site.index <- c(prox.in$site.index)
n.sites <- as.numeric(length(unique(site.index)))

# list the input sites and their associated site index
sites <- data.frame((sort(unique(prox.in$site), decreasing = FALSE)),seq(1:length(sort(unique(prox.in$site), decreasing = FALSE))))
names(sites) <- c("site", "site.index")
print(sites)
print("Make sure the site indexes in the temp offset .csv files are ordered as they are here! (Numerical increasing)")


# generate a 'flattened' data frame which has all interpolated data for existing combinations of 'ai' and 'site.index'
flattened <- unique(prox.in[c("ai", "site.index")])
flattened <- flattened[order(flattened$ai, flattened$site.index), ]
ai.flat <- flattened$ai

# add 'ages' column using age.indices (assumes age.indices has 'ai' and 'age')
flattened$ages <- age.indices$age[match(flattened$ai, age.indices$ai)]

# interpolate each column at the specified 'ages'
flattened$GMST_PhanDA_interp <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = flattened$ages, rule = 2)$y
flattened$GMST_Li22_interp <- approx(prox.in$age, prox.in$GMST_Li22,    xout = flattened$ages, rule = 2)$y
flattened$GMST_PhanDA_sd_interp <- approx(prox.in$age, PhanDA_sd, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_interp <- approx(prox.in$age, prox.in$temp_offset, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_PhanDA_interp <- approx(prox.in$age, prox.in$temp_offset_PhanDA, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_sd_interp <- approx(prox.in$age, prox.in$temp_offset_sd,  xout = flattened$ages, rule = 2)$y
flattened <- flattened[order(flattened$ai, flattened$site.index), ]
flattened$row.index <- 1:nrow(flattened)
rownames(flattened) <- NULL


# clean and prepare proxy data 
clean.d13C <- prox.in[complete.cases(prox.in$d13C), ]
clean.d13Cbf <- clean.d13C[clean.d13C$material == "bf",]
clean.d13Cbulk <- clean.d13C[clean.d13C$material == "bulk",]
clean.d13Cbrach <- clean.d13C[clean.d13C$material == "brach",]

# benthic foraminifera - determine index vectors and rank matrix
ai.d13Cbf <- sort(c(as.integer(clean.d13Cbf$ai)), decreasing = FALSE)    
si.d13Cbf <- clean.d13Cbf$site.index
d13Cbf.data <- clean.d13Cbf$d13C
n.d13Cbf = length(d13Cbf.data)


# bulk carbonate - determine index vectors and rank matrix
ai.d13Cbulk <- sort(c(as.integer(clean.d13Cbulk$ai)), decreasing = FALSE)    
si.d13Cbulk <- clean.d13Cbulk$site.index
d13Cbulk.data <- clean.d13Cbulk$d13C
n.d13Cbulk = length(d13Cbulk.data)


# brachiopod - determine index vectors and rank matrix
ai.d13Cbrach <- sort(c(as.integer(clean.d13Cbrach$ai)), decreasing = FALSE)    
si.d13Cbrach <- clean.d13Cbrach$site.index
d13Cbrach.data <- clean.d13Cbrach$d13C
n.d13Cbrach = length(d13Cbrach.data)

# index each row of data to 'flattened' data frame combinations 
ri.d13Cbf <- flattened$row.index[
ri.d13Cbf <- match(
  interaction(clean.d13Cbf$ai, clean.d13Cbf$site.index),
  interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbulk <- flattened$row.index[
ri.d13Cbulk <- match(
  interaction(clean.d13Cbulk$ai, clean.d13Cbulk$site.index),
  interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbrach <- flattened$row.index[
ri.d13Cbrach <- match(
  interaction(clean.d13Cbrach$ai, clean.d13Cbrach$site.index),
  interaction(flattened$ai, flattened$site.index))
]

# BWT interpolated for time steps
BWT.Cen <- as.data.frame(read.csv(file = "PhanData/CenozoicBWT.csv"))
names(BWT.Cen) <- c("age", "BWT", "BWT_2sd")
BWT.Cen$age <- BWT.Cen$age*1e3
BWT.Cen <- BWT.Cen[order(BWT.Cen$age, decreasing = TRUE),]
BWT.Cen_last <- cbind(age.max.spinup, BWT.Cen[1,2:3])
names(BWT.Cen_last) <- c("age", "BWT", "BWT_2sd")
BWT <- rbind(BWT.Cen, BWT.Cen_last)
BWT <- BWT[order(BWT$age, decreasing = TRUE),]
BWT.m <- approx(BWT$age, BWT$BWT, xout=ages, method="linear")
BWT.sd <- approx(BWT$age, BWT$BWT_2sd/2, xout=ages, method="linear")
BWT.m <- BWT.m[["y"]]
BWT.sd <- BWT.sd[["y"]]


# Environmental prior (d13C of atm CO2) 
############################################################################################

d13CO2.m <- -6
d13CO2.sd <- 3
############################################################################################


if (temp_offset_model == "Li22"){
toff.m <- flattened$temp_offset_interp
toff.sd <- flattened$temp_offset_sd_interp
} else if (temp_offset_model == "PhanDA"){
  toff.m <- flattened$temp_offset_PhanDA_interp
  toff.sd <- flattened$temp_offset_sd_interp
 }

# Select objects to pass to jags 
############################################################################################

data.pass = list("n.steps" = n.steps,
                 "dt" = dt,
                 "n.sites" = n.sites,
                 "ai.flat" = ai.flat,
                 "GMST.m" = GMST.m,
                 "GMST.sd" = GMST.sd,
                 "BWT.m" = BWT.m, 
                 "BWT.sd" = BWT.sd, 
                 "toff_sd_uniform_bot" = toff_sd_uniform_bot,
                 "toff.m" = toff.m,
                 "toff.sd" = toff.sd,
                 "d13CO2.m" = d13CO2.m,
                 "d13CO2.sd" = d13CO2.sd) 

data.pass.bf = list("d13Cbf.data" = d13Cbf.data,   
                    "ai.d13Cbf" = ai.d13Cbf,
                    "ri.d13Cbf" = ri.d13Cbf,
                    "n.d13Cbf" = n.d13Cbf,
                    "bf.alteration" =bf.alteration)

data.pass.bulk = list("d13Cbulk.data" = d13Cbulk.data,   
                      "ai.d13Cbulk" = ai.d13Cbulk,
                      "ri.d13Cbulk" = ri.d13Cbulk,
                      "n.d13Cbulk" = n.d13Cbulk,
                      "bulk.alteration" = bulk.alteration)

data.pass.brach = list("d13Cbrach.data" = d13Cbrach.data,   
                       "ai.d13Cbrach" = ai.d13Cbrach,
                       "ri.d13Cbrach" = ri.d13Cbrach,
                       "n.d13Cbrach" = n.d13Cbrach,
                       "brach.alteration" = brach.alteration)

data.pass <- c(data.pass, data.pass.bf, data.pass.bulk, data.pass.brach)

############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("d13CO2", "GMST", "BWT", "tempC", "tempC_bot", "toff", "toff_bot", "d13Cbf", "d13Cbulk")
############################################################################################


# Run the inversion using jags 
############################################################################################
n.chains <- 6
n.iter <- 5e4
n.burnin <- 1e4
n.thin <- 10

system.time({inv.out = jags.parallel(data = data.pass, model.file = "Phan/d13CO2_PSM.R", 
                                     parameters.to.save = parms, inits = NULL, n.chains = 3, 
                                     n.iter = 1e4, n.burnin = 5e3, n.thin = 1)})

############################################################################################
# 110 time steps, 3 chains and 1e4 iteration takes 4 minutes 
# 110 time steps, 6 chains and 5e4 iteration takes 24 minutes 
# 550 time steps, 6 chains and 5e4 iteration takes 54 minutes 


