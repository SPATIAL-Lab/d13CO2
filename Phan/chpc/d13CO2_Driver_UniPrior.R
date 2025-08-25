

# Driver script for PSM used to reconstruct d13C of Phanerozoic atmospheric CO2 
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
GMST_model <- "Li22" 

# uniform standard deviation applied to GMST from Li et al (2022); they do not provide uncertainty estimates
GMST_sd_Li22 <- 5

# select 'PhanDA' for Li22 MAT minus PhanDA GMST; select 'Li22' for Li22 MAT minus Li22 GMST; 
temp_offset_model <- "Li22" 

# uniform standard deviation applied to all surface temperature site offset values - no uncertainty in Li22 estimates
toff_sd_uniform <- 2 

# uniform standard deviation applied to all bottom temperature site offset values 
toff_sd_uniform_bot <- 1 


# non-secular bias of archival d13C from high-fidelity (true) value (per mille), standard deviation w/ fixed mean = 0 

# Benthic forams
bf.nsb.m <- 0
bf.nsb.sd <- 0.5

# Bulk carbonate - open ocean
bulk.nsb.m <- 0
bulk.nsb.sd <- 0.5

# micrite - open ocean
micrite.nsb.m <- 0
micrite.nsb.sd <- 0.25

# Bulk carbonate - semi-restricted
bulk_sr.nsb.m <- 1
bulk_sr.nsb.sd <- 1

# Bulk carbonate - marginal sea
bulk_marg.nsb.m <- 0.5
bulk_marg.nsb.sd <- 1

############################################################################################


# load and prepare proxy and climate data, indexing vectors and matrices 
############################################################################################

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
ai.all <- c(c(1:n.spinup), sort(unique(prox.in$ai), decreasing = FALSE))
age.indices <- as.data.frame(cbind(ai.all, ages))
names(age.indices) <- c("ai", "age")

dt <- abs(diff(unique(ages), lag=1))
n.steps <- as.numeric(length(dt)+1)
age.max.spinup <- age.max + step.int*n.spinup 


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
clean.d13Cbf <- clean.d13C[clean.d13C$category == "bf",]
clean.d13Cmicrite <- clean.d13C[clean.d13C$category == "micrite open ocean",]
clean.d13Cbulk <- clean.d13C[clean.d13C$category %in% c("bulk", "bulk open water", "bulk open ocean"), ]
clean.d13Cbulk_sr <- clean.d13C[clean.d13C$category == "bulk semi restricted",]
clean.d13Cbulk_marg <- clean.d13C[clean.d13C$category %in% c("bulk marginal sea", "bulk marginal sea restricting up section"), ]

# benthic foraminifera - determine index vectors and rank matrix
ai.d13Cbf <- sort(c(as.integer(clean.d13Cbf$ai)), decreasing = FALSE)    
si.d13Cbf <- clean.d13Cbf$site.index
d13Cbf.data <- clean.d13Cbf$d13C
n.d13Cbf = length(d13Cbf.data)

# micrite - determine index vectors and rank matrix
ai.d13Cmicrite <- sort(c(as.integer(clean.d13Cmicrite$ai)), decreasing = FALSE)
si.d13Cmicrite <- clean.d13Cmicrite$site.index
d13Cmicrite.data <- clean.d13Cmicrite$d13C
n.d13Cmicrite = length(d13Cmicrite.data)

# bulk carbonate - determine index vectors and rank matrix
ai.d13Cbulk <- sort(c(as.integer(clean.d13Cbulk$ai)), decreasing = FALSE)    
si.d13Cbulk <- clean.d13Cbulk$site.index
d13Cbulk.data <- clean.d13Cbulk$d13C
n.d13Cbulk = length(d13Cbulk.data)

# bulk carbonate semi restricted - determine index vectors and rank matrix
ai.d13Cbulk_sr <- sort(c(as.integer(clean.d13Cbulk_sr$ai)), decreasing = FALSE)
si.d13Cbulk_sr <- clean.d13Cbulk_sr$site.index
d13Cbulk_sr.data <- clean.d13Cbulk_sr$d13C
n.d13Cbulk_sr = length(d13Cbulk_sr.data)

# bulk carbonate marginal sea - determine index vectors and rank matrix
ai.d13Cbulk_marg <- sort(c(as.integer(clean.d13Cbulk_marg$ai)), decreasing = FALSE)
si.d13Cbulk_marg <- clean.d13Cbulk_marg$site.index
d13Cbulk_marg.data <- clean.d13Cbulk_marg$d13C
n.d13Cbulk_marg = length(d13Cbulk_marg.data)

# index each row of data to 'flattened' data frame combinations 
ri.d13Cbf <- flattened$row.index[
ri.d13Cbf <- match(
  interaction(clean.d13Cbf$ai, clean.d13Cbf$site.index),
  interaction(flattened$ai, flattened$site.index))
]

ri.d13Cmicrite <- flattened$row.index[
  ri.d13Cmicrite <- match(
    interaction(clean.d13Cmicrite$ai, clean.d13Cmicrite$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbulk <- flattened$row.index[
  ri.d13Cbulk <- match(
    interaction(clean.d13Cbulk$ai, clean.d13Cbulk$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbulk_sr <- flattened$row.index[
ri.d13Cbulk_sr <- match(
  interaction(clean.d13Cbulk_sr$ai, clean.d13Cbulk_sr$site.index),
  interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbulk_marg <- flattened$row.index[
  ri.d13Cbulk_marg <- match(
    interaction(clean.d13Cbulk_marg$ai, clean.d13Cbulk_marg$site.index),
    interaction(flattened$ai, flattened$site.index))
]

# All mapped?
stopifnot(!any(is.na(ri.d13Cbf)), !any(is.na(ri.d13Cbulk)))

# Indices in range?
stopifnot(all(ri.d13Cbf >= 1 & ri.d13Cbf <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk >= 1 & ri.d13Cbulk <= nrow(flattened)))
stopifnot(all(ri.d13Cmicrite >= 1 & ri.d13Cmicrite <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk_sr >= 1 & ri.d13Cbulk_sr <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk_marg >= 1 & ri.d13Cbulk_marg <= nrow(flattened)))

# No silent duplication problems?
stopifnot(length(unique(flattened$row.index)) == nrow(flattened))

# BWT interpolated for time steps
BWT.Cen <- as.data.frame(read.csv(file = "Phan/PhanData/CenozoicBWT.csv"))
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

d13CO2.u <- 0
d13CO2.l <- -12
############################################################################################


if (temp_offset_model == "Li22"){
toff.m <- flattened$temp_offset_interp
toff.sd <- flattened$temp_offset_sd_interp
} else if (temp_offset_model == "PhanDA"){
  toff.m <- flattened$temp_offset_PhanDA_interp
  toff.sd <- flattened$temp_offset_sd_interp
}

stopifnot(!any(is.na(GMST.m)), !any(is.na(GMST.sd)))
stopifnot(!any(is.na(BWT.m)), !any(is.na(BWT.sd)))


# Catch indexing misalignment

plot(age.indices$age, GMST.m, type='l', main='GMST.m on JAGS grid', xlab='age (kyr)', ylab='°C')
lines(age.indices$age, GMST.sd, lty=2)

plot(age.indices$age, BWT.m,  type='l', main='BWT.m on JAGS grid', xlab='age (kyr)', ylab='°C')
lines(age.indices$age, BWT.sd, lty=2)

# Where are observations in (ai, site) space?
plot(flattened$ai, flattened$site.index, pch='.', main='Flattened (ai, site)')
points(flattened$ai[ri.d13Cbf],   flattened$site.index[ri.d13Cbf],   col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk], flattened$site.index[ri.d13Cbulk], col=4, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk_sr], flattened$site.index[ri.d13Cbulk_sr], col=6, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk_marg], flattened$site.index[ri.d13Cbulk_marg], col=8, pch=19, cex=.4)
points(flattened$ai[ri.d13Cmicrite], flattened$site.index[ri.d13Cmicrite], col=10, pch=19, cex=.4)
legend('topright', c('bf','bulk', 'micrite', 'sr', 'marg'), col=c(2,4,6,8,10), pch=19, bty='n')

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
                 "d13CO2.u" = d13CO2.u,
                 "d13CO2.l" = d13CO2.l) 

data.pass.bf = list("d13Cbf.data" = d13Cbf.data,   
                    "ai.d13Cbf" = ai.d13Cbf,
                    "ri.d13Cbf" = ri.d13Cbf,
                    "n.d13Cbf" = n.d13Cbf,
                    "bf.nsb.m" = bf.nsb.m,
                    "bf.nsb.sd" = bf.nsb.sd)

data.pass.micrite = list("d13Cmicrite.data" = d13Cmicrite.data,
                      "ai.d13Cmicrite" = ai.d13Cmicrite,
                      "ri.d13Cmicrite" = ri.d13Cmicrite,
                      "n.d13Cmicrite" = n.d13Cmicrite,
                      "micrite.nsb.m" = micrite.nsb.m,
                      "micrite.nsb.sd" = micrite.nsb.sd)

data.pass.bulk = list("d13Cbulk.data" = d13Cbulk.data,   
                      "ai.d13Cbulk" = ai.d13Cbulk,
                      "ri.d13Cbulk" = ri.d13Cbulk,
                      "n.d13Cbulk" = n.d13Cbulk,
                      "bulk.nsb.m" = bulk.nsb.m,
                      "bulk.nsb.sd" = bulk.nsb.sd)

data.pass.bulk_sr = list("d13Cbulk_sr.data" = d13Cbulk_sr.data,
                      "ai.d13Cbulk_sr" = ai.d13Cbulk_sr,
                      "ri.d13Cbulk_sr" = ri.d13Cbulk_sr,
                      "n.d13Cbulk_sr" = n.d13Cbulk_sr,
                      "bulk_sr.nsb.m" = bulk_sr.nsb.m,
                      "bulk_sr.nsb.sd" = bulk_sr.nsb.sd)

data.pass.bulk_marg = list("d13Cbulk_marg.data" = d13Cbulk_marg.data,
                      "ai.d13Cbulk_marg" = ai.d13Cbulk_marg,
                      "ri.d13Cbulk_marg" = ri.d13Cbulk_marg,
                      "n.d13Cbulk_marg" = n.d13Cbulk_marg,
                      "bulk_marg.nsb.m" = bulk_marg.nsb.m,
                      "bulk_marg.nsb.sd" = bulk_marg.nsb.sd)


data.pass <- c(data.pass, data.pass.bf, data.pass.micrite, data.pass.bulk, data.pass.bulk_sr, data.pass.bulk_marg)

############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("d13CO2", "GMST", "BWT", "tempC", "tempC_bot", "toff", "toff_bot", "d13Cbf", "d13Cbulk",
        "d13Cbulk_sr", "d13Cbulk_marg", "d13Cmicrite")

############################################################################################


# Run the inversion using jags 
############################################################################################

system.time({inv.out = jags.parallel(data = data.pass, model.file = "Phan/d13CO2_PSM_UniPrior.R", 
                                     parameters.to.save = parms, inits = NULL, n.chains = 6, 
                                     n.iter = 3e5, n.burnin = 1e5, n.thin = 100)})

############################################################################################
# 550 time steps, 3 chains and 1e4 iteration takes 6.5 minutes 
# 550 time steps, 6 chains and 5e4 iteration takes 38 minutes 

save(inv.out, file = "inv.out_UniPrior.rda")
save(ages, file = "ages_UniPrior.rda")
save(prox.in, file = "prox.in_UniPrior.rda")
save(flattened, file = "flattened_UniPrior.rda")
save(BWT, file = "BWT_UniPrior.rda")
save(sites, file = "sites_UniPrior.rda")

# Posterior draws for model-predicted bulk d13C at every (ai, site) in `flattened`
bulk_fit <- inv.out$BUGSoutput$sims.list$d13Cbulk  # matrix: [iterations, nrow(flattened)]
stopifnot(ncol(bulk_fit) == nrow(flattened))

# Posterior mean at each flattened row
bulk_mu_all <- colMeans(bulk_fit)

# Match to the observation locations only
bulk_mu_obs <- bulk_mu_all[ri.d13Cbulk]
bulk_obs <- d13Cbulk.data

# Sanity check lengths
stopifnot(length(bulk_mu_obs) == length(bulk_obs))

# Plot observed vs posterior mean prediction (1:1 if well fit)
plot(bulk_obs, bulk_mu_obs, pch = 19, cex = 0.4,
     xlab = "Observed bulk d13C",
     ylab = "Posterior mean prediction")
abline(0, 1, lty = 2)




