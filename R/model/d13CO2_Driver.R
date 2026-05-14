

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

# select from 'PhanDA' or 'Scotese21' 
GMST_model <- "PhanDA" 

# uniform standard deviation applied to GMST from Scotese et al. (2021) used in Li et al (2022); 
# they do not provide uncertainty estimates
GMST_sd_Scotese21 <- 5

# select 'PhanDA' for Li22 MAT minus PhanDA GMST; 
# select 'Li22' for Li22 MAT minus Scotese21 GMST (same as Li22 GMST); 
temp_offset_model <- "Li22" 

# uniform standard deviation applied to all surface temperature site offset values - 
# no uncertainty in Scotese21 estimates
toff_sd_uniform <- 2 

# uniform standard deviation applied to all bottom temperature site offset values 
toff_sd_uniform_bot <- 1 


# residual offset in non-secular bias of archival d13C from high-fidelity (true) value (per mille), 
# standard deviation w/ fixed mean = 0 

# Benthic forams
bf.nsb.m <- 0
bf.nsb.sd <- 0.25

# Planktic forams
pf.nsb.m <- 0
pf.nsb.sd <- 0.25

# Brachiopods
brach.nsb.m <- 0
brach.nsb.sd <- 1

# Bivalves
bivalve.nsb.m <- 0
bivalve.nsb.sd <- 1

# Ammonites
amm.nsb.m <- 0
amm.nsb.sd <- 1

# Belemnites
bel.nsb.m <- 0
bel.nsb.sd <- 1

# micrite - open ocean
micrite.nsb.m <- 0
micrite.nsb.sd <- 0.25

# Bulk carbonate - open ocean
bulk.nsb.m <- 0
bulk.nsb.sd <- 0.5

# Bulk carbonate - semi-restricted
bulk_sr.nsb.m <- 0
bulk_sr.nsb.sd <- 1

# Bulk carbonate - marginal sea
bulk_marg.nsb.m <- 0
bulk_marg.nsb.sd <- 0.75

############################################################################################


# load and prepare proxy and climate data, indexing vectors and matrices 
############################################################################################

# load proxy data
prox.in <- as.data.frame(read.csv(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv"))
prox.in <- cbind(prox.in[,1:7], prox.in[,9:10], prox.in[,21:27],rep(x = toff_sd_uniform, times = nrow(prox.in)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category", 
                    "paleolon","paleolat", "MAT", "GMST_Scotese21", "GMST_PhanDA", "GMST_PhanDA_hi",
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
} else if (GMST_model == "Scotese21"){
  GMST.m <- approx(prox.in$age, prox.in$GMST_Scotese21, xout = ages, rule = 2)$y
  GMST.sd <- rep(x = GMST_sd_Scotese21, times = length(ages)) 
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
flattened$GMST_Scotese21_interp <- approx(prox.in$age, prox.in$GMST_Scotese21,    xout = flattened$ages, rule = 2)$y
flattened$GMST_PhanDA_sd_interp <- approx(prox.in$age, PhanDA_sd, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_interp <- approx(prox.in$age, prox.in$temp_offset, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_PhanDA_interp <- approx(prox.in$age, prox.in$temp_offset_PhanDA, xout = flattened$ages, rule = 2)$y
flattened$temp_offset_sd_interp <- approx(prox.in$age, prox.in$temp_offset_sd,  xout = flattened$ages, rule = 2)$y
flattened <- flattened[order(flattened$ai, flattened$site.index), ]
flattened$row.index <- 1:nrow(flattened)
rownames(flattened) <- NULL
si.flat <- flattened$site.index

# clean and prepare proxy data 
clean.d13C <- prox.in[complete.cases(prox.in$d13C), ]
clean.d13Cbf <- clean.d13C[clean.d13C$category == "bf",]
clean.d13Cpf <- clean.d13C[clean.d13C$category == "Planktonic foraminifera",]
clean.d13Cbrach <- clean.d13C[clean.d13C$category == "Brachiopod calcite",]
clean.d13Cbivalve <- clean.d13C[clean.d13C$category == "Bivalve",]
clean.d13Camm <- clean.d13C[clean.d13C$category == "Ammonite",]
clean.d13Cbel <- clean.d13C[clean.d13C$category == "Belemnite",]
clean.d13Cmicrite <- clean.d13C[clean.d13C$category == "micrite open ocean",]
clean.d13Cbulk <- clean.d13C[clean.d13C$category %in% c("bulk", "bulk open water", "bulk open ocean"), ]
clean.d13Cbulk_sr <- clean.d13C[clean.d13C$category == "bulk semi restricted",]
clean.d13Cbulk_marg <- clean.d13C[clean.d13C$category %in% c("bulk marginal sea", "bulk marginal sea restricting up section"), ]


# benthic foraminifera - determine index vectors and rank matrix
ai.d13Cbf <- sort(c(as.integer(clean.d13Cbf$ai)), decreasing = FALSE)    
si.d13Cbf <- clean.d13Cbf$site.index
d13Cbf.data <- clean.d13Cbf$d13C
n.d13Cbf = length(d13Cbf.data)

# planktic foraminifera - determine index vectors and rank matrix
ai.d13Cpf <- sort(c(as.integer(clean.d13Cpf$ai)), decreasing = FALSE)    
si.d13Cpf <- clean.d13Cpf$site.index
d13Cpf.data <- clean.d13Cpf$d13C
n.d13Cpf = length(d13Cpf.data)

# brachiopod - determine index vectors and rank matrix
ai.d13Cbrach <- sort(c(as.integer(clean.d13Cbrach$ai)), decreasing = FALSE)    
si.d13Cbrach <- clean.d13Cbrach$site.index
d13Cbrach.data <- clean.d13Cbrach$d13C
n.d13Cbrach = length(d13Cbrach.data)

# bivalve - determine index vectors and rank matrix
ai.d13Cbivalve <- sort(c(as.integer(clean.d13Cbivalve$ai)), decreasing = FALSE)    
si.d13Cbivalve <- clean.d13Cbivalve$site.index
d13Cbivalve.data <- clean.d13Cbivalve$d13C
n.d13Cbivalve = length(d13Cbivalve.data)

# ammonite - determine index vectors and rank matrix
ai.d13Camm <- sort(c(as.integer(clean.d13Camm$ai)), decreasing = FALSE)    
si.d13Camm <- clean.d13Camm$site.index
d13Camm.data <- clean.d13Camm$d13C
n.d13Camm = length(d13Camm.data)

# belemnite - determine index vectors and rank matrix
ai.d13Cbel <- sort(c(as.integer(clean.d13Cbel$ai)), decreasing = FALSE)    
si.d13Cbel <- clean.d13Cbel$site.index
d13Cbel.data <- clean.d13Cbel$d13C
n.d13Cbel = length(d13Cbel.data)

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

ri.d13Cpf <- flattened$row.index[
  ri.d13Cpf <- match(
    interaction(clean.d13Cpf$ai, clean.d13Cpf$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbrach <- flattened$row.index[
  ri.d13Cbrach <- match(
    interaction(clean.d13Cbrach$ai, clean.d13Cbrach$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbivalve <- flattened$row.index[
  ri.d13Cbivalve <- match(
    interaction(clean.d13Cbivalve$ai, clean.d13Cbivalve$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Camm <- flattened$row.index[
  ri.d13Camm <- match(
    interaction(clean.d13Camm$ai, clean.d13Camm$site.index),
    interaction(flattened$ai, flattened$site.index))
]

ri.d13Cbel <- flattened$row.index[
  ri.d13Cbel <- match(
    interaction(clean.d13Cbel$ai, clean.d13Cbel$site.index),
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
stopifnot(all(ri.d13Cpf >= 1 & ri.d13Cpf <= nrow(flattened)))
stopifnot(all(ri.d13Cbrach >= 1 & ri.d13Cbrach <= nrow(flattened)))
stopifnot(all(ri.d13Cbivalve >= 1 & ri.d13Cbivalve <= nrow(flattened)))
stopifnot(all(ri.d13Camm >= 1 & ri.d13Camm <= nrow(flattened)))
stopifnot(all(ri.d13Cbel >= 1 & ri.d13Cbel <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk >= 1 & ri.d13Cbulk <= nrow(flattened)))
stopifnot(all(ri.d13Cmicrite >= 1 & ri.d13Cmicrite <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk_sr >= 1 & ri.d13Cbulk_sr <= nrow(flattened)))
stopifnot(all(ri.d13Cbulk_marg >= 1 & ri.d13Cbulk_marg <= nrow(flattened)))

# No silent duplication problems?
stopifnot(length(unique(flattened$row.index)) == nrow(flattened))

# BWT interpolated for time steps
BWT.Cen <- as.data.frame(read.csv(file = "data/raw/CenozoicBWT.csv"))
names(BWT.Cen) <- c("age", "BWT", "BWT_2sd")
BWT.Cen$age <- BWT.Cen$age*1e3
BWT.Cen <- BWT.Cen[order(BWT.Cen$age, decreasing = TRUE),]
BWT.Cen_last <- cbind(age.max.spinup, BWT.Cen[1,2:3])
names(BWT.Cen_last) <- c("age", "BWT", "BWT_2sd")
BWT <- rbind(BWT.Cen, BWT.Cen_last)
BWT <- BWT[order(BWT$age, decreasing = TRUE),]
BWT.m <- approx(BWT$age, BWT$BWT, xout=ages, method="linear", rule =2)
BWT.sd <- approx(BWT$age, BWT$BWT_2sd/2, xout=ages, method="linear", rule =2)
BWT.m <- BWT.m[["y"]]
BWT.sd <- BWT.sd[["y"]]


# Environmental prior (d13C of atm CO2) 
############################################################################################

d13CO2.l <- -12
d13CO2.u <- 0
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
points(flattened$ai[ri.d13Cbf], flattened$site.index[ri.d13Cbf], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cpf], flattened$site.index[ri.d13Cpf], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbrach], flattened$site.index[ri.d13Cbrach], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbivalve], flattened$site.index[ri.d13Cbivalve], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Camm], flattened$site.index[ri.d13Camm], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbel], flattened$site.index[ri.d13Cbel], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbf], flattened$site.index[ri.d13Cbf], col=2, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk], flattened$site.index[ri.d13Cbulk], col=4, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk_sr], flattened$site.index[ri.d13Cbulk_sr], col=6, pch=19, cex=.4)
points(flattened$ai[ri.d13Cbulk_marg], flattened$site.index[ri.d13Cbulk_marg], col=8, pch=19, cex=.4)
points(flattened$ai[ri.d13Cmicrite], flattened$site.index[ri.d13Cmicrite], col=10, pch=19, cex=.4)
legend('topright', c('bf','pf', 'brach', 'bivavle', 'ammonite', 'belemnite', ' bulk', 'micrite', 'sr', 'marg'), col=c(2,2,2,2,2,2,2,4,6,8,10), pch=19, bty='n')

# Select objects to pass to jags 
############################################################################################

data.pass = list("n.steps" = n.steps,
                 "dt" = dt,
                 "n.sites" = n.sites,
                 "si.flat" = si.flat,
                 "ai.flat" = ai.flat,
                 "GMST.obs" = GMST.m,
                 "GMST.sd" = GMST.sd,
                 "BWT.obs" = BWT.m, 
                 "BWT.sd" = BWT.sd, 
                 "toff_sd_uniform_bot" = toff_sd_uniform_bot,
                 "toff.m" = toff.m,
                 "toff.sd" = toff.sd,
                 "d13CO2.l" = d13CO2.l,
                 "d13CO2.u" = d13CO2.u) 

data.pass.bf = list("d13Cbf.data" = d13Cbf.data,   
                    "ai.d13Cbf" = ai.d13Cbf,
                    "ri.d13Cbf" = ri.d13Cbf,
                    "n.d13Cbf" = n.d13Cbf,
                    "bf.nsb.m" = bf.nsb.m,
                    "bf.nsb.sd" = bf.nsb.sd)

data.pass.pf = list("d13Cpf.data" = d13Cpf.data,   
                    "ai.d13Cpf" = ai.d13Cpf,
                    "ri.d13Cpf" = ri.d13Cpf,
                    "n.d13Cpf" = n.d13Cpf,
                    "pf.nsb.m" = pf.nsb.m,
                    "pf.nsb.sd" = pf.nsb.sd)

data.pass.brach = list("d13Cbrach.data" = d13Cbrach.data,   
                    "ai.d13Cbrach" = ai.d13Cbrach,
                    "ri.d13Cbrach" = ri.d13Cbrach,
                    "n.d13Cbrach" = n.d13Cbrach,
                    "brach.nsb.m" = brach.nsb.m,
                    "brach.nsb.sd" = brach.nsb.sd)

data.pass.bivalve = list("d13Cbivalve.data" = d13Cbivalve.data,   
                    "ai.d13Cbivalve" = ai.d13Cbivalve,
                    "ri.d13Cbivalve" = ri.d13Cbivalve,
                    "n.d13Cbivalve" = n.d13Cbivalve,
                    "bivalve.nsb.m" = bivalve.nsb.m,
                    "bivalve.nsb.sd" = bivalve.nsb.sd)

data.pass.amm = list("d13Camm.data" = d13Camm.data,   
                    "ai.d13Camm" = ai.d13Camm,
                    "ri.d13Camm" = ri.d13Camm,
                    "n.d13Camm" = n.d13Camm,
                    "amm.nsb.m" = amm.nsb.m,
                    "amm.nsb.sd" = amm.nsb.sd)

data.pass.bel = list("d13Cbel.data" = d13Cbel.data,   
                    "ai.d13Cbel" = ai.d13Cbel,
                    "ri.d13Cbel" = ri.d13Cbel,
                    "n.d13Cbel" = n.d13Cbel,
                    "bel.nsb.m" = bel.nsb.m,
                    "bel.nsb.sd" = bel.nsb.sd)

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


data.pass <- c(data.pass, data.pass.bf, data.pass.pf, data.pass.brach, data.pass.bivalve, data.pass.amm, 
               data.pass.bel, data.pass.micrite, data.pass.bulk, data.pass.bulk_sr, data.pass.bulk_marg)

############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("d13CO2", "GMST", "BWT", "BWT_GMST_beta", "tempC", "tempC_bot", "toff", "toff_bot", "d13Cbf", "d13Cpf", 
          "d13Cbrach", "d13Cbivalve", "d13Camm", "d13Cbel", "d13Cbulk", "d13Cbulk_sr", "d13Cbulk_marg", 
          "d13Cmicrite", "bf.nsb_site", "pf.nsb_site", "brach.nsb_site", "bivalve.nsb_site", "amm.nsb_site", 
          "bel.nsb_site", "bulk.nsb_site",  "micrite.nsb_site", "bulk_sr.nsb_site", "bulk_marg.nsb_site", "Abf", 
          "Asurf")

############################################################################################


# Run the inversion using jags 
############################################################################################
system.time({inv.out = jags.parallel(data = data.pass, model.file = "R/model/d13CO2_PSM.R", 
                                     parameters.to.save = parms, inits = NULL, n.chains = 6, 
                                     n.iter = 1e5, n.burnin = 3e4, n.thin = 100)})


############################################################################################
# 550 time steps, 3 chains and 1e4 iteration takes 6.5 minutes 
# 550 time steps, 6 chains and 5e4 iteration takes 38 minutes 
# 550 time steps, 6 chains and 3e5 iteration takes 3.5 hours

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




