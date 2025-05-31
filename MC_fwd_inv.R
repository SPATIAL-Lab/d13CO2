
# monte carlo style forward model for d13Cbulk from d13CO2 and temp


# ----- monte carlo style inverse model for d13bulk from d13CO2 and temp -----

n.iterations <- 10000
tempC.m <- 32
tempC.sd <- 3
d13CO2.m <- -9
d13CO2.sd <- 1

init.vec <- vector("numeric", n.iterations)
coeff1 <- init.vec
coeff2 <- init.vec
eps.int <- init.vec
f_co2 <- init.vec
tempC <- init.vec
d13CO2 <- init.vec
eps.dic_co2_surf <- init.vec
d13Cbulk <- init.vec

for (i in 1:n.iterations){
  
  coeff1[i] <- rnorm(1, mean = 0.0144, sd = 0.01)
  coeff2[i] <- rnorm(1, mean = 0.107, sd = 0.002)
  eps.int[i] <- rnorm(1, mean = 10.53, sd = 0.05)
  f_co2[i] <- rnorm(1, mean = 0.12, sd = 0.04)
  
  tempC[i] <- rnorm(1, mean = tempC.m, sd = tempC.sd)
  d13CO2[i] <- rnorm(1, mean = d13CO2.m, sd = d13CO2.sd)
  
  eps.dic_co2_surf[i] <- coeff1[i]*tempC[i]*f_co2[i] - coeff2[i]*tempC[i] + eps.int[i]
  d13Cbulk[i] <- d13CO2[i] + eps.dic_co2_surf[i]
}

print(mean(d13Cbulk))
print(sd(d13Cbulk))

plot(x = c(-5,6), y =c(0,0.6), col= NULL, xlab = expression(paste(delta^{13}, "C bulk (\u2030 VPDB)")),
     ylab = "Probability Density")
polygon(density(d13Cbulk),
        col = adjustcolor("coral1", alpha.f = 0.80), border = NA, title = NULL)




# ----- monte carlo style inverse model for d13CO2 from d13Cbulk and temp -----

n.iterations <- 10000
tempC.m <- 32
tempC.sd <- 3
d13Cbulk.m <- -1
d13Cbulk.sd <- 0.5

init.vec <- vector("numeric", n.iterations)
coeff1 <- init.vec
coeff2 <- init.vec
eps.int <- init.vec
f_co2 <- init.vec
tempC <- init.vec
d13CO2 <- init.vec
eps.dic_co2_surf <- init.vec
d13Cbulk <- init.vec

for (i in 1:n.iterations){
  
  coeff1[i] <- rnorm(1, mean = 0.0144, sd = 0.01)
  coeff2[i] <- rnorm(1, mean = 0.107, sd = 0.002)
  eps.int[i] <- rnorm(1, mean = 10.53, sd = 0.05)
  f_co2[i] <- rnorm(1, mean = 0.12, sd = 0.04)
  
  tempC[i] <- rnorm(1, mean = tempC.m, sd = tempC.sd)
  d13Cbulk[i] <- rnorm(1, mean = d13Cbulk.m, sd = d13bulk.sd)
  
  eps.dic_co2_surf[i] <- coeff1[i]*tempC[i]*f_co2[i] - coeff2[i]*tempC[i] + eps.int[i]
  d13CO2[i] <- d13Cbulk[i] - eps.dic_co2_surf[i] 
}

print(mean(d13Cbulk))
print(sd(d13Cbulk))

plot(x = c(-10,-2), y =c(0,0.6), col= NULL, xlab = expression(paste(delta^{13}, 
    "C of ",CO[2]," (\u2030 VPDB)")), ylab = "Probability Density")
polygon(density(d13CO2),
        col = adjustcolor("dodgerblue", alpha.f = 0.80), border = NA, title = NULL)


