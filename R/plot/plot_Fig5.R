

####################################################################################################
####################################################################################################
# Script to generate Figure 5 in manuscript: multipanel d13CO2, d13Corg, Delta, correlation, CO2, 
# and timescale plot

# d13CO2 posterior draws, plant d13Corg, d13Corg - d13CO2, windowed correlations,
# Foster et al. CO2, GSA period strip, and shared bottom x-axis 

# Assumes objects already exist in the environment:
# inv.out
# ages

load("output/model_runs/chpc_260511/output/inv.out_main.rda")
load("output/model_runs/chpc_260511/output/ages_main.rda")

####################################################################################################
####################################################################################################

library(readxl)

####################################################################################################
# d13CO2 from inv.out
####################################################################################################
d13CO2_draws <- as.matrix(inv.out$BUGSoutput$sims.list$d13CO2)

if (ncol(d13CO2_draws) == length(ages)) {
  d13CO2_draws <- t(d13CO2_draws)
}

d13CO2_med <- apply(d13CO2_draws, 1, median)
d13CO2_q025 <- apply(d13CO2_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
d13CO2_q975 <- apply(d13CO2_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
d13CO2_sd <- apply(d13CO2_draws, 1, sd, na.rm = TRUE)

co2_hw95 <- (d13CO2_q975 - d13CO2_q025) / 2
co2_sigma <- co2_hw95 / 1.96

n_draws <- min(100, ncol(d13CO2_draws))
ix_draws <- round(seq(1, ncol(d13CO2_draws), length.out = n_draws))

ages_Ma <- ages / 1000


####################################################################################################
# Plant d13Corg data
####################################################################################################
iso23 <- read_excel("data/raw/ISOORG23.xlsx", sheet = 1)
iso16 <- read_excel("data/raw/ISOORG16.xlsx", sheet = 1)

d23 <- data.frame(
  age = suppressWarnings(as.numeric(iso23[["1 Myr bin"]])),
  d13C = suppressWarnings(as.numeric(iso23[["d13Cp"]]))
)

d23 <- d23[is.finite(d23$age) &
             is.finite(d23$d13C) &
             d23$age >= 0 &
             d23$age <= 540, , drop = FALSE]

d16 <- data.frame(
  age = suppressWarnings(as.numeric(iso16[["Age (Ma)"]])),
  d13C = suppressWarnings(as.numeric(iso16[["d13Corganic"]]))
)

d16 <- d16[is.finite(d16$age) &
             is.finite(d16$d13C) &
             d16$age > 250 &
             d16$age <= 540, , drop = FALSE]

dorg <- rbind(d23, d16)
dorg$bin <- floor(dorg$age + 0.5)

t_mu <- tapply(dorg$d13C, dorg$bin, mean)
t_sd <- tapply(dorg$d13C, dorg$bin, sd)
t_n <- tapply(dorg$d13C, dorg$bin, length)

org_age <- as.numeric(names(t_mu))
ord <- order(org_age)

org_age <- org_age[ord]
org_mu <- as.numeric(t_mu)[ord]
org_sd <- as.numeric(t_sd)[ord]
org_n <- as.numeric(t_n)[ord]

org_sd[org_n < 2] <- NA_real_

d13CO2_at_org <- approx(x = ages_Ma, y = d13CO2_med, xout = org_age, rule = 2)$y
d13CO2_sd_at_org <- approx(x = ages_Ma, y = co2_sigma, xout = org_age, rule = 2)$y

Dorg_atm <- org_mu - d13CO2_at_org
diff_series <- Dorg_atm

sigma_delta <- sqrt(d13CO2_sd_at_org^2 + org_sd^2)
delta_hw95 <- 1.96 * sigma_delta
delta_hw2sdorg <- 2 * org_sd


####################################################################################################
# Foster et al. CO2 data
####################################################################################################
fos <- read_excel("data/raw/Fosteretal17.xlsx", sheet = 1)

co2 <- data.frame(
  age = suppressWarnings(as.numeric(fos[["age"]])),
  mu = suppressWarnings(as.numeric(fos[["co2"]])),
  u68 = suppressWarnings(as.numeric(fos[["up68"]])),
  l68 = suppressWarnings(as.numeric(fos[["lw68"]])),
  u95 = suppressWarnings(as.numeric(fos[["up95"]])),
  l95 = suppressWarnings(as.numeric(fos[["lw95"]]))
)

co2 <- co2[is.finite(co2$age) &
             is.finite(co2$mu) &
             is.finite(co2$u68) &
             is.finite(co2$l68) &
             is.finite(co2$u95) &
             is.finite(co2$l95) &
             co2$age >= 0 &
             co2$age <= 540, , drop = FALSE]

co2$l95 <- pmax(100, pmin(3500, co2$l95))
co2$u95 <- pmax(100, pmin(3500, co2$u95))
co2$l68 <- pmax(100, pmin(3500, co2$l68))
co2$u68 <- pmax(100, pmin(3500, co2$u68))
co2$mu <- pmax(100, pmin(3500, co2$mu))

co2 <- co2[order(co2$age), , drop = FALSE]

CO2_at_org <- approx(x = co2$age, y = co2$mu, xout = org_age, rule = 2)$y
log_CO2_at_org <- log(CO2_at_org)


####################################################################################################
# Windowed correlations
####################################################################################################
window_cor <- function(age, x, y, window_width = 50, step = 1, min_n = 15, method = "pearson") {
  
  centers <- seq(min(age, na.rm = TRUE) + window_width / 2,
                 max(age, na.rm = TRUE) - window_width / 2,
                 by = step)
  
  out <- data.frame(
    age = centers,
    n = NA_integer_,
    r = NA_real_,
    p = NA_real_
  )
  
  for (i in seq_along(centers)) {
    ok <- is.finite(age) &
      is.finite(x) &
      is.finite(y) &
      age >= centers[i] - window_width / 2 &
      age <= centers[i] + window_width / 2
    
    out$n[i] <- sum(ok)
    
    if (sum(ok) >= min_n) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
      out$r[i] <- unname(ct$estimate)
      out$p[i] <- ct$p.value
    }
  }
  
  out
}

window_width <- 50
window_step <- 1
min_n <- 15

cor_linear <- window_cor(
  age = org_age,
  x = Dorg_atm,
  y = CO2_at_org,
  window_width = window_width,
  step = window_step,
  min_n = min_n,
  method = "pearson"
)

cor_log <- window_cor(
  age = org_age,
  x = Dorg_atm,
  y = log_CO2_at_org,
  window_width = window_width,
  step = window_step,
  min_n = min_n,
  method = "pearson"
)

cor_out <- cor_linear
names(cor_out)[names(cor_out) == "r"] <- "r_Dorg_CO2"
names(cor_out)[names(cor_out) == "p"] <- "p_Dorg_CO2"

cor_out$r_Dorg_logCO2 <- cor_log$r
cor_out$p_Dorg_logCO2 <- cor_log$p


####################################################################################################
# Layout 
####################################################################################################
quartz(width = 11, height = 13)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

layout(matrix(1:7, ncol = 1), heights = c(3.6, 3.0, 3.2, 3.0, 3.0, 0.9, 1.2))
par(xaxs = "i", yaxs = "i")

xlim <- c(358, 0)
ticks10 <- seq(0, 360, by = 10)
ticks50 <- seq(0, 350, by = 50)

y_sides <- c(2, 4, 2, 4, 2)

y_labs <- c(expression(delta^13*C[CO[2]]~"(" * "\u2030 VPDB" * ")"),
            expression(delta^13*C[org]~"(" * "\u2030 VPDB" * ")"),
            expression(Delta[org-atm]~"(" * "\u2030 VPDB" * ")"),
            "windowed r",
            expression(CO[2]~"(ppm)"))


####################################################################################################
## Top panel: d13CO2 with posterior draws 
####################################################################################################

par(mar = c(0.8, 5.2, 2.2, 5.2))   
yl1 <- range(d13CO2_med, na.rm = TRUE); yl1 <- yl1 + c(-0.2, 0.2)

plot(NA, xlim = xlim, ylim = yl1, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# Posterior draws 
for (k in ix_draws) {
  lines(ages_Ma, d13CO2_draws[, k],
        col = adjustcolor("black", 0.15), lwd = 0.6)
}

# Median line
lines(ages_Ma, d13CO2_med, col = "dodgerblue", lwd = 3)

# Axes
axis(3, at = ticks10, labels = FALSE)                
axis(3, at = ticks50, labels = ticks50, tick = FALSE) 
axis(y_sides[1])
mtext(y_labs[1], side = y_sides[1], line = 2.8)


####################################################################################################
## Second panel: plant d13Corg
####################################################################################################
par(mar = c(0.6, 5.2, 0.6, 5.2))

yl2 <- c(-30, -20)

plot(NA, xlim = xlim, ylim = yl2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

points(dorg$age, dorg$d13C, pch = 16, cex = 0.20, col = adjustcolor("grey20", 0.25))
lines(org_age, org_mu, col = "darkgreen", lwd = 3)

axis(y_sides[2])
mtext(y_labs[2], side = y_sides[2], line = 2.8)


####################################################################################################
## Third panel: d13Corg - d13CO2 with uncertainty
####################################################################################################
par(mar = c(0.6, 5.2, 0.6, 5.2))

yl3_core <- range(diff_series, na.rm = TRUE)
yl3_95hi <- diff_series + delta_hw95
yl3_95lo <- diff_series - delta_hw95
yl3_2shi <- diff_series + delta_hw2sdorg
yl3_2slo <- diff_series - delta_hw2sdorg
yl3_all <- range(c(yl3_core, yl3_95hi, yl3_95lo, yl3_2shi, yl3_2slo), na.rm = TRUE)
pad <- diff(range(yl3_all)) * 0.08
yl3 <- c(-27, -10) #yl3_all + c(-pad, pad)

plot(NA, xlim = xlim, ylim = yl3, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# ±2 sd org envelope
ok2s <- is.finite(delta_hw2sdorg)
x2s <- org_age[ok2s]
y2sL <- (diff_series - delta_hw2sdorg)[ok2s]
y2sU <- (diff_series + delta_hw2sdorg)[ok2s]

if (length(x2s)) {
  polygon(c(x2s, rev(x2s)), c(y2sU, rev(y2sL)),
          col = adjustcolor("grey70", 0.45), border = NA)
}

# 95% CI envelope
ok95 <- is.finite(delta_hw95)
x95 <- org_age[ok95]
y95L <- (diff_series - delta_hw95)[ok95]
y95U <- (diff_series + delta_hw95)[ok95]

if (length(x95)) {
  polygon(c(x95, rev(x95)), c(y95U, rev(y95L)),
          col = adjustcolor("coral", 0.45), border = NA)
}

# Mean Delta
lines(org_age, diff_series, lwd = 2, col = "darkred")

axis(y_sides[3])
mtext(y_labs[3], side = y_sides[3], line = 2.8)


####################################################################################################
## Fourth panel: windowed correlations
####################################################################################################
par(mar = c(0.6, 5.2, 0.6, 5.2))

yl4 <- c(-1, 1)

plot(NA, xlim = xlim, ylim = yl4, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

abline(h = 0, lwd = 1, lty = 2)
abline(h = c(-0.5, 0.5), lwd = 1, lty = 3, col = adjustcolor("grey30", 0.6))

lines(cor_out$age, cor_out$r_Dorg_CO2, lwd = 2, col = "black")
lines(cor_out$age, cor_out$r_Dorg_logCO2, lwd = 2, col = "firebrick")

axis(y_sides[4])
mtext(y_labs[4], side = y_sides[4], line = 2.8)

legend("bottomleft",
       inset = 0.02,
       legend = c(expression(Delta[org-atm]~"vs."~CO[2]),
                  expression(Delta[org-atm]~"vs. log"~CO[2])),
       lwd = 2,
       col = c("black", "firebrick"),
       bty = "n",
       cex = 0.9)


####################################################################################################
## Fifth panel: Foster et al. CO2
####################################################################################################
par(mar = c(0.2, 5.2, 0.6, 5.2))

yl5 <- c(100, 3500)

plot(NA, xlim = xlim, ylim = yl5, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

x <- co2$age

polygon(c(x, rev(x)), c(co2$u95, rev(co2$l95)),
        col = adjustcolor("grey70", 0.6), border = NA)

polygon(c(x, rev(x)), c(co2$u68, rev(co2$l68)),
        col = adjustcolor("grey50", 0.6), border = NA)

lines(x, co2$mu, lwd = 2)

axis(y_sides[5])
mtext(y_labs[5], side = y_sides[5], line = 2.8)


####################################################################################################
## GSA timescale strip
####################################################################################################
par(mar = c(0.2, 5.2, 0.2, 5.2))

gsa <- data.frame(
  from = c(0.00, 2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40),
  to = c(2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40, 540.00),
  lab = c("Q", "Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Є"),
  col = c("#FFF2A6", "#FFD37F", "#FFAB70", "#7FC64E", "#34B2C9", "#812B92", "#F04028", "#67A599", "#CB8C37", "#B3E1B6", "#009270", "#7FA056")
)

plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")

for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
  text((gsa$from[i] + gsa$to[i]) / 2, 0.5, gsa$lab[i], cex = 0.9)
}

box()


####################################################################################################
## Shared bottom x axis
####################################################################################################
par(mar = c(3.6, 5.2, 0.2, 5.2))

plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")

axis(1, at = ticks10, labels = FALSE)
axis(1, at = ticks50, labels = ticks50, tick = FALSE)

mtext("Age (Ma)", side = 1, line = 2.4)


