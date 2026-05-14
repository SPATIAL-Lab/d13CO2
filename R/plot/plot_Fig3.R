

####################################################################################################
####################################################################################################
# Script to generate Figure 3 in the manuscript - wavelet analysis for d13Ca

# d13CO2 (with draws), raw d13C colored by paleolat/paleolon + symbol by category; GSA period strip
# and shared bottom x-axis 

# Assumes objects already exist in the environment:
# inv.out
# ages

load("output/model_runs/chpc_260511/output/inv.out_main.rda")
load("output/model_runs/chpc_260511/output/ages_main.rda")

####################################################################################################
####################################################################################################

library(biwavelet)

############################################################
# Extract d13Ca
############################################################

d13CO2_draws <- as.matrix(inv.out$BUGSoutput$sims.list$d13CO2)

if (ncol(d13CO2_draws) == length(ages)) {
  d13CO2_draws <- t(d13CO2_draws)
}

d13CO2_med <- apply(d13CO2_draws, 1, median)

dat <- data.frame(
  time_Myr = ages / 1000,
  d13Ca = d13CO2_med
)

dat <- dat[order(dat$time_Myr), ]

dt <- mean(diff(dat$time_Myr))
wt_input <- cbind(dat$time_Myr, dat$d13Ca)

############################################################
# Wavelet setup
############################################################

s0 <- 2 * dt
dj <- 0.1
max_period <- diff(range(dat$time_Myr)) / 2
J1 <- floor(log2(max_period / s0) / dj)

wt_res <- wt(wt_input,
             dj = dj,
             s0 = s0,
             J1 = J1,
             mother = "morlet",
             pad = TRUE,
             do.sig = TRUE)

############################################################
# Periods and global spectrum
############################################################

target_periods <- c(45, 60, 70, 80, 180)

periods <- wt_res$period
global_power <- apply(wt_res$power, 1, mean, na.rm = TRUE)

band_5_200 <- which(periods >= 5 & periods <= 200)
thresh <- quantile(global_power[band_5_200], 0.90, na.rm = TRUE)

target_idx <- sapply(target_periods, function(p) {
  which.min(abs(periods - p))
})

cols <- c("red", "orange", "darkgreen", "purple", "brown")

############################################################
# Plot 1: Wavelet power spectrum + global wavelet spectrum
############################################################

quartz(width = 8.5, height = 9)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 2))

############################################################
# Wavelet power spectrum
############################################################

plot(wt_res,
     type = "power",
     main = "Wavelet power spectrum",
     xlab = "Age (Ma)",
     ylab = "Period (Myr)")

usr <- par("usr")

for (p in target_periods) {
  half_width <- if (p >= 150) 10 else 5
  
  rect(usr[1], p - half_width, usr[2], p + half_width,
       col = adjustcolor("white", 0.15), border = NA)
  
  segments(usr[1], p, usr[2], p,
           col = adjustcolor("white", 0.7), lty = 2, lwd = 1.5)
}

text(usr[2], 60, "60 Myr", pos = 2, xpd = NA, cex = 1.0, font = 2, col = "white")
text(usr[2], 180, "180 Myr", pos = 2, xpd = NA, cex = 1.0, font = 2, col = "white")

############################################################
# Global wavelet spectrum
############################################################

plot(global_power, periods, type = "l",
     xlab = "Global wavelet power",
     ylab = "Period (Myr)",
     main = "Global wavelet spectrum",
     log = "y",
     lwd = 1.5)

grid(col = "grey85", lty = "dotted")

usr <- par("usr")

for (p in target_periods) {
  half_width <- if (p >= 150) 10 else 5
  
  rect(usr[1], p - half_width, usr[2], p + half_width,
       col = adjustcolor("grey", 0.12), border = NA)
}

lines(global_power, periods, lwd = 1.5)
abline(v = thresh, col = "blue", lty = 2, lwd = 1.5)

for (i in seq_along(target_periods)) {
  idx <- target_idx[i]
  
  points(global_power[idx], periods[idx],
         pch = 19, col = cols[i], cex = if (target_periods[i] == 60) 1.4 else 1.1)
  
  text(global_power[idx], periods[idx],
       labels = paste0(target_periods[i], " Myr"),
       pos = 4, col = cols[i], cex = 0.8, font = if (target_periods[i] == 60) 2 else 1)
}


############################################################
# Plot 2: d13Ca time series + power at selected periods
############################################################

quartz(width = 8.5, height = 9)

par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 2))

############################################################
# Time series
############################################################

plot(dat$time_Myr, dat$d13Ca, type = "l",
     xlab = "Age (Ma)",
     ylab = expression(delta^13*C[a]~"(‰ VPDB)"),
     main = expression(delta^13*C[a]*" time series"),
     lwd = 1.5)

grid(col = "grey85", lty = "dotted")

############################################################
# Power at target periods
############################################################

ridge_power <- lapply(target_idx, function(k) {
  wt_res$power[k, ]
})

plot(dat$time_Myr, ridge_power[[1]], type = "l",
     xlab = "Age (Ma)",
     ylab = "Wavelet power",
     main = "Power at selected periods",
     lwd = 1.5,
     col = cols[1])

grid(col = "grey85", lty = "dotted")

for (i in 2:length(target_periods)) {
  lines(dat$time_Myr, ridge_power[[i]], col = cols[i], lwd = 1.3)
}

legend("topright",
       legend = paste(target_periods, "Myr"),
       col = cols,
       lwd = 1.5,
       cex = 0.75,
       bg = "white")


