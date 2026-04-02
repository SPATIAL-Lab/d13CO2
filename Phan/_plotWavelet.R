# Wavelet analysis for d13Ca
# Focus on 45, 60, 70, 80, 180 Myr; shaded bands; ridge traces; threshold
############################################################

library(biwavelet)

# 1. Load data --------------------------------------------------------------

dat_raw <- read.csv("d13Ca.csv", stringsAsFactors = FALSE)

# Sort by age and simplify
dat <- dat_raw[order(dat_raw$ages), c("ages", "d13Ca")]

# Ages are in kyr; convert to Myr
dat$time_Myr <- dat$ages / 1000

# Sampling interval (Myr)
dt <- mean(diff(dat$time_Myr))

# Prepare input matrix (time, value)
wt_input <- cbind(dat$time_Myr, dat$d13Ca)

# 2. Wavelet setup ----------------------------------------------------------

s0 <- 2 * dt
total_length <- max(dat$time_Myr) - min(dat$time_Myr)
max_period   <- total_length / 2
dj <- 0.1
J1 <- floor(log2(max_period / s0) / dj)

# 3. Run continuous wavelet transform ---------------------------------------

wt_res <- wt(
  wt_input,
  dj     = dj,
  s0     = s0,
  J1     = J1,
  mother = "morlet",
  pad    = TRUE,
  do.sig = TRUE
)

# 4. Target periods and global spectrum -------------------------------------

# Periods we want to emphasize (in Myr)
target_periods <- c(45, 60, 70, 80, 180)

# Global wavelet spectrum
global_power <- apply(wt_res$power, 1, mean, na.rm = TRUE)
periods      <- wt_res$period   # in Myr (because time was in Myr)

# Threshold for "significant" power in 5–200 Myr band (e.g. 90th percentile)
band_5_200 <- which(periods >= 5 & periods <= 200)
thresh <- if (length(band_5_200) > 0) {
  stats::quantile(global_power[band_5_200], 0.90, na.rm = TRUE)
} else {
  NA_real_
}

# Get indices for target periods (nearest scale)
target_idx <- sapply(target_periods, function(p) which.min(abs(periods - p)))

# 5. Plot -------------------------------------------------------------------

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

## (A) Time series ----------------------------------------------------------

plot(dat$time_Myr, dat$d13Ca, type = "l",
     xlab = "Age (Myr)",
     ylab = expression(delta^{13}*C[a]~("\u2030")),
     main = expression("d"^{13}*"C"[a]*" time series"),
     lwd = 1.5)
grid(col = "grey85", lty = "dotted")

## (B) Wavelet power spectrum with shaded bands + labels --------------------

plot(
  wt_res,
  type = "power",
  main = "Wavelet power spectrum",
  xlab = "Age (Myr)",
  ylab = "Period (Myr)"
)

usr <- par("usr")

# Shaded bands around each target period
for (p in target_periods) {
  # Narrower band for shorter periods, wider for 180 Myr
  half_width <- if (p >= 150) 10 else 5
  rect(xleft = usr[1], xright = usr[2],
       ybottom = p - half_width, ytop = p + half_width,
       col = grDevices::adjustcolor("white", alpha.f = 0.15),
       border = NA)
}

# Re-draw thin center lines for each target period
for (p in target_periods) {
  segments(x0 = usr[1], x1 = usr[2], y0 = p,
           col = grDevices::adjustcolor("white", alpha.f = 0.7),
           lty = 2, lwd = 1.5)
}

# Explicit labels at 60 and 180 Myr
text(usr[2], 60,
     labels = "60 Myr",
     pos = 2, xpd = NA, cex = 1.0, font = 2, col = "white")
text(usr[2], 180,
     labels = "180 Myr",
     pos = 2, xpd = NA, cex = 1.0, font = 2, col = "white")

## (C) Global wavelet spectrum with shaded bands + threshold ----------------

plot(global_power, periods, type = "l",
     xlab = "Global wavelet power",
     ylab = "Period (Myr)",
     main = "Global wavelet spectrum",
     log  = "y",
     lwd  = 1.5)
grid(col = "grey85", lty = "dotted")

usrC <- par("usr")

# Shaded bands in global spectrum for each target period
for (p in target_periods) {
  half_width <- if (p >= 150) 10 else 5
  rect(xleft = usrC[1], xright = usrC[2],
       ybottom = p - half_width, ytop = p + half_width,
       col = grDevices::adjustcolor("grey", alpha.f = 0.12),
       border = NA)
}

# Re-plot line so it appears over shading
lines(global_power, periods, lwd = 1.5)

# Threshold line
if (!is.na(thresh)) {
  abline(v = thresh, col = "blue", lty = 2, lwd = 1.5)
  text(thresh, periods[min(band_5_200)],
       labels = "90% threshold",
       pos = 4, col = "blue", cex = 0.8)
}

# Mark and label target periods
cols <- c("red", "orange", "darkgreen", "purple", "brown")
for (i in seq_along(target_periods)) {
  idx <- target_idx[i]
  p   <- target_periods[i]
  points(global_power[idx], periods[idx],
         pch = 19, col = cols[i], cex = if (p == 60) 1.4 else 1.1)
  text(global_power[idx], periods[idx],
       labels = paste0(p, " Myr"),
       pos = 4, col = cols[i], cex = 0.8, font = if (p == 60) 2 else 1)
}

## (D) Power ridge traces at each target period -----------------------------

# For each target period, extract power vs time at nearest scale
ridge_power <- lapply(target_idx, function(k) wt_res$power[k, ])

# Plot the first ridge to set axes
plot(dat$time_Myr, ridge_power[[1]], type = "l",
     xlab = "Age (Myr)",
     ylab = "Wavelet power",
     main = "Power at selected periods",
     lwd  = 1.5, col = cols[1])
grid(col = "grey85", lty = "dotted")

# Add remaining ridges
for (i in 2:length(target_periods)) {
  lines(dat$time_Myr, ridge_power[[i]],
        col = cols[i], lwd = 1.3)
}

legend("topright",
       legend = paste(target_periods, "Myr"),
       col    = cols,
       lwd    = 1.5,
       cex    = 0.75,
       bg     = "white")

par(old_par)
