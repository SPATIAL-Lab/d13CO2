####################################################################################################
# Figure 3: wavelet analysis of atmospheric d13C

library(biwavelet)

if (!exists("model.output.root", inherits = FALSE)) {
  source("R/model/d13CO2_RunPaths.R", local = TRUE)
  selected.run <- if (exists("model.run", inherits = FALSE)) model.run else NULL
  model.output.root <- d13CO2_model_run_dir(selected.run, must.exist = TRUE)
}
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"

model.file <- file.path(model.output.root, "inv_out_main.rda")
figure.file <- file.path(figure.output.root, "Figure3.pdf")

load(model.file)

keep <- ages >= 0 & ages <= 540000
age <- ages[keep]/1000
d13CO2.draws <- inv.out$BUGSoutput$sims.list$d13CO2

if (ncol(d13CO2.draws) == length(ages)) {
  d13CO2.draws <- d13CO2.draws[, keep, drop = FALSE]
} else if (nrow(d13CO2.draws) == length(ages)) {
  d13CO2.draws <- t(d13CO2.draws[keep, , drop = FALSE])
} else {
  stop("The atmospheric d13C draws do not match the age vector")
}

d13CO2.median <- apply(d13CO2.draws, 2, median)
o <- order(age)
wavelet.data <- cbind(age[o], d13CO2.median[o])

dt <- mean(diff(wavelet.data[, 1]))
s0 <- 2*dt
dj <- 0.1
J1 <- floor(log2((diff(range(wavelet.data[, 1]))/2)/s0)/dj)

wavelet.result <- wt(wavelet.data, dj = dj, s0 = s0, J1 = J1,
                     mother = "morlet", pad = TRUE, do.sig = TRUE,
                     sig.level = 0.95)

period <- wavelet.result$period
global.power <- rowMeans(wavelet.result$power, na.rm = TRUE)
global.signif <- wt.sig(wavelet.data, dt = dt, scale = wavelet.result$scale,
                        sig.test = 1, sig.level = 0.95,
                        dof = nrow(wavelet.data), mother = "morlet",
                        sigma2 = var(wavelet.data[, 2]))$signif

set.seed(26080703)
n.wavelet.draws <- min(400L, nrow(d13CO2.draws))
draw.rows <- sort(sample(seq_len(nrow(d13CO2.draws)), n.wavelet.draws))

posterior.power <- vapply(draw.rows, function(i) {
  z <- cbind(age[o], d13CO2.draws[i, o])
  w <- wt(z, dj = dj, s0 = s0, J1 = J1, mother = "morlet",
          pad = TRUE, do.sig = FALSE)
  rowMeans(w$power, na.rm = TRUE)
}, numeric(length(period)))

power.median <- apply(posterior.power, 1, median)
power.lower <- apply(posterior.power, 1, quantile, 0.025)
power.upper <- apply(posterior.power, 1, quantile, 0.975)

peak.range <- period >= 55 & period <= 75
peak.index <- which(peak.range)[which.max(power.median[peak.range])]
peak.period <- period[peak.index]

dir.create(dirname(figure.file), recursive = TRUE, showWarnings = FALSE)
pdf(figure.file, width = 7.2, height = 8.5)

op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 1.5, 1.5))

plot(wavelet.result, type = "power", xlab = "Age (Ma)",
     ylab = "Period (Myr)", main = "", plot.coi = TRUE,
     plot.sig = TRUE, lwd.sig = 1.2, col.sig = "white")
abline(h = log2(peak.period), col = "white", lty = 2, lwd = 1.2)
text(max(wavelet.data[, 1]), log2(peak.period),
     paste0(round(peak.period, 1), " Myr"), pos = 2,
     col = "white", font = 2)
mtext("A", side = 3, adj = 0, line = 0.2, font = 2)

x.range <- range(c(0, power.lower, power.upper, global.power,
                   global.signif[peak.index]), finite = TRUE)
x.range[2] <- x.range[2]*1.08
plot(power.median, period, type = "n", log = "y", ylim = c(5, 200),
     xlim = x.range, xlab = "Global wavelet power", ylab = "Period (Myr)")
polygon(c(power.lower, rev(power.upper)), c(period, rev(period)),
        col = adjustcolor("#0072B2", alpha.f = 0.20), border = NA)
lines(power.median, period, col = "#0072B2", lwd = 2)
lines(global.power, period, col = "black", lwd = 1.2)
abline(v = global.signif[peak.index], col = "grey40", lty = 2, lwd = 1.2)
segments(global.signif[peak.index], peak.period,
         power.median[peak.index], peak.period,
         col = "#0072B2", lwd = 2)
points(power.median[peak.index], peak.period, pch = 19,
       col = "#0072B2", cex = 1.1)
text(power.median[peak.index], peak.period,
     paste0(round(peak.period, 1), " Myr"), pos = 4, font = 2)
legend("bottomright",
       c("Posterior median", "95% posterior interval",
         "Median-series power", "95% threshold at 66.1 Myr"),
       col = c("#0072B2", adjustcolor("#0072B2", alpha.f = 0.20),
               "black", "grey40"),
       lwd = c(2, 8, 1.2, 1.2), lty = c(1, 1, 1, 2), bty = "n", cex = 0.82)
mtext("B", side = 3, adj = 0, line = 0.2, font = 2)

par(op)
dev.off()

cat("Dominant 55-75 Myr period:", round(peak.period, 1), "Myr\n")
cat("Posterior median power:", power.median[peak.index], "\n")
cat("95% red-noise threshold:", global.signif[peak.index], "\n")
