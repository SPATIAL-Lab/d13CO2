## =====================  TWO-PANEL FIGURE (+GSA STRIP & X-AXIS)  =====================

ensure_device(11, 9)

# Layout: 1 = d13CO2 (with draws)
#         2 = raw d13C colored by site.index
#         3 = GSA period strip (utility band, not counted as a panel)
#         4 = shared bottom x-axis (utility band)
layout(matrix(1:4, ncol = 1),
       heights = c(3.8, 3.6, 0.9, 1.2))

op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)

par(xaxs = "i", yaxs = "i")

# Keep the “older → younger” direction consistent with your prior figure
xlim   <- c(540, 0)

# X-axis tick logic: ticks every 10 Ma; show numeric labels every 50 Ma
ticks10 <- seq(0, 540, by = 10)
ticks50 <- seq(0, 540, by = 50)

# y-axis sides/labels to match your earlier style
y_side_top <- 2
y_lab_top  <- expression(delta^13*C[CO[2]]~"\u2030")
y_side_bot <- 4
y_lab_bot  <- expression(delta^13*C~"\u2030")

## ===== (1) Top panel: d13CO2 with posterior draws (unchanged look) =====
par(mar = c(0.6, 5.2, 2.2, 5.2))   # top x-axis here; no bottom axis on this panel
yl1 <- range(d13CO2_med, na.rm = TRUE); yl1 <- yl1 + c(-0.2, 0.2)

plot(NA, xlim = xlim, ylim = yl1, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# Posterior draws (thin, semi-transparent)
for (k in ix_draws) {
  lines(ages_Ma, d13CO2_draws[, k], col = adjustcolor("black", 0.15), lwd = 0.6)
}
# Median line
lines(ages_Ma, d13CO2_med, col = "dodgerblue", lwd = 3)

# Axes: keep the x-axis on TOP (x2); y-axis on the right (2)
# Ticks every 10 Ma (no labels), plus labels every 50 Ma
axis(3, at = ticks10, labels = FALSE)                 # tick marks every 10 Ma
axis(3, at = ticks50, labels = ticks50, tick = FALSE) # labeled every 50 Ma
axis(y_side_top)
mtext(y_lab_top, side = y_side_top, line = 2.8)


## ===== (2) Second panel: raw d13C points colored by site.index (match original) =====
par(mar = c(0.2, 5.2, 0.8, 5.2))   # bottom axis still goes in panel 4

# Finite mask only; do NOT restrict ages to [0,540]
ok <- is.finite(prox.in$age) & is.finite(prox.in$d13C) & is.finite(prox.in$site.index)

# xlim exactly like your original line (ensures nothing is culled if ages extend beyond 540)
xlim_raw <- c(max(prox.in$age, na.rm = TRUE), min(prox.in$age, na.rm = TRUE))

# Let y-limits be from the data (like base defaults would do)
yl2 <- range(prox.in$d13C[ok], na.rm = TRUE)

plot(NA, xlim = xlim_raw, ylim = yl2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# Plot points exactly like your example (cex=0.3; use provided colors)
points(prox.in$age[ok], prox.in$d13C[ok],
       col = prox.in$site.index[ok], cex = 0.3)

axis(y_side_bot)
mtext(y_lab_bot, side = y_side_bot, line = 2.8)


## ===== (3) GSA period strip (utility band; not counted as a panel) =====
par(mar = c(0.1, 5.2, 0.1, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
  text((gsa$from[i] + gsa$to[i]) / 2, 0.5, gsa$lab[i], cex = 0.9)
}
box()

## ===== (4) Shared bottom x-axis (x1): ticks every 10 Ma, labels every 50 Ma =====
par(mar = c(3.4, 5.2, 0.1, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

# Ticks every 10 Ma (unlabeled), plus labels every 50 Ma
axis(1, at = ticks10, labels = FALSE)                 # tick marks every 10 Ma
axis(1, at = ticks50, labels = ticks50, tick = FALSE) # labeled every 50 Ma

mtext("Age (Ma)", side = 1, line = 2.6)

# macOS (Quartz window already open):
quartz.save("figure_2panel_raw_d13C.pdf", type = "pdf", device = dev.cur())
