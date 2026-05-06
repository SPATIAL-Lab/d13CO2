## =====================  STANDALONE HELPERS  =====================

# Need viridisLite for viridis colors
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  stop("Please install 'viridisLite': install.packages('viridisLite')")
}

# Safe par() reset (avoid restoring device-size internals)
safe_par_reset <- function(op) {
  bad <- c("cin","cra","csi","cxy","din","page","usr","mfg",
           "plt","pin","fin","fig","omd")
  keep <- setdiff(names(op), bad)
  par(op[keep])
}

# Open a device reliably; fallback to PDF if headless
ensure_device <- function(w = 11, h = 9, pdf_file = "figure_2panel_raw_d13C_paleolatlon.pdf") {
  if (dev.cur() > 1L) {
    while (dev.cur() > 1L) try(dev.off(), silent = TRUE)
  }
  opened <- FALSE
  try({
    if (.Platform$OS.type == "windows") {
      windows(width = w, height = h); opened <- TRUE
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      quartz(width = w, height = h); opened <- TRUE
    } else {
      X11(width = w, height = h); opened <- TRUE
    }
  }, silent = TRUE)
  if (!opened) {
    pdf(pdf_file, width = w, height = h)
    on.exit(dev.off(), add = TRUE)
  }
}

## =====================  TWO-PANEL FIGURE (+GSA STRIP & X-AXIS)  =====================

# Assumes objects already exist in the environment:
#   d13CO2_med, d13CO2_draws, ages_Ma, ix_draws, prox.in, gsa

ensure_device(11, 9)

# Layout: 1 = d13CO2 (with draws)
#         2 = raw d13C colored by paleolat/paleolon + symbol by category
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

# y-axis sides/labels
y_side_top <- 2
y_lab_top  <- expression(delta^13*C[CO[2]]~"(‰ VPDB)")
y_side_bot <- 4
y_lab_bot  <- expression(delta^13*C[carb]~"(‰ VPDB)")

## ===== (1) Top panel: d13CO2 with posterior draws =====
par(mar = c(0.6, 5.2, 2.2, 5.2))   # top x-axis here; no bottom axis on this panel
yl1 <- range(d13CO2_med, na.rm = TRUE); yl1 <- yl1 + c(-0.2, 0.2)

plot(NA, xlim = xlim, ylim = yl1, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# Posterior draws (thin, semi-transparent)
for (k in ix_draws) {
  lines(ages_Ma, d13CO2_draws[, k],
        col = adjustcolor("black", 0.15), lwd = 0.6)
}
# Median line
lines(ages_Ma, d13CO2_med, col = "dodgerblue", lwd = 3)

# Axes: keep the x-axis on TOP (x2); y-axis on the right (2)
axis(3, at = ticks10, labels = FALSE)                 # tick marks every 10 Ma
axis(3, at = ticks50, labels = ticks50, tick = FALSE) # labeled every 50 Ma
axis(y_side_top)
mtext(y_lab_top, side = y_side_top, line = 2.8)


## ===== (2) Second panel: raw d13C colored by paleolat + paleolon, pch by (re-coded) category =====
par(mar = c(0.2, 5.2, 0.8, 5.2))   # bottom axis still goes in panel 4

# Finite mask only
ok <- is.finite(prox.in$age) &
  is.finite(prox.in$d13C) &
  is.finite(prox.in$site.index) &
  is.finite(prox.in$paleolat) &
  is.finite(prox.in$paleolon)

# ---- Recode and combine categories as requested ----
prox.in$category_rec <- as.character(prox.in$category)

prox.in$category_rec[prox.in$category_rec == "bf"]                      <- "benthic foraminifera"
prox.in$category_rec[prox.in$category_rec == "Ammonite"]                <- "ammonite"
prox.in$category_rec[prox.in$category_rec == "Belemnite"]               <- "belemnite"
prox.in$category_rec[prox.in$category_rec == "Bivalve"]                 <- "bivalve"
prox.in$category_rec[prox.in$category_rec == "Brachiopod calcite"]      <- "brachiopod"

prox.in$category_rec[prox.in$category_rec %in% c("bulk",
                                                 "bulk open ocean",
                                                 "bulk open water")]    <- "bulk carbonate (open ocean)"

prox.in$category_rec[prox.in$category_rec %in% c("bulk marginal sea",
                                                 "bulk marginal sea restricting up section")] <-
  "bulk carbonate (marginal sea)"

prox.in$category_rec[prox.in$category_rec == "bulk semi restricted"]    <- "bulk carbonate (semi-restricted)"

prox.in$category_rec[prox.in$category_rec == "micrite open ocean"]      <- "isolated micrite (open ocean)"

prox.in$category_rec[prox.in$category_rec == "Planktonic foraminifera"] <- "planktonic foraminifera"

# xlim exactly like your original line (ensures nothing is culled if ages extend beyond 540)
xlim_raw <- c(max(prox.in$age, na.rm = TRUE),
              min(prox.in$age, na.rm = TRUE))  # e.g. c(540, 0)

# y-limits: force lower to -10, upper from data
yl2 <- range(prox.in$d13C[ok], na.rm = TRUE)
yl2[1] <- -10

plot(NA, xlim = xlim_raw, ylim = yl2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

## ---- Color by site-mean paleolat & paleolon using viridis (2D via bins) ----

# Compute mean paleolat/paleolon per site
site_lat <- tapply(prox.in$paleolat[ok], prox.in$site.index[ok], mean, na.rm = TRUE)
site_lon <- tapply(prox.in$paleolon[ok], prox.in$site.index[ok], mean, na.rm = TRUE)

# Define bins for legend and for mapping
lat_bins <- c(0, 30, 60, 90)                 # absolute latitude
lon_bins <- c(-180, -90, 0, 90, 180)         # longitude

n_lat <- length(lat_bins) - 1
n_lon <- length(lon_bins) - 1

# Precompute a viridis color grid over (lat_bin, lon_bin)
# low lat (equator) → yellow-ish; high lat (poles) → blue-ish
col_vec   <- viridisLite::viridis(n_lat * n_lon, option = "D")
grid_cols <- matrix(col_vec, nrow = n_lat, ncol = n_lon, byrow = TRUE)

# Function to map continuous lat/lon to nearest bin and get color
latlon_to_col <- function(lat, lon) {
  a_lat <- abs(lat)
  # bin indices
  lat_idx <- findInterval(a_lat, lat_bins, rightmost.closed = TRUE)
  lon_idx <- findInterval(lon,   lon_bins, rightmost.closed = TRUE)
  # clamp to 1..n
  lat_idx[lat_idx < 1] <- 1;    lat_idx[lat_idx > n_lat] <- n_lat
  lon_idx[lon_idx < 1] <- 1;    lon_idx[lon_idx > n_lon] <- n_lon
  # flip latitude index so low lat uses high viridis (yellow), high lat uses low (blue)
  lat_idx2 <- n_lat - lat_idx + 1
  grid_cols[cbind(lat_idx2, lon_idx)]
}

# Site-level colors
site_ids  <- names(site_lat)
site_cols <- latlon_to_col(site_lat, site_lon)
names(site_cols) <- site_ids

# Assign color and pch to each observation via its site.index and recoded category
col_pts <- site_cols[as.character(prox.in$site.index[ok])]

# Symbol by recoded category
cat_vec    <- as.character(prox.in$category_rec[ok])
cat_levels <- sort(unique(cat_vec))

# Ensure we have enough pch values for all categories (recycle if needed)
base_pch <- c(16, 17, 15, 3, 4, 8, 1, 2, 6, 7)
pch_seq  <- rep(base_pch, length.out = length(cat_levels))
cat_pch_map <- setNames(pch_seq, cat_levels)

# Swap symbols for planktonic foraminifera and ammonite
if (all(c("planktonic foraminifera", "ammonite") %in% names(cat_pch_map))) {
  tmp_pch <- cat_pch_map["planktonic foraminifera"]
  cat_pch_map["planktonic foraminifera"] <- cat_pch_map["ammonite"]
  cat_pch_map["ammonite"]                <- tmp_pch
}

pch_pts <- cat_pch_map[cat_vec]

# Plot points (colored by lat/lon, shaped by recoded category)
points(prox.in$age[ok], prox.in$d13C[ok],
       col = col_pts, cex = 0.3, pch = pch_pts)

axis(y_side_bot)
mtext(y_lab_bot, side = y_side_bot, line = 2.8)

## ---- 2D legend for paleolat / paleolon mapping (true lower-left, ultra-compact) ----

lat_labels <- c("0–30°", "30–60°", "60–90°")
lon_labels <- c("180–90W", "90–0W", "0–90E", "90–180E")  # W on left, E on right

# Because xlim_raw is in descending order (e.g., c(540, 0)),
# "left" on the plot corresponds to the first element.
x_left  <- xlim_raw[1]   # older age, left edge in plot
x_right <- xlim_raw[2]   # younger age, right edge
y_min   <- yl2[1]
y_max   <- yl2[2]

x_leg_left   <- x_left + 0.02 * (x_right - x_left)
x_leg_right  <- x_left + 0.26 * (x_right - x_left)
y_leg_bottom <- y_min  + 0.02 * (y_max  - y_min)
y_leg_top    <- y_min  + 0.24 * (y_max  - y_min)

dx <- (x_leg_right - x_leg_left) / (n_lon + 0.1)
dy <- (y_leg_top   - y_leg_bottom) / (n_lat + 0.1)

# Draw colored grid cells (no outer box)
for (i in seq_len(n_lat)) {
  for (j in seq_len(n_lon)) {
    lat_mid <- mean(lat_bins[i:(i+1)])
    lon_mid <- mean(lon_bins[j:(j+1)])
    col_ij  <- latlon_to_col(lat_mid, lon_mid)
    
    x0 <- x_leg_left  + dx * (j - 0.0)
    x1 <- x_leg_left  + dx * (j + 0.85)
    y0 <- y_leg_bottom + dy * (i - 0.0)
    y1 <- y_leg_bottom + dy * (i + 0.85)
    
    rect(x0, y0, x1, y1, col = col_ij, border = NA)
  }
}

# Lat labels (rows) on left of grid – right up next to blocks
for (i in seq_len(n_lat)) {
  y_mid <- y_leg_bottom + dy * (i + 0.32)
  text(x_leg_left - dx * 0.015, y_mid, lat_labels[i],
       adj = c(1, 0.5), cex = 1.0)
}

# Lon labels (columns) below grid – right up under the blocks
for (j in seq_len(n_lon)) {
  x_mid <- x_leg_left + dx * (j + 0.32)
  text(x_mid, y_leg_bottom - dy * 0.02, lon_labels[j],
       adj = c(0.5, 1), cex = 1.0, srt = 30)
}

## ---- Small symbol legend for recoded category (lower-right, slightly larger, no title) ----

legend("bottomright",
       inset = 0.02,
       legend = cat_levels,
       pch    = cat_pch_map[cat_levels],
       col    = "black",
       pt.cex = 1.05,
       cex    = 1.05,
       bty    = "n",
       title  = NULL)

## ===== (3) GSA period strip (utility band; not counted as a panel) =====
par(mar = c(0.1, 5.2, 0.1, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")
for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
  text((gsa$from[i] + gsa$to[i]) / 2, 0.5, gsa$lab[i], cex = 0.9)
}
box()

## ===== (4) Shared bottom x-axis (x1): ticks every 10 Ma, labels every 50 Ma =====
# Move this band closer to the GSA strip by reducing the top margin
par(mar = c(1.0, 5.2, 0.1, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")

axis(1, at = ticks10, labels = FALSE)                 # tick marks every 10 Ma
axis(1, at = ticks50, labels = ticks50, tick = FALSE) # labeled every 50 Ma

mtext("Age (Ma)", side = 1, line = 2.6)

# macOS (Quartz window already open):
quartz.save("figure_2panel_raw_d13C_paleolatlon_final.pdf",
            type = "pdf", device = dev.cur())
