## =====================  THREE-PANEL: d13CO2 post | GMST post | GMST priors  =====================

## --- Helpers (base R) ---
safe_par_reset <- function(op) {
  bad <- c("cin","cra","csi","cxy","din","page","usr","mfg","plt","pin","fin","fig","omd")
  par(op[setdiff(names(op), bad)])
}
qband <- function(mat) {  # mat: iterations x time
  q <- apply(as.matrix(mat), 2, quantile, probs = c(0.025, 0.975, 0.5), na.rm = TRUE)
  list(q025 = q[1, ], q975 = q[2, ], med = q[3, ])
}
interp_to <- function(x_src, y_src, x_tgt) {
  # Interpolate to common target grid; return NA outside source domain
  approx(x_src, y_src, xout = x_tgt, rule = 1)$y
}

## =====================  BOTTOM PANEL DATA (PRIORS)  =====================
age.min <- 0
age.max <- 540000
step.int <- 1000
n.spinup <- 10

# uniform SDs
GMST_sd_Li22 <- 5
toff_sd_uniform <- 2
toff_sd_uniform_bot <- 1

# load proxy data
prox.in <- as.data.frame(read.csv(file = "Phan/PhanData/PhanCompWithTemp_PALEOMAP.csv"))
prox.in <- cbind(prox.in[,1:7], prox.in[,9:10], prox.in[,21:27],
                 rep(x = toff_sd_uniform, times = nrow(prox.in)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category",
                    "paleolon","paleolat", "MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                    "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA", "temp_offset_sd")

# ages in kyr; window; build regular age grids
prox.in$age <- prox.in$age * 1e3
prox.in <- prox.in[prox.in$age >= (age.min) & prox.in$age <= (age.max), ]
prox.in <- transform(prox.in,
                     ai = n.spinup + as.numeric(1 + floor((max(prox.in$age) - prox.in$age) / step.int)))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE), ]
ages.short <- seq(from = max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
ages       <- seq(from = n.spinup*step.int + max(prox.in$age), to = min(prox.in$age),
                  by = -1*step.int) - 0.5*step.int
ages_Ma    <- ages / 1e3

# symmetric PhanDA sd
PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                (prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo)) / 2

# interpolate GMST priors to model ages
GMST.m    <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = ages, rule = 2)$y
GMST.sd   <- approx(prox.in$age, PhanDA_sd,          xout = ages, rule = 2)$y
GMST.low  <- GMST.m - 2*GMST.sd
GMST.hi   <- GMST.m + 2*GMST.sd

GMST.m.Li   <- approx(prox.in$age, prox.in$GMST_Li22, xout = ages, rule = 2)$y
GMST.sd.Li  <- rep(GMST_sd_Li22, times = length(ages))
GMST.low.Li <- GMST.m.Li - 2*GMST.sd.Li
GMST.hi.Li  <- GMST.m.Li + 2*GMST.sd.Li

## =====================  LOAD MCMC OBJECTS FOR TOP/MIDDLE  =====================
load("Phan/chpc_output2/ages_PALEOMAP.rda")      # may set/confirm 'ages'
load("Phan/chpc_output2/inv.out_LiGMST.rda"); mcmcout_Li <- inv.out; rm(inv.out)
load("Phan/chpc_output2/inv.out_PALEOMAP.rda")   # loads 'inv.out' for PhanDA-driven inversion

## Target window (common to all panels)
idx_d13CO2 <- 10:549
x_common   <- ages_Ma[idx_d13CO2]                # Ma, decreasing (older→younger reversed in xlim)

## =====================  QUANTILES: d13CO2 (top) and GMST (middle)  =====================
# d13CO2 (two models)
d13_Li_mat   <- mcmcout_Li[["BUGSoutput"]][["sims.list"]][["d13CO2"]]   # iters x time
d13_Phan_mat <- inv.out     [["BUGSoutput"]][["sims.list"]][["d13CO2"]]
d13_Li   <- qband(d13_Li_mat)
d13_Phan <- qband(d13_Phan_mat)

# GMST (two models) – align to common x grid
gmst_Li_mat   <- mcmcout_Li[["BUGSoutput"]][["sims.list"]][["GMST"]]
gmst_Phan_mat <- inv.out     [["BUGSoutput"]][["sims.list"]][["GMST"]]
gmst_Li_raw   <- qband(gmst_Li_mat)
gmst_Phan_raw <- qband(gmst_Phan_mat)

# Native age vectors for each GMST series (assume alignment with head of ages_Ma)
age_gmst_Li   <- ages_Ma[seq_len(length(gmst_Li_raw$med))]
age_gmst_Phan <- ages_Ma[seq_len(length(gmst_Phan_raw$med))]

# Interpolate each GMST quantile to the common x grid (no extrapolation)
gmst_Li <- list(
  q025 = interp_to(age_gmst_Li,   gmst_Li_raw$q025, x_common),
  q975 = interp_to(age_gmst_Li,   gmst_Li_raw$q975, x_common),
  med  = interp_to(age_gmst_Li,   gmst_Li_raw$med,  x_common)
)
gmst_Phan <- list(
  q025 = interp_to(age_gmst_Phan, gmst_Phan_raw$q025, x_common),
  q975 = interp_to(age_gmst_Phan, gmst_Phan_raw$q975, x_common),
  med  = interp_to(age_gmst_Phan, gmst_Phan_raw$med,  x_common)
)

## =====================  PLOTTING  =====================
# Colors
# Top & middle (shared palette)
col_phan_ts_fill <- adjustcolor("#E66101", 0.35)  # amber (PhanDA)
col_phan_ts_line <- "#E66101"
col_li_ts_fill   <- adjustcolor("#5E3C99", 0.35)  # deep violet (Li)
col_li_ts_line   <- "#5E3C99"

# Bottom (priors): darker red (PhanDA) and blue/navy (Li)
col_phan_prior_fill <- adjustcolor("#A50026", 0.35)  # dark red
col_phan_prior_line <- "#A50026"
col_li_prior_fill   <- adjustcolor("#4575B4", 0.35)  # navy-ish blue
col_li_prior_line   <- "#4575B4"

op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)
par(xaxs = "i", yaxs = "i")
layout(matrix(1:3, ncol = 1), heights = c(3.2, 3.2, 3.0))

xlim_all <- rev(range(x_common, na.rm = TRUE))   # reversed Ma

## ----- (1) TOP: d13C_CO2 posteriors -----
par(mar = c(0.5, 5.0, 2.2, 1.8))
ylab_d13 <- expression(delta^13*C[CO[2]]~"(" * "\u2030" * ")")

x_ma <- x_common

ylim_top <- range(d13_Phan$q025[idx_d13CO2], d13_Phan$q975[idx_d13CO2],
                  d13_Li$q025[idx_d13CO2],   d13_Li$q975[idx_d13CO2], na.rm = TRUE)

plot(NA, xlim = xlim_all, ylim = ylim_top, xlab = "", ylab = ylab_d13, xaxt = "n", yaxt = "s")
grid(nx = NA, ny = NULL)

# PhanDA
polygon(c(x_ma, rev(x_ma)),
        c(d13_Phan$q025[idx_d13CO2], rev(d13_Phan$q975[idx_d13CO2])),
        col = col_phan_ts_fill, border = NA)
lines(x_ma, d13_Phan$med[idx_d13CO2], col = col_phan_ts_line, lwd = 1.2)

# Li22
polygon(c(x_ma, rev(x_ma)),
        c(d13_Li$q025[idx_d13CO2], rev(d13_Li$q975[idx_d13CO2])),
        col = col_li_ts_fill, border = NA)
lines(x_ma, d13_Li$med[idx_d13CO2], col = col_li_ts_line, lwd = 1.2)

legend("topright",
       legend = c("PhanDA GMST", "Li22 GMST"),
       fill   = c(col_phan_ts_fill, col_li_ts_fill),
       border = NA, bty = "n")

## ----- (2) MIDDLE: GMST posteriors (same colors as top; y >= 0) -----
par(mar = c(0.5, 5.0, 0.8, 1.8))
ylab_gmst <- expression(paste("GMST (", degree, "C)"))
x_ma_g <- x_common

# Force y >= 0
ymax_mid <- max(gmst_Phan$q975, gmst_Li$q975, na.rm = TRUE)
ylim_mid <- c(0, ymax_mid)

plot(NA, xlim = xlim_all, ylim = ylim_mid, xlab = "", ylab = ylab_gmst, xaxt = "n", yaxt = "s")
grid(nx = NA, ny = NULL)

# PhanDA (amber)
polygon(c(x_ma_g, rev(x_ma_g)),
        c(gmst_Phan$q025, rev(gmst_Phan$q975)),
        col = col_phan_ts_fill, border = NA)
lines(x_ma_g, gmst_Phan$med, col = col_phan_ts_line, lwd = 1.2)

# Li22 (violet)
polygon(c(x_ma_g, rev(x_ma_g)),
        c(gmst_Li$q025, rev(gmst_Li$q975)),
        col = col_li_ts_fill, border = NA)
lines(x_ma_g, gmst_Li$med, col = col_li_ts_line, lwd = 1.2)

legend("topright",
       legend = c("PhanDA GMST", "Li22 GMST"),
       fill   = c(col_phan_ts_fill, col_li_ts_fill),
       border = NA, bty = "n")

## ----- (3) BOTTOM: GMST priors (PhanDA vs Li; new colors, legend placed in blank area) -----
par(mar = c(3.6, 5.0, 0.8, 1.8))
ylim_bot <- range(GMST.low[idx_d13CO2], GMST.hi[idx_d13CO2],
                  GMST.low.Li[idx_d13CO2], GMST.hi.Li[idx_d13CO2], na.rm = TRUE)

plot(x_common, GMST.m[idx_d13CO2], type = "n",
     xlab = "age (Ma)", ylab = ylab_gmst, xlim = xlim_all, ylim = ylim_bot)
grid(nx = NA, ny = NULL)

# PhanDA prior (dark red)
polygon(x = c(x_common, rev(x_common)),
        y = c(GMST.low[idx_d13CO2], rev(GMST.hi[idx_d13CO2])),
        col = col_phan_prior_fill, border = NA)
lines(x_common, GMST.m[idx_d13CO2], col = col_phan_prior_line, lwd = 1.5)

# Li22 prior (navy)
polygon(x = c(x_common, rev(x_common)),
        y = c(GMST.low.Li[idx_d13CO2], rev(GMST.hi.Li[idx_d13CO2])),
        col = col_li_prior_fill, border = NA)
lines(x_common, GMST.m.Li[idx_d13CO2], col = col_li_prior_line, lwd = 1.5)

# Legend placed in a likely blank space (top-left) with a tiny inset
legend("topleft",
       legend = c("PhanDA GMST prior", "Li22 GMST prior"),
       fill   = c(col_phan_prior_fill, col_li_prior_fill),
       border = NA, bty = "n", inset = 0.01)



