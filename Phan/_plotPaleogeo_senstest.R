###################################################################################################
# Sensitivity test: influence of paleogeographic model — FOUR posteriors in one panel

# Helper to get 95% band + median (mat: iterations x time)
qband <- function(mat) {
  q <- apply(as.matrix(mat), 2, quantile, probs = c(0.025, 0.975, 0.5), na.rm = TRUE)
  list(q025 = q[1, ], q975 = q[2, ], med = q[3, ])
}

# Load ages grid (assumes 'ages' compatible with these runs)
load("Phan/chpc_output3/ages_PALEOMAP.rda")

# Load three existing runs
load("Phan/chpc_output3/inv.out_PALEOMAP.rda");           mcmcout_Sc16 <- inv.out; rm(inv.out)
load("Phan/chpc_output3/inv.out_TorsvikCocks2017.rda");   mcmcout_TC17 <- inv.out; rm(inv.out)
load("Phan/chpc_output3/inv.out_MERDITH2021.rda");        mcmcout_MU22 <- inv.out; rm(inv.out)

# Load the new fourth run (Cao et al., 2024)
load("Phan/chpc_output3/inv.out_CAO2024.rda");            mcmcout_CAO24 <- inv.out; rm(inv.out)

# Extract draws for δ13C_CO2
parm <- "d13CO2"
get_draws  <- function(obj, name) obj[["BUGSoutput"]][["sims.list"]][[name]]
get_median <- function(obj, name) obj[["BUGSoutput"]][["median"]][[name]]

sl_Sc16  <- get_draws(mcmcout_Sc16,  parm)
sl_TC17  <- get_draws(mcmcout_TC17,  parm)
sl_MU22  <- get_draws(mcmcout_MU22,  parm)
sl_CAO24 <- get_draws(mcmcout_CAO24, parm)

med_Sc16  <- get_median(mcmcout_Sc16,  parm)
med_TC17  <- get_median(mcmcout_TC17,  parm)
med_MU22  <- get_median(mcmcout_MU22,  parm)
med_CAO24 <- get_median(mcmcout_CAO24, parm)

# Quantile bands (vectorized)
qb_Sc16  <- qband(sl_Sc16)
qb_TC17  <- qband(sl_TC17)
qb_MU22  <- qband(sl_MU22)
qb_CAO24 <- qband(sl_CAO24)

# Common plotting window
idx   <- 10:549
x_ma  <- ages[idx] / 1e3

# Y limits spanning ALL four series’ 95% CIs over idx
ylim_all <- range(
  qb_Sc16$q025[idx],  qb_Sc16$q975[idx],
  qb_TC17$q025[idx],  qb_TC17$q975[idx],
  qb_MU22$q025[idx],  qb_MU22$q975[idx],
  qb_CAO24$q025[idx], qb_CAO24$q975[idx],
  na.rm = TRUE
)

# Color-blind–safe colors (Okabe–Ito) with semi-transparent envelopes
col_Sc16_fill  <- adjustcolor("#999999", alpha.f = 0.35); col_Sc16_line  <- "#000000"  # grey band, black line
col_TC17_fill  <- adjustcolor("#E69F00", alpha.f = 0.35); col_TC17_line  <- "#E69F00"  # orange
col_MU22_fill  <- adjustcolor("#56B4E9", alpha.f = 0.35); col_MU22_line  <- "#0072B2"  # blue
col_CAO24_fill <- adjustcolor("#009E73", alpha.f = 0.35); col_CAO24_line <- "#009E73"  # bluish green

# Axis label: δ^13C with CO[2] as subscript, and per-mille symbol
ylab_expr <- expression(delta^13*C[CO[2]]~"(" * "‰" * ")")

# Plot
plot(x_ma, med_Sc16[idx], type = "n",
     xlab = "age (Ma)", ylab = ylab_expr,
     xlim = rev(range(x_ma)), ylim = ylim_all)
grid(nx = NA, ny = NULL)

# Scotese (2016)
polygon(c(x_ma, rev(x_ma)), c(qb_Sc16$q025[idx],  rev(qb_Sc16$q975[idx])),
        col = col_Sc16_fill, border = NA)
lines(x_ma, med_Sc16[idx], col = col_Sc16_line, lwd = 2)

# Torsvik & Cocks (2017)
polygon(c(x_ma, rev(x_ma)), c(qb_TC17$q025[idx],  rev(qb_TC17$q975[idx])),
        col = col_TC17_fill, border = NA)
lines(x_ma, med_TC17[idx], col = col_TC17_line, lwd = 2)

# Müller et al. (2022)
polygon(c(x_ma, rev(x_ma)), c(qb_MU22$q025[idx],  rev(qb_MU22$q975[idx])),
        col = col_MU22_fill, border = NA)
lines(x_ma, med_MU22[idx], col = col_MU22_line, lwd = 2)

# Cao et al. (2024) — NEW 4th series
polygon(c(x_ma, rev(x_ma)), c(qb_CAO24$q025[idx], rev(qb_CAO24$q975[idx])),
        col = col_CAO24_fill, border = NA)
lines(x_ma, med_CAO24[idx], col = col_CAO24_line, lwd = 2)

legend("topright",
       legend = c("Scotese (2016)", "Torsvik & Cocks (2017)",
                  "Müller et al. (2022)", "Cao et al. (2024)"),
       col    = c(col_Sc16_line, col_TC17_line, col_MU22_line, col_CAO24_line),
       lwd    = 2,
       bty    = "n",
       seg.len = 3)
