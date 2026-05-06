

## =====================  SETUP / HELPERS  =====================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Please install 'readxl': install.packages('readxl')")
}
library(readxl)

clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

# Robust column picker (handles spaces, NBSP, punctuation, case)
pick_col <- function(df, pattern) {
  nms <- names(df)
  norm_spaces <- function(s) {
    s <- tolower(s)
    gsub("[[:space:]\u00A0]+", "", s, perl = TRUE)
  }
  norm_alnum <- function(s) gsub("[^[:alnum:]]+", "", tolower(s), perl = TRUE)
  key <- norm_spaces(nms); pat <- norm_spaces(pattern)
  hit <- which(key == pat)
  if (!length(hit)) {
    key2 <- norm_alnum(nms); pat2 <- norm_alnum(pattern)
    hit <- which(key2 == pat2)
  }
  if (!length(hit)) stop("Column matching '", pattern, "' not found. Available: ", paste(nms, collapse = ", "))
  nms[hit[1]]
}

as_num <- function(x) suppressWarnings(as.numeric(x))

# Extractor: returns list(ages = Ma, draws = [time x draws])
get_d13CO2 <- function(inv.out, ages_obj = NULL) {
  sims <- NULL
  if (!is.null(inv.out$BUGSoutput$sims.list$d13CO2)) {
    sims <- inv.out$BUGSoutput$sims.list$d13CO2                 # draws x time (typical)
  } else if (!is.null(inv.out$BUGSoutput$sims.matrix)) {
    sm   <- inv.out$BUGSoutput$sims.matrix
    pick <- grepl("^d13CO2\\[\\d+\\]$", colnames(sm))
    if (!any(pick)) stop("No columns like d13CO2[<i>] in sims.matrix.")
    sm   <- sm[, pick, drop = FALSE]
    idx  <- as.integer(sub("^d13CO2\\[(\\d+)\\]$", "\\1", colnames(sm)))
    sims <- sm[, order(idx), drop = FALSE]                      # draws x time
  } else stop("Could not find d13CO2 in inv.out.")
  sims <- as.matrix(sims)
  
  # Ages (prefer provided; else global 'ages' or 'age.indices$age')
  if (!is.null(ages_obj))        ages <- as.numeric(ages_obj)
  else if (exists("ages", envir = .GlobalEnv)) ages <- get("ages", envir = .GlobalEnv)
  else if (exists("age.indices", envir = .GlobalEnv)) ages <- get("age.indices", envir = .GlobalEnv)$age
  else stop("No 'ages' or 'age.indices$age' found.")
  ages <- as.numeric(ages)
  
  ages_Ma <- ages
  if (max(ages_Ma, na.rm = TRUE) > 2e3) ages_Ma <- ages_Ma/1e3   # kyr → Ma
  if (max(ages_Ma, na.rm = TRUE) > 1e4) ages_Ma <- ages_Ma/1e6   # yr  → Ma
  
  # Orient to time x draws
  if (ncol(sims) == length(ages_Ma)) {
    draws <- t(sims)
  } else if (nrow(sims) == length(ages_Ma)) {
    draws <- sims
  } else {
    stop(sprintf("d13CO2 dims = %d x %d, length(ages) = %d.",
                 nrow(sims), ncol(sims), length(ages_Ma)))
  }
  
  keep    <- is.finite(ages_Ma) & ages_Ma >= 0 & ages_Ma <= 540
  ages_Ma <- ages_Ma[keep]; draws <- draws[keep, , drop = FALSE]
  ord     <- order(ages_Ma)
  list(ages = ages_Ma[ord], draws = draws[ord, , drop = FALSE])
}

# Safe par() reset (avoid restoring device-size internals)
safe_par_reset <- function(op) {
  bad <- c("cin","cra","csi","cxy","din","page","usr","mfg","plt","pin","fin","fig","omd")
  keep <- setdiff(names(op), bad)
  par(op[keep])
}

# Open a device reliably; fallback to PDF if headless
ensure_device <- function(w = 11, h = 11, pdf_file = "multiseries_panels.pdf") {
  if (dev.cur() > 1L) {
    while (dev.cur() > 1L) try(dev.off(), silent = TRUE)
  }
  opened <- FALSE
  try({
    if (.Platform$OS.type == "windows") { windows(width = w, height = h); opened <- TRUE }
    else if (Sys.info()[["sysname"]] == "Darwin") { quartz(width = w, height = h); opened <- TRUE }
    else { X11(width = w, height = h); opened <- TRUE }
  }, silent = TRUE)
  if (!opened) { pdf(pdf_file, width = w, height = h); on.exit(dev.off(), add = TRUE) }
}

## =====================  LOAD / PREP DATA  =====================

## --- d13CO2 from inv.out ---
co2_obj       <- get_d13CO2(inv.out)                     # uses global 'ages' if present
ages_Ma       <- co2_obj$ages
d13CO2_draws  <- co2_obj$draws                           # [time x draws]
d13CO2_med    <- apply(d13CO2_draws, 1, median)

## === NEW/UPDATED: CO2 uncertainty summaries (per time) ===
d13CO2_q025   <- apply(d13CO2_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
d13CO2_q975   <- apply(d13CO2_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
d13CO2_sd     <- apply(d13CO2_draws, 1, sd, na.rm = TRUE)         # for propagation
co2_hw95      <- (d13CO2_q975 - d13CO2_q025) / 2                  # half-width (95%)
co2_sigma     <- co2_hw95 / 1.96                                  # approx sd from 95% CI

set.seed(1)
ix_draws      <- sample(seq_len(ncol(d13CO2_draws)), size = min(100, ncol(d13CO2_draws)))

## --- Plant d13Corg (ISOORG23 ≤~251 Ma; ISOORG16 >250 Ma) ---
iso23 <- read_excel("Phan/PhanData/ISOORG23.xlsx", sheet = 1)
iso16 <- read_excel("Phan/PhanData/ISOORG16.xlsx", sheet = 1)

c_age23 <- pick_col(iso23, "1 Myr bin");  c_val23 <- pick_col(iso23, "d13Cp")
d23 <- data.frame(age = as_num(iso23[[c_age23]]), d13C = as_num(iso23[[c_val23]]))
d23 <- d23[is.finite(rowSums(d23)) & d23$age >= 0 & d23$age <= 540, , drop = FALSE]

c_age16 <- pick_col(iso16, "Age (Ma)");   c_val16 <- pick_col(iso16, "d13Corganic")
d16 <- data.frame(age = as_num(iso16[[c_age16]]), d13C = as_num(iso16[[c_val16]]))
d16 <- d16[is.finite(rowSums(d16)) & d16$age > 250 & d16$age <= 540, , drop = FALSE]

dorg <- rbind(d23, d16)
if (!nrow(dorg)) stop("No valid d13Corg rows in the ISOORG files.")


# Bin to 1-Myr and take mean (and sd) — robust against shape issues
dorg$bin <- floor(dorg$age + 0.5)

t_mu <- tapply(dorg$d13C, dorg$bin, mean)
t_sd <- tapply(dorg$d13C, dorg$bin, sd)      # NA when n==1 (that’s fine)
t_n  <- tapply(dorg$d13C, dorg$bin, length)

# Ensure consistent ordering and plain numeric vectors
org_age <- as.numeric(names(t_mu))
ord     <- order(org_age)
org_age <- org_age[ord]
org_mu  <- as.numeric(t_mu)[ord]
org_sd  <- as.numeric(t_sd)[ord]
org_n   <- as.numeric(t_n)[ord]

# Explicitly blank sd where there’s only one point in the bin
org_sd[org_n < 2] <- NA_real_


# Difference: d13CO2_median − d13Corg_mean  (on org_age)  ← CHANGED SIGN HERE
co2_at_org     <- approx(x = ages_Ma, y = d13CO2_med,  xout = org_age, rule = 2)$y

## === NEW/UPDATED: interpolate CO2 sd to org ages for propagation ===
co2_sd_at_org  <- approx(x = ages_Ma, y = co2_sigma,   xout = org_age, rule = 2)$y

diff_series    <- co2_at_org - org_mu

## === NEW/UPDATED: propagate to Δ 95% CI and 2·sd(org) envelope ===
# σΔ^2 = σ_CO2^2 + σ_org^2   (assumed independent)
sigma_delta    <- sqrt(co2_sd_at_org^2 + org_sd^2)
delta_hw95     <- 1.96 * sigma_delta              # half-width for 95% CI of Δ
delta_hw2sdorg <- 2 * org_sd                      # half-width for ±2·sd(org) envelope

## --- CO2 (Foster et al., 2017) with 68/95 bands ---
fos <- read_excel("Phan/PhanData/Fosteretal17.xlsx", sheet = 1)
c_age <- pick_col(fos, "age"); c_mu <- pick_col(fos, "co2")
c_u68 <- pick_col(fos, "up68"); c_l68 <- pick_col(fos, "lw68")
c_u95 <- pick_col(fos, "up95"); c_l95 <- pick_col(fos, "lw95")

co2 <- data.frame(
  age = as_num(fos[[c_age]]),
  mu  = as_num(fos[[c_mu]]),
  u68 = as_num(fos[[c_u68]]), l68 = as_num(fos[[c_l68]]),
  u95 = as_num(fos[[c_u95]]), l95 = as_num(fos[[c_l95]])
)
co2 <- co2[is.finite(rowSums(co2)) & co2$age >= 0 & co2$age <= 540, , drop = FALSE]
co2$l95 <- clamp(co2$l95, 100, 3500); co2$u95 <- clamp(co2$u95, 100, 3500)
co2$l68 <- clamp(co2$l68, 100, 3500); co2$u68 <- clamp(co2$u68, 100, 3500)
co2$mu  <- clamp(co2$mu , 100, 3500)
co2 <- co2[order(co2$age), , drop = FALSE]

## --- GSA-style period strip (embedded once; edit here if you want to tweak) ---
gsa_hex <- c(
  Q  = "#FFF2A6",  Ng = "#FFD37F",  Pg = "#FFAB70",  K  = "#7FC64E",
  J  = "#34B2C9",  Tr = "#812B92",  P  = "#F04028",  C  = "#67A599",
  D  = "#CB8C37",  S  = "#B3E1B6",  O  = "#009270",  "Є"= "#7FA056"
)
gsa <- data.frame(
  from = c(0.00, 2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40),
  to   = c(2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40, 540.00),
  lab  = c("Q","Ng","Pg","K","J","Tr","P","C","D","S","O","Є"),
  col  = NA_character_
)
gsa$col <- unname(gsa_hex[gsa$lab])

## =====================  PLOTTING  =====================

ensure_device(11, 12)

# Layout: 1=d13CO2 (expanded), 2=d13Corg, 3=Δ (slightly smaller), 4=CO2, 5=GSA strip, 6=bottom x-axis
layout(matrix(1:6, ncol = 1),
       heights = c(3.6, 3.0, 3.2, 3.0, 0.9, 1.2))

op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)

par(xaxs = "i", yaxs = "i")
xlim <- c(358, 0)                                 # shared time axis (older → younger)
xticks <- pretty(c(0, 358))

# Alternating y-axis sides (2=right, 1=left) and labels
y_sides <- c(2, 4, 2, 4)                          # panels 1..4
y_labs  <- c(expression(delta^13*C[CO[2]]~"\u2030"),
             expression(delta^13*C[org]~"\u2030"),
             expression(Delta~"\u2030"),
             expression(CO[2]~"(ppm)"))

# ===== (1) d13CO2 =====
par(mar = c(0.8, 5.2, 2.2, 5.2))                  # widen L/R to avoid chopped labels
yl1 <- range(d13CO2_med, na.rm = TRUE); yl1 <- yl1 + c(-0.2, 0.2)
plot(NA, xlim = xlim, ylim = yl1, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
for (k in ix_draws) {
  lines(ages_Ma, d13CO2_draws[, k], col = adjustcolor("black", 0.15), lwd = 0.6)
}
lines(ages_Ma, d13CO2_med, col = "dodgerblue", lwd = 3)
axis(3, at = xticks, labels = xticks)             # x2 axis only here
axis(y_sides[1]); mtext(y_labs[1], side = y_sides[1], line = 2.8)

# ===== (2) d13Corg (truncated to [-30,-20]) =====
par(mar = c(0.6, 5.2, 0.6, 5.2))
yl2 <- c(-30, -20)
plot(NA, xlim = xlim, ylim = yl2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
points(dorg$age, dorg$d13C, pch = 16, cex = 0.20, col = adjustcolor("grey20", 0.25))
lines(org_age, org_mu, col = "darkgreen", lwd = 3)
axis(y_sides[2]); mtext(y_labs[2], side = y_sides[2], line = 2.8)

# ===== (3) Difference with uncertainty (95% CI + ±2·sd(org)) =====
par(mar = c(0.6, 5.2, 0.6, 5.2))

## === NEW/UPDATED: y-limits include both envelopes ===
yl3_core <- range(diff_series, na.rm = TRUE)
yl3_95hi <- diff_series + delta_hw95
yl3_95lo <- diff_series - delta_hw95
yl3_2shi <- diff_series + delta_hw2sdorg
yl3_2slo <- diff_series - delta_hw2sdorg
yl3_all  <- range(c(yl3_core, yl3_95hi, yl3_95lo, yl3_2shi, yl3_2slo), na.rm = TRUE)
pad <- diff(range(yl3_all)) * 0.08
yl3  <- c(10, 27) #yl3_all + c(-pad, pad)

plot(NA, xlim = xlim, ylim = yl3, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)

# (a) ±2·sd(org) envelope (lighter)
ok2s <- is.finite(delta_hw2sdorg)
x2s  <- org_age[ok2s]
y2sL <- (diff_series - delta_hw2sdorg)[ok2s]
y2sU <- (diff_series + delta_hw2sdorg)[ok2s]
if (length(x2s)) {
  polygon(c(x2s, rev(x2s)), c(y2sU, rev(y2sL)),
          col = adjustcolor("grey70", 0.45), border = NA)
}

# (b) 95% CI envelope from propagated σΔ (darker)
ok95 <- is.finite(delta_hw95)
x95  <- org_age[ok95]
y95L <- (diff_series - delta_hw95)[ok95]
y95U <- (diff_series + delta_hw95)[ok95]
if (length(x95)) {
  polygon(c(x95, rev(x95)), c(y95U, rev(y95L)),
          col = adjustcolor("coral", 0.45), border = NA)
}

# Mean Δ
lines(org_age, diff_series, lwd = 2, col = "darkred")

axis(y_sides[3]); mtext(y_labs[3], side = y_sides[3], line = 2.8)

# ===== (4) CO2 Foster+17 (bands + mean) =====
par(mar = c(0.2, 5.2, 0.6, 5.2))
yl4 <- c(100, 3500)
plot(NA, xlim = xlim, ylim = yl4, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
x <- co2$age
polygon(c(x, rev(x)), c(co2$u95, rev(co2$l95)), col = adjustcolor("grey70", 0.6), border = NA)
polygon(c(x, rev(x)), c(co2$u68, rev(co2$l68)), col = adjustcolor("grey50", 0.6), border = NA)
lines(x, co2$mu, lwd = 2)
axis(y_sides[4]); mtext(y_labs[4], side = y_sides[4], line = 2.8)

# ===== (5) GSA period strip (just above x-axis labels) =====
par(mar = c(0.2, 5.2, 0.2, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
  text((gsa$from[i]+gsa$to[i])/2, 0.5, gsa$lab[i], cex = 0.9)
}
box()

# ===== (6) Shared bottom x-axis only =====
par(mar = c(3.6, 5.2, 0.2, 5.2))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
axis(1, at = xticks, labels = xticks)
mtext("Age (Ma)", side = 1, line = 2.4)

# macOS (Quartz window already open):
quartz.save("figure.pdf", type = "pdf", device = dev.cur())
