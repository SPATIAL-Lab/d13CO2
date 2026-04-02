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

# Safe par() reset (avoid restoring 'pin','fig','fin','plt','usr','mfg','omd','din','page', etc)
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

# Bin to 1-Myr and take mean
dorg$bin <- floor(dorg$age + 0.5)
org_mean <- aggregate(d13C ~ bin, data = dorg, FUN = mean)
org_age  <- org_mean$bin
org_mu   <- org_mean$d13C

# Difference: d13Corg_mean − d13CO2_median (on org_age)
co2_at_org  <- approx(x = ages_Ma, y = d13CO2_med, xout = org_age, rule = 2)$y
diff_series <- org_mu - co2_at_org

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

## --- GSA-style period strip (example palette; swap if you have official) ---
gsa <- data.frame(
  from = c(  0.00,   2.58,  23.03,  66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40),
  to   = c(  2.58,  23.03,  66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40, 540.00),
  lab  = c("Q", "Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Є"),
  col  = c("#FCD69C", "#F89B5C", "#6C4934", "#8A3EA4", "#33BDE9", "#80CE5C", "#83AD6A", "#00A990", "#B3E4C2", "#CA9547", "#49706F", "#E74D40")
)

## =====================  PLOTTING  =====================

ensure_device(11, 12)

# Layout: 1=d13CO2, 2=d13Corg, 3=Δ (taller), 4=CO2, 5=GSA strip, 6=bottom x-axis
layout(matrix(1:6, ncol = 1),
       heights = c(3.0, 3.0, 4.0, 3.0, 0.9, 1.2))

op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)

par(xaxs = "i", yaxs = "i")
xlim <- c(540, 0)                                 # shared time axis (older → younger)
xticks <- pretty(c(0, 540))

# Alternating y-axis sides (2=right, 1=left) and labels
y_sides <- c(2, 4, 2, 4)                          # panels 1..4
y_labs  <- c(expression(delta^13*C[CO[2]]~"\u2030"),
             expression(delta^13*C[org]~"\u2030"),
             expression(Delta~"\u2030"),
             expression(CO[2]~"(ppm)"))

# ===== (1) d13CO2 =====
par(mar = c(0.5, 4, 2.0, 4))                      # top panel, no bottom axis
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
par(mar = c(0.5, 4, 0.5, 4))
yl2 <- c(-30, -20)
plot(NA, xlim = xlim, ylim = yl2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
points(dorg$age, dorg$d13C, pch = 16, cex = 0.20, col = adjustcolor("grey20", 0.25))
lines(org_age, org_mu, col = "darkgreen", lwd = 3)
axis(y_sides[2]); mtext(y_labs[2], side = y_sides[2], line = 2.8)

# ===== (3) Difference (taller) =====
par(mar = c(0.5, 4, 0.5, 4))
yl3 <- range(diff_series, na.rm = TRUE); pad <- diff(range(yl3))*0.08; yl3 <- yl3 + c(-pad, pad)
plot(NA, xlim = xlim, ylim = yl3, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
lines(org_age, diff_series, lwd = 2)
axis(y_sides[3]); mtext(y_labs[3], side = y_sides[3], line = 2.8)

# ===== (4) CO2 Foster+17 (bands + mean) =====
par(mar = c(0.1, 4, 0.5, 4))                      # small bottom margin; axis goes in panel 6
yl4 <- c(100, 3500)
plot(NA, xlim = xlim, ylim = yl4, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid(nx = NA, ny = NULL)
x <- co2$age
polygon(c(x, rev(x)), c(co2$u95, rev(co2$l95)), col = adjustcolor("grey70", 0.6), border = NA)
polygon(c(x, rev(x)), c(co2$u68, rev(co2$l68)), col = adjustcolor("grey50", 0.6), border = NA)
lines(x, co2$mu, lwd = 2)
axis(y_sides[4]); mtext(y_labs[4], side = y_sides[4], line = 2.8)

# ===== (5) GSA period strip (just above x-axis labels) =====
par(mar = c(0.1, 4, 0.1, 4))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
  text((gsa$from[i]+gsa$to[i])/2, 0.5, gsa$lab[i], cex = 0.9)
}
box()

# ===== (6) Shared bottom x-axis only =====
par(mar = c(3.5, 4, 0.2, 4))
plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
axis(1, at = xticks, labels = xticks)
mtext("Age (Ma)", side = 1, line = 2.4)
