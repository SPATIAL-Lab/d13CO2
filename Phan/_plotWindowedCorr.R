## =====================  SETUP / HELPERS  =====================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Please install 'readxl': install.packages('readxl')")
}
library(readxl)

clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

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

get_d13CO2 <- function(inv.out, ages_obj = NULL) {
  sims <- NULL
  
  if (!is.null(inv.out$BUGSoutput$sims.list$d13CO2)) {
    sims <- inv.out$BUGSoutput$sims.list$d13CO2
  } else if (!is.null(inv.out$BUGSoutput$sims.matrix)) {
    sm <- inv.out$BUGSoutput$sims.matrix
    pick <- grepl("^d13CO2\\[\\d+\\]$", colnames(sm))
    if (!any(pick)) stop("No columns like d13CO2[<i>] in sims.matrix.")
    sm <- sm[, pick, drop = FALSE]
    idx <- as.integer(sub("^d13CO2\\[(\\d+)\\]$", "\\1", colnames(sm)))
    sims <- sm[, order(idx), drop = FALSE]
  } else {
    stop("Could not find d13CO2 in inv.out.")
  }
  
  sims <- as.matrix(sims)
  
  if (!is.null(ages_obj)) {
    ages <- as.numeric(ages_obj)
  } else if (exists("ages", envir = .GlobalEnv)) {
    ages <- get("ages", envir = .GlobalEnv)
  } else if (exists("age.indices", envir = .GlobalEnv)) {
    ages <- get("age.indices", envir = .GlobalEnv)$age
  } else {
    stop("No 'ages' or 'age.indices$age' found.")
  }
  
  ages <- as.numeric(ages)
  ages_Ma <- ages
  
  if (max(ages_Ma, na.rm = TRUE) > 2e3) ages_Ma <- ages_Ma / 1e3
  if (max(ages_Ma, na.rm = TRUE) > 1e4) ages_Ma <- ages_Ma / 1e6
  
  if (ncol(sims) == length(ages_Ma)) {
    draws <- t(sims)
  } else if (nrow(sims) == length(ages_Ma)) {
    draws <- sims
  } else {
    stop(sprintf("d13CO2 dims = %d x %d, length(ages) = %d.",
                 nrow(sims), ncol(sims), length(ages_Ma)))
  }
  
  keep <- is.finite(ages_Ma) & ages_Ma >= 0 & ages_Ma <= 540
  ages_Ma <- ages_Ma[keep]
  draws <- draws[keep, , drop = FALSE]
  
  ord <- order(ages_Ma)
  list(ages = ages_Ma[ord], draws = draws[ord, , drop = FALSE])
}

safe_par_reset <- function(op) {
  bad <- c("cin","cra","csi","cxy","din","page","usr","mfg","plt","pin","fin","fig","omd")
  keep <- setdiff(names(op), bad)
  par(op[keep])
}

ensure_device <- function(w = 11, h = 9, pdf_file = "windowed_correlation.pdf") {
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

window_cor <- function(age, x, y, window_width = 50, step = 1,
                       min_n = 15, method = "pearson") {
  
  centers <- seq(min(age, na.rm = TRUE), max(age, na.rm = TRUE), by = step)
  
  out <- data.frame(
    age_center = centers,
    age_min = centers - window_width / 2,
    age_max = centers + window_width / 2,
    n = NA_integer_,
    r = NA_real_,
    p = NA_real_
  )
  
  for (i in seq_along(centers)) {
    keep <- is.finite(age) & is.finite(x) & is.finite(y) &
      age >= out$age_min[i] & age <= out$age_max[i]
    
    out$n[i] <- sum(keep)
    
    if (sum(keep) >= min_n &&
        length(unique(x[keep])) > 2 &&
        length(unique(y[keep])) > 2) {
      
      ct <- suppressWarnings(cor.test(x[keep], y[keep], method = method))
      out$r[i] <- unname(ct$estimate)
      out$p[i] <- ct$p.value
    }
  }
  
  out
}

## =====================  LOAD / PREP DATA  =====================

co2_obj <- get_d13CO2(inv.out)

ages_Ma <- co2_obj$ages
d13CO2_draws <- co2_obj$draws
d13CO2_med <- apply(d13CO2_draws, 1, median, na.rm = TRUE)

## --- Plant d13Corg ---
iso23 <- read_excel("Phan/PhanData/ISOORG23.xlsx", sheet = 1)
iso16 <- read_excel("Phan/PhanData/ISOORG16.xlsx", sheet = 1)

c_age23 <- pick_col(iso23, "1 Myr bin")
c_val23 <- pick_col(iso23, "d13Cp")

d23 <- data.frame(
  age = as_num(iso23[[c_age23]]),
  d13C = as_num(iso23[[c_val23]])
)
d23 <- d23[is.finite(rowSums(d23)) & d23$age >= 0 & d23$age <= 540, , drop = FALSE]

c_age16 <- pick_col(iso16, "Age (Ma)")
c_val16 <- pick_col(iso16, "d13Corganic")

d16 <- data.frame(
  age = as_num(iso16[[c_age16]]),
  d13C = as_num(iso16[[c_val16]])
)
d16 <- d16[is.finite(rowSums(d16)) & d16$age > 250 & d16$age <= 540, , drop = FALSE]

dorg <- rbind(d23, d16)
if (!nrow(dorg)) stop("No valid d13Corg rows in the ISOORG files.")

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

## --- Δorg-atm ---
d13CO2_at_org <- approx(
  x = ages_Ma,
  y = d13CO2_med,
  xout = org_age,
  rule = 2
)$y

Dorg_atm <- d13CO2_at_org - org_mu

## --- Foster et al. CO2 ---
fos <- read_excel("Phan/PhanData/Fosteretal17.xlsx", sheet = 1)

c_age <- pick_col(fos, "age")
c_mu <- pick_col(fos, "co2")
c_u68 <- pick_col(fos, "up68")
c_l68 <- pick_col(fos, "lw68")
c_u95 <- pick_col(fos, "up95")
c_l95 <- pick_col(fos, "lw95")

co2 <- data.frame(
  age = as_num(fos[[c_age]]),
  mu = as_num(fos[[c_mu]]),
  u68 = as_num(fos[[c_u68]]),
  l68 = as_num(fos[[c_l68]]),
  u95 = as_num(fos[[c_u95]]),
  l95 = as_num(fos[[c_l95]])
)

co2 <- co2[is.finite(rowSums(co2)) & co2$age >= 0 & co2$age <= 540, , drop = FALSE]

co2$l95 <- clamp(co2$l95, 100, 3500)
co2$u95 <- clamp(co2$u95, 100, 3500)
co2$l68 <- clamp(co2$l68, 100, 3500)
co2$u68 <- clamp(co2$u68, 100, 3500)
co2$mu <- clamp(co2$mu, 100, 3500)

co2 <- co2[order(co2$age), , drop = FALSE]

co2_at_org <- approx(
  x = co2$age,
  y = co2$mu,
  xout = org_age,
  rule = 2
)$y

log_co2_at_org <- log10(co2_at_org)

## =====================  WINDOWED CORRELATION  =====================

window_width <- 50
window_step <- 1
min_n <- 15

cor_linear <- window_cor(
  age = org_age,
  x = Dorg_atm,
  y = co2_at_org,
  window_width = window_width,
  step = window_step,
  min_n = min_n,
  method = "pearson"
)

cor_log <- window_cor(
  age = org_age,
  x = Dorg_atm,
  y = log_co2_at_org,
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

write.csv(cor_out, "Dorg_atm_CO2_windowed_correlation.csv", row.names = FALSE)

## =====================  PLOTTING  =====================

ensure_device(11, 9)

op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)

layout(matrix(1:4, ncol = 1),
       heights = c(2.8, 2.8, 3.2, 1.0))

par(xaxs = "i", yaxs = "i")

xlim <- c(358, 0)
xticks <- pretty(c(0, 358))

## ===== (1) Δorg-atm =====
par(mar = c(0.8, 5.2, 2.2, 5.2))

yl1 <- range(Dorg_atm, na.rm = TRUE)
yl1 <- yl1 + c(-1, 1)

plot(NA, xlim = xlim, ylim = yl1,
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

grid(nx = NA, ny = NULL)

points(org_age, Dorg_atm,
       pch = 16,
       cex = 0.45,
       col = adjustcolor("darkred", 0.55))

lines(org_age, Dorg_atm,
      col = "darkred",
      lwd = 2)

axis(3, at = xticks, labels = xticks)
axis(2)

mtext(expression(Delta[org-atm]~"\u2030"),
      side = 2, line = 3)

mtext(expression(delta^13*C[CO[2]] - delta^13*C[org]),
      side = 3, line = 0.4, cex = 0.9)

## ===== (2) CO2 =====
par(mar = c(0.6, 5.2, 0.6, 5.2))

yl2 <- c(100, 3500)

plot(NA, xlim = xlim, ylim = yl2,
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

grid(nx = NA, ny = NULL)

x <- co2$age

polygon(c(x, rev(x)),
        c(co2$u95, rev(co2$l95)),
        col = adjustcolor("grey70", 0.6),
        border = NA)

polygon(c(x, rev(x)),
        c(co2$u68, rev(co2$l68)),
        col = adjustcolor("grey50", 0.6),
        border = NA)

lines(co2$age, co2$mu, lwd = 2)

axis(4)
mtext(expression(CO[2]~"(ppm)"), side = 4, line = 3)

## ===== (3) Windowed correlation =====
par(mar = c(0.8, 5.2, 0.8, 5.2))

plot(NA, xlim = xlim, ylim = c(-1, 1),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

grid(nx = NA, ny = NULL)
abline(h = 0, lty = 2)

sig <- is.finite(cor_out$p_Dorg_logCO2) & cor_out$p_Dorg_logCO2 < 0.05

if (any(sig)) {
  points(cor_out$age_center[sig],
         cor_out$r_Dorg_logCO2[sig],
         pch = 16,
         cex = 0.65,
         col = adjustcolor("black", 0.7))
}

lines(cor_out$age_center,
      cor_out$r_Dorg_CO2,
      lwd = 1.8,
      lty = 2,
      col = "grey40")

lines(cor_out$age_center,
      cor_out$r_Dorg_logCO2,
      lwd = 2.5,
      col = "black")

axis(2)

mtext(expression("Windowed correlation, " * r),
      side = 2, line = 3)

legend("bottomleft",
       legend = c(expression(Delta[org-atm]~"vs CO"[2]),
                  expression(Delta[org-atm]~"vs log"[10]*"(CO"[2]*")"),
                  "p < 0.05 for log CO2"),
       lwd = c(1.8, 2.5, NA),
       lty = c(2, 1, NA),
       pch = c(NA, NA, 16),
       col = c("grey40", "black", "black"),
       bty = "n",
       cex = 0.9)

mtext(paste0(window_width, " Myr moving window; minimum n = ", min_n),
      side = 3, line = -1.2, adj = 0.98, cex = 0.85)

## ===== (4) Bottom axis =====
par(mar = c(3.6, 5.2, 0.2, 5.2))

plot(NA, xlim = xlim, ylim = c(0, 1),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "",
     bty = "n")

axis(1, at = xticks, labels = xticks)
mtext("Age (Ma)", side = 1, line = 2.4)

## Save on macOS if using Quartz
if (Sys.info()[["sysname"]] == "Darwin") {
  quartz.save("Dorg_atm_CO2_windowed_correlation.pdf",
              type = "pdf",
              device = dev.cur())
}