

####################################################################################################
####################################################################################################
# Script to generate Figure 6 in manuscript: late Cenozoic plant wax, atmospheric d13C, 
# fractionation, CO2, and GMST plot

# Plant wax d13C, atmospheric d13C from model posterior, leaf-atmosphere fractionation,
# atmospheric CO2 reconstruction, and GMST anomaly reconstruction

# Assumes objects already exist in the environment:
# inv.out
# ages

if (!exists("model.output.root", inherits = FALSE)) {
  source("R/model/d13CO2_RunPaths.R", local = TRUE)
  selected.run <- if (exists("model.run", inherits = FALSE)) model.run else NULL
  model.output.root <- d13CO2_model_run_dir(selected.run, must.exist = TRUE)
}
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"
load(file.path(model.output.root, "inv_out_cenozoic.rda"))

####################################################################################################
####################################################################################################

library(readxl)

age.max <- 31.5

Delta_C3 <- -24.7
Delta_wsC3 <- -21.8
Delta_C4 <- -15.1

####################################################################################################
####################################################################################################


####################################################################################################
# Plant wax d13C data
####################################################################################################
dat <- read_excel("data/raw/lateCenozoic_plantwax_d13C.xlsx", sheet = "my_comp", col_names = TRUE)

dat$age <- as.numeric(dat$age)
dat$d13C <- as.numeric(dat$d13C)
dat$category <- as.character(dat$category)
dat$reference <- as.character(dat$reference)

dat <- dat[complete.cases(dat[, c("age", "d13C", "category")]), ]
dat <- dat[dat$age >= 0 & dat$age <= age.max, ]

# Tipple et al. (2010) Figure 4 benthic-based atmospheric d13C reconstruction
tipple.raw <- read.csv("data/processed/Tipple2010_Fig4_benthic_d13CO2.csv",
                       stringsAsFactors = FALSE)
tipple <- data.frame(age = as.numeric(tipple.raw$age_Ma),
                     d13Ca = as.numeric(tipple.raw$d13CO2_permil))
tipple <- tipple[is.finite(tipple$age) & is.finite(tipple$d13Ca), ]
tipple <- tipple[order(tipple$age), ]


####################################################################################################
# Atmospheric d13C from inv.out
####################################################################################################
d13CO2_draws <- as.matrix(inv.out$BUGSoutput$sims.list$d13CO2)

if (ncol(d13CO2_draws) == length(ages)) {
  d13CO2_draws <- t(d13CO2_draws)
}

if (nrow(d13CO2_draws) != length(ages)) {
  stop("Dimensions of d13CO2 draws do not match ages.")
}

ages_Ma <- ages / 1000

atm <- data.frame(
  age = ages_Ma,
  d13Ca_2p5 = apply(d13CO2_draws, 1, quantile, probs = 0.025, na.rm = TRUE),
  d13Ca_25 = apply(d13CO2_draws, 1, quantile, probs = 0.25, na.rm = TRUE),
  d13Ca_50 = apply(d13CO2_draws, 1, median, na.rm = TRUE),
  d13Ca_75 = apply(d13CO2_draws, 1, quantile, probs = 0.75, na.rm = TRUE),
  d13Ca_97p5 = apply(d13CO2_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
)

atm <- atm[complete.cases(atm[, c("age", "d13Ca_2p5", "d13Ca_25", "d13Ca_50", "d13Ca_75", "d13Ca_97p5")]), ]
atm <- atm[atm$age >= 0 & atm$age <= age.max, ]
atm <- atm[order(atm$age), ]


####################################################################################################
# Atmospheric CO2 from CenCO2PIP
####################################################################################################
co2 <- read.csv("data/raw/CenCO2PIP_500kyrCO2.csv", stringsAsFactors = FALSE)

co2$age <- as.numeric(co2[[1]])
co2$CO2_2p5 <- exp(as.numeric(co2[[2]]))
co2$CO2_25 <- exp(as.numeric(co2[[3]]))
co2$CO2_50 <- exp(as.numeric(co2[[4]]))
co2$CO2_75 <- exp(as.numeric(co2[[5]]))
co2$CO2_97p5 <- exp(as.numeric(co2[[6]]))

co2 <- co2[complete.cases(co2[, c("age", "CO2_2p5", "CO2_25", "CO2_50", "CO2_75", "CO2_97p5")]), ]
co2 <- co2[co2$age >= 0 & co2$age <= age.max, ]
co2 <- co2[order(co2$age), ]


####################################################################################################
# GMST anomaly from CenCO2PIP
####################################################################################################
gmst <- read.csv("data/raw/CenCO2PIP_500kyrTemp.csv", stringsAsFactors = FALSE)

gmst$age <- as.numeric(gmst[[1]])
gmst$GMST_2p5 <- as.numeric(gmst[[2]])
gmst$GMST_25 <- as.numeric(gmst[[3]])
gmst$GMST_50 <- as.numeric(gmst[[4]])
gmst$GMST_75 <- as.numeric(gmst[[5]])
gmst$GMST_97p5 <- as.numeric(gmst[[6]])

gmst <- gmst[complete.cases(gmst[, c("age", "GMST_2p5", "GMST_25", "GMST_50", "GMST_75", "GMST_97p5")]), ]
gmst <- gmst[gmst$age >= 0 & gmst$age <= age.max, ]
gmst <- gmst[order(gmst$age), ]


####################################################################################################
# Assign subgroup and depositional environment
####################################################################################################
dat$group <- NA_character_
dat$env <- NA_character_

dat$group[grepl("^Africa", dat$category)] <- "Africa"
dat$group[grepl("^lacustrine", dat$category)] <- "lacustrine"
dat$group[grepl("^DSDP94_C31", dat$category)] <- "DSDP94_C31"
dat$group[grepl("^DSDP94_C29", dat$category)] <- "DSDP94_C29"
dat$group[grepl("^SouthChinaSea", dat$category)] <- "SouthChinaSea"
dat$group[grepl("^marine", dat$category)] <- "marine"
dat$group[grepl("^terrestrial", dat$category)] <- "terrestrial"

dat$env[dat$group %in% c("marine", "DSDP94_C31", "DSDP94_C29", "SouthChinaSea")] <- "marine"
dat$env[dat$group %in% c("terrestrial", "Africa")] <- "terrestrial"
dat$env[dat$group %in% c("lacustrine")] <- "lacustrine"

dat <- dat[!is.na(dat$group) & !is.na(dat$env), ]

group_levels <- c("marine", "lacustrine", "terrestrial", "DSDP94_C31", "DSDP94_C29", "Africa", "SouthChinaSea")
dat$group <- factor(dat$group, levels = group_levels)

env_levels <- c("marine", "terrestrial", "lacustrine")
dat$env <- factor(dat$env, levels = env_levels)


####################################################################################################
# Interpolate atmospheric d13C to plant wax ages
####################################################################################################
dat$d13C_atm <- approx(
  x = atm$age,
  y = atm$d13Ca_50,
  xout = dat$age,
  rule = 1
)$y

dat$Delta <- dat$d13C - dat$d13C_atm
dat$d13C_atm_Tipple <- approx(tipple$age, tipple$d13Ca, xout = dat$age, rule = 1)$y
dat$Delta_Tipple <- dat$d13C - dat$d13C_atm_Tipple


####################################################################################################
# Plotting styles
####################################################################################################
group_pchs <- c(
  marine = 16,
  lacustrine = 17,
  terrestrial = 15,
  DSDP94_C31 = 0,
  DSDP94_C29 = 2,
  Africa = 18,
  SouthChinaSea = 1
)

group_cols_plot <- c(
  marine = "#2C7FB8",
  lacustrine = "#31A354",
  terrestrial = "#B35806",
  DSDP94_C31 = "#2C7FB8",
  DSDP94_C29 = "#2C7FB8",
  Africa = "#B35806",
  SouthChinaSea = "#2C7FB8"
)

group_labels <- c(
  marine = "marine - Peppe et al. (2023) comp",
  lacustrine = "lacustrine - Peppe et al. (2023) comp",
  terrestrial = "terrestrial - Peppe et al. (2023) comp",
  DSDP94_C31 = "DSDP 94 (C31) - Tipple and Pagani (2010)",
  DSDP94_C29 = "DSDP 94 (C29) - Tipple and Pagani (2010)",
  Africa = "African sites - Peppe et al. (2023)",
  SouthChinaSea = "South China Sea sites - Zhou et al. (2017)"
)

cats <- levels(dat$group)
cats <- cats[cats %in% as.character(unique(dat$group))]

major_tick <- 0.03
minor_tick <- 0.015

cex_axis <- 0.9
cex_lab <- 0.88
cex_panel <- 1.22
cex_ref <- 0.78

ylab_a <- expression(delta^13*C[wax]~"(‰ VPDB)")
ylab_b <- expression(delta^13*C[CO[2]]~"(‰ VPDB)")
ylab_c <- expression(Delta[leaf-atm] ~ "(" * "\u2030" ~ VPDB * ")")
ylab_d <- expression(paleo-CO[2] ~ "(ppm)")
ylab_e <- expression("GMST anomaly (" * degree * "C)")


####################################################################################################
# Shared x-axis settings
####################################################################################################
xlim <- c(age.max, 0)
x_major <- seq(0, age.max, by = 5)
x_minor <- seq(0, age.max, by = 1)


####################################################################################################
# Panel label helper
####################################################################################################
panel_lab <- function(lab, xlim, ylim) {
  x <- xlim[1] - diff(range(xlim)) * 0.035
  y <- ylim[1] + diff(range(ylim)) * 0.08
  text(x, y, lab, font = 2, cex = cex_panel, adj = c(0, 0))
}


####################################################################################################
# Plot
####################################################################################################
dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure6.pdf"), width = 5.388889, height = 6.555556)

op <- par(no.readonly = TRUE)

layout(matrix(1:5, ncol = 1), heights = c(1, 1, 1, 1, 1))

par(
  mgp = c(2.4, 0.7, 0),
  tcl = 0,
  xaxs = "i",
  yaxs = "i",
  cex.axis = cex_axis,
  cex.lab = cex_lab
)


####################################################################################################
# Panel a: Plant wax d13C
####################################################################################################
par(mar = c(0.4, 5.2, 0.2, 5.2))

ylim_a <- range(dat$d13C, na.rm = TRUE)
y_major_a <- pretty(ylim_a)
y_minor_a <- seq(floor(min(ylim_a, na.rm = TRUE)),
                 ceiling(max(ylim_a, na.rm = TRUE)),
                 by = 1)

plot(
  NA,
  NA,
  xlim = xlim,
  ylim = ylim_a,
  xlab = "",
  ylab = "",
  axes = FALSE
)

for (i in seq_along(cats)) {
  ii <- dat$group == cats[i]
  points(
    dat$age[ii],
    dat$d13C[ii],
    pch = group_pchs[cats[i]],
    col = group_cols_plot[cats[i]],
    bg = group_cols_plot[cats[i]],
    cex = 1.1
  )
}

axis(1, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(2, at = y_major_a, las = 1, tck = major_tick, cex.axis = cex_axis)
axis(3, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(4, at = y_major_a, labels = FALSE, tck = major_tick, cex.axis = cex_axis)

axis(1, at = x_minor, labels = FALSE, tck = minor_tick)
axis(2, at = y_minor_a, labels = FALSE, tck = minor_tick)
axis(3, at = x_minor, labels = FALSE, tck = minor_tick)
axis(4, at = y_minor_a, labels = FALSE, tck = minor_tick)

box(lwd = 1)

mtext(ylab_a, side = 2, line = 3, cex = cex_lab)

legend(
  "topleft",
  legend = group_labels[cats],
  pch = group_pchs[cats],
  col = group_cols_plot[cats],
  pt.bg = group_cols_plot[cats],
  bty = "n",
  cex = 0.78
)

panel_lab("a.", xlim, ylim_a)


####################################################################################################
# Panel b: d13CCO2 with credible intervals
####################################################################################################
par(mar = c(0.4, 5.2, 0.2, 5.2))

ylim_b <- range(c(atm$d13Ca_2p5, atm$d13Ca_97p5), na.rm = TRUE)
y_major_b <- pretty(ylim_b)
y_minor_b <- seq(floor(min(ylim_b, na.rm = TRUE) / 0.1) * 0.1,
                 ceiling(max(ylim_b, na.rm = TRUE) / 0.1) * 0.1,
                 by = 0.1)

plot(
  NA,
  NA,
  xlim = xlim,
  ylim = ylim_b,
  xlab = "",
  ylab = "",
  axes = FALSE
)

polygon(
  x = c(atm$age, rev(atm$age)),
  y = c(atm$d13Ca_2p5, rev(atm$d13Ca_97p5)),
  col = adjustcolor("dodgerblue", 0.18),
  border = NA
)

polygon(
  x = c(atm$age, rev(atm$age)),
  y = c(atm$d13Ca_25, rev(atm$d13Ca_75)),
  col = adjustcolor("dodgerblue", 0.35),
  border = NA
)

lines(atm$age, atm$d13Ca_50, col = "dodgerblue", lwd = 2)

axis(1, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(2, at = y_major_b, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(3, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(4, at = y_major_b, las = 1, tck = major_tick, cex.axis = cex_axis)

axis(1, at = x_minor, labels = FALSE, tck = minor_tick)
axis(2, at = y_minor_b, labels = FALSE, tck = minor_tick)
axis(3, at = x_minor, labels = FALSE, tck = minor_tick)
axis(4, at = y_minor_b, labels = FALSE, tck = minor_tick)

box(lwd = 1)

mtext(ylab_b, side = 4, line = 3, cex = cex_lab)

panel_lab("b.", xlim, ylim_b)


####################################################################################################
# Panel c: Photosynthetic fractionation
####################################################################################################
par(mar = c(0.4, 5.2, 0.2, 5.2))

ylim_c <- range(c(dat$Delta, dat$Delta_Tipple, Delta_C3, Delta_wsC3, Delta_C4), na.rm = TRUE)
y_major_c <- pretty(ylim_c)
y_minor_c <- seq(floor(min(ylim_c, na.rm = TRUE)),
                 ceiling(max(ylim_c, na.rm = TRUE)),
                 by = 1)

plot(
  NA,
  NA,
  xlim = xlim,
  ylim = ylim_c,
  xlab = "",
  ylab = "",
  axes = FALSE
)

points(dat$age, dat$Delta_Tipple, pch = 1, cex = 0.65,
       col = adjustcolor(gray(0.4), 0.7))

abline(h = Delta_C3, lty = 2, lwd = 1.15, col = gray(0.25))
abline(h = Delta_wsC3, lty = 2, lwd = 1.15, col = gray(0.25))
abline(h = Delta_C4, lty = 2, lwd = 1.15, col = gray(0.25))

text(
  x = age.max - 0.6,
  y = Delta_C3 + diff(range(ylim_c)) * 0.025,
  labels = "C3",
  adj = c(0, 0),
  cex = cex_ref,
  col = gray(0.25)
)

text(
  x = age.max - 0.6,
  y = Delta_wsC3 + diff(range(ylim_c)) * 0.025,
  labels = "wsC3",
  adj = c(0, 0),
  cex = cex_ref,
  col = gray(0.25)
)

text(
  x = age.max - 0.6,
  y = Delta_C4 + diff(range(ylim_c)) * 0.025,
  labels = "C4",
  adj = c(0, 0),
  cex = cex_ref,
  col = gray(0.25)
)

for (i in seq_along(cats)) {
  ii <- dat$group == cats[i] & !is.na(dat$Delta)
  points(
    dat$age[ii],
    dat$Delta[ii],
    pch = group_pchs[cats[i]],
    col = group_cols_plot[cats[i]],
    bg = group_cols_plot[cats[i]],
    cex = 1.1
  )
}

legend("bottomleft", legend = "Tipple et al. (2010) Fig. 4-based",
       col = gray(0.4), pch = 1, pt.cex = 0.75,
       inset = c(0.055, 0.01), bty = "n", cex = 0.7)

axis(1, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(2, at = y_major_c, las = 1, tck = major_tick, cex.axis = cex_axis)
axis(3, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(4, at = y_major_c, labels = FALSE, tck = major_tick, cex.axis = cex_axis)

axis(1, at = x_minor, labels = FALSE, tck = minor_tick)
axis(2, at = y_minor_c, labels = FALSE, tck = minor_tick)
axis(3, at = x_minor, labels = FALSE, tck = minor_tick)
axis(4, at = y_minor_c, labels = FALSE, tck = minor_tick)

box(lwd = 1)

mtext(ylab_c, side = 2, line = 3, cex = cex_lab)

panel_lab("c.", xlim, ylim_c)


####################################################################################################
# Panel d: Atmospheric CO2 with credible intervals
####################################################################################################
par(mar = c(0.4, 5.2, 0.2, 5.2))

ylim_d <- range(c(co2$CO2_2p5, co2$CO2_97p5), na.rm = TRUE)
y_major_d <- pretty(ylim_d)
y_minor_d <- seq(floor(min(ylim_d, na.rm = TRUE) / 50) * 50,
                 ceiling(max(ylim_d, na.rm = TRUE) / 50) * 50,
                 by = 50)

plot(
  NA,
  NA,
  xlim = xlim,
  ylim = ylim_d,
  xlab = "",
  ylab = "",
  axes = FALSE
)

polygon(
  x = c(co2$age, rev(co2$age)),
  y = c(co2$CO2_2p5, rev(co2$CO2_97p5)),
  col = gray(0.85),
  border = NA
)

polygon(
  x = c(co2$age, rev(co2$age)),
  y = c(co2$CO2_25, rev(co2$CO2_75)),
  col = gray(0.65),
  border = NA
)

lines(co2$age, co2$CO2_50, col = gray(0.2), lwd = 2)

axis(1, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(2, at = y_major_d, las = 1, tck = major_tick, cex.axis = cex_axis)
axis(3, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(4, at = y_major_d, labels = FALSE, tck = major_tick, cex.axis = cex_axis)

axis(1, at = x_minor, labels = FALSE, tck = minor_tick)
axis(2, at = y_minor_d, labels = FALSE, tck = minor_tick)
axis(3, at = x_minor, labels = FALSE, tck = minor_tick)
axis(4, at = y_minor_d, labels = FALSE, tck = minor_tick)

box(lwd = 1)

mtext(ylab_d, side = 2, line = 3, cex = cex_lab)

panel_lab("d.", xlim, ylim_d)


####################################################################################################
# Panel e: GMST anomaly with credible intervals
####################################################################################################
par(mar = c(3.0, 5.2, 0.2, 5.2))

ylim_e <- range(c(gmst$GMST_2p5, gmst$GMST_97p5), na.rm = TRUE)
y_major_e <- pretty(ylim_e)
y_minor_e <- seq(floor(min(ylim_e, na.rm = TRUE)),
                 ceiling(max(ylim_e, na.rm = TRUE)),
                 by = 1)

plot(
  NA,
  NA,
  xlim = xlim,
  ylim = ylim_e,
  xlab = "",
  ylab = "",
  axes = FALSE
)

polygon(
  x = c(gmst$age, rev(gmst$age)),
  y = c(gmst$GMST_2p5, rev(gmst$GMST_97p5)),
  col = adjustcolor("firebrick", 0.18),
  border = NA
)

polygon(
  x = c(gmst$age, rev(gmst$age)),
  y = c(gmst$GMST_25, rev(gmst$GMST_75)),
  col = adjustcolor("firebrick", 0.35),
  border = NA
)

lines(gmst$age, gmst$GMST_50, col = "firebrick", lwd = 2)

axis(1, at = x_major, tck = major_tick, cex.axis = cex_axis)
axis(2, at = y_major_e, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(3, at = x_major, labels = FALSE, tck = major_tick, cex.axis = cex_axis)
axis(4, at = y_major_e, las = 1, tck = major_tick, cex.axis = cex_axis)

axis(1, at = x_minor, labels = FALSE, tck = minor_tick)
axis(2, at = y_minor_e, labels = FALSE, tck = minor_tick)
axis(3, at = x_minor, labels = FALSE, tck = minor_tick)
axis(4, at = y_minor_e, labels = FALSE, tck = minor_tick)

box(lwd = 1)

mtext("Age (Ma)", side = 1, line = 2.4, cex = cex_lab)
mtext(ylab_e, side = 4, line = 3, cex = cex_lab)

panel_lab("e.", xlim, ylim_e)

invisible(dev.off())
