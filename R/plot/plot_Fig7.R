# Figure 7: GMST sensitivity test

safe_par_reset <- function(op) {
  bad <- c("cin","cra","csi","cxy","din","page","usr","mfg","plt","pin","fin","fig","omd")
  par(op[setdiff(names(op), bad)])
}
load_run <- function(file) {
  run <- readRDS(file)
  required <- c("ages", "GMST.m", "GMST.sd", "run.metadata", "d13CO2", "GMST")
  if (!all(required %in% names(run))) stop("Run summary is incomplete: ", file)
  run
}

## =====================  LOAD RUN BUNDLES  =====================
if (!exists("model.output.root", inherits = FALSE)) {
  source("R/model/d13CO2_RunPaths.R", local = TRUE)
  selected.run <- if (exists("model.run", inherits = FALSE)) model.run else NULL
  model.output.root <- d13CO2_model_run_dir(selected.run, must.exist = TRUE)
}
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"
phan <- load_run(file.path(model.output.root, "posterior_summary_main.rds"))
li <- load_run(file.path(model.output.root, "posterior_summary_gmst_scotese.rds"))

if (phan$run.metadata$GMST_model != "PhanDA") stop("The main run is not the PhanDA profile")
if (li$run.metadata$GMST_model != "Scotese21") stop("The sensitivity run is not the Scotese21/Li22 profile")
if (!isTRUE(all.equal(phan$ages, li$ages))) stop("The two runs do not use the same age grid")

ages <- phan$ages
idx <- which(ages >= 0 & ages <= 540000)
x_common <- ages[idx]/1000
xlim_all <- rev(range(x_common, na.rm = TRUE))

## =====================  QUANTILES  =====================
d13_Phan <- phan$d13CO2
d13_Li <- li$d13CO2
gmst_Phan <- phan$GMST
gmst_Li <- li$GMST

GMST.low <- phan$GMST.m - 1.96*phan$GMST.sd
GMST.hi <- phan$GMST.m + 1.96*phan$GMST.sd
GMST.low.Li <- li$GMST.m - 1.96*li$GMST.sd
GMST.hi.Li <- li$GMST.m + 1.96*li$GMST.sd

## =====================  PLOTTING  =====================
col_phan_ts_fill <- adjustcolor("#E66101", 0.35)
col_phan_ts_line <- "#E66101"
col_li_ts_fill <- adjustcolor("#5E3C99", 0.35)
col_li_ts_line <- "#5E3C99"
col_phan_prior_fill <- adjustcolor("#A50026", 0.35)
col_phan_prior_line <- "#A50026"
col_li_prior_fill <- adjustcolor("#4575B4", 0.35)
col_li_prior_line <- "#4575B4"

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure7_GMST_sensitivity.pdf"), width = 5.735804, height = 3.639535)
op <- par(no.readonly = TRUE)
par(xaxs = "i", yaxs = "i")
layout(matrix(1:3, ncol = 1), heights = c(3.2, 3.2, 3.0))

panel_lab <- function(lab) {
  usr <- par("usr")
  text(usr[1] - 0.025*diff(range(usr[1:2])), usr[4] - 0.08*diff(usr[3:4]),
       lab, adj = c(0, 1), font = 2)
}

## ----- (1) TOP: d13C_CO2 posteriors -----
par(mar = c(0.5, 5.0, 2.2, 1.8))
ylab_d13 <- expression(delta^13*C[CO[2]]~"(" * "‰" * ")")
ylim_top <- range(d13_Phan$q025[idx], d13_Phan$q975[idx],
                  d13_Li$q025[idx], d13_Li$q975[idx], na.rm = TRUE)
plot(NA, xlim = xlim_all, ylim = ylim_top, xlab = "", ylab = ylab_d13, xaxt = "n", yaxt = "s")
grid(nx = NA, ny = NULL)
polygon(c(x_common, rev(x_common)), c(d13_Phan$q025[idx], rev(d13_Phan$q975[idx])),
        col = col_phan_ts_fill, border = NA)
lines(x_common, d13_Phan$med[idx], col = col_phan_ts_line, lwd = 1.2)
polygon(c(x_common, rev(x_common)), c(d13_Li$q025[idx], rev(d13_Li$q975[idx])),
        col = col_li_ts_fill, border = NA)
lines(x_common, d13_Li$med[idx], col = col_li_ts_line, lwd = 1.2)
legend("topright", legend = c("PhanDA GMST", "Scotese21/Li22 GMST"),
       fill = c(col_phan_ts_fill, col_li_ts_fill), border = NA, bty = "n")
panel_lab("a.")

## ----- (2) MIDDLE: GMST posteriors -----
par(mar = c(0.5, 5.0, 0.8, 1.8))
ylab_gmst <- expression(paste("GMST (", degree, "C)"))
ymax_mid <- max(gmst_Phan$q975[idx], gmst_Li$q975[idx], na.rm = TRUE)
plot(NA, xlim = xlim_all, ylim = c(0, ymax_mid), xlab = "", ylab = ylab_gmst,
     xaxt = "n", yaxt = "s")
grid(nx = NA, ny = NULL)
polygon(c(x_common, rev(x_common)), c(gmst_Phan$q025[idx], rev(gmst_Phan$q975[idx])),
        col = col_phan_ts_fill, border = NA)
lines(x_common, gmst_Phan$med[idx], col = col_phan_ts_line, lwd = 1.2)
polygon(c(x_common, rev(x_common)), c(gmst_Li$q025[idx], rev(gmst_Li$q975[idx])),
        col = col_li_ts_fill, border = NA)
lines(x_common, gmst_Li$med[idx], col = col_li_ts_line, lwd = 1.2)
legend("topright", legend = c("PhanDA GMST", "Scotese21/Li22 GMST"),
       fill = c(col_phan_ts_fill, col_li_ts_fill), border = NA, bty = "n")
panel_lab("b.")

## ----- (3) BOTTOM: GMST priors -----
par(mar = c(3.6, 5.0, 0.8, 1.8))
ylim_bot <- range(GMST.low[idx], GMST.hi[idx], GMST.low.Li[idx], GMST.hi.Li[idx], na.rm = TRUE)
plot(x_common, phan$GMST.m[idx], type = "n", xlab = "age (Ma)", ylab = ylab_gmst,
     xlim = xlim_all, ylim = ylim_bot)
grid(nx = NA, ny = NULL)
polygon(c(x_common, rev(x_common)), c(GMST.low[idx], rev(GMST.hi[idx])),
        col = col_phan_prior_fill, border = NA)
lines(x_common, phan$GMST.m[idx], col = col_phan_prior_line, lwd = 1.5)
polygon(c(x_common, rev(x_common)), c(GMST.low.Li[idx], rev(GMST.hi.Li[idx])),
        col = col_li_prior_fill, border = NA)
lines(x_common, li$GMST.m[idx], col = col_li_prior_line, lwd = 1.5)
legend("topleft", legend = c("PhanDA GMST prior", "Scotese21/Li22 GMST prior"),
       fill = c(col_phan_prior_fill, col_li_prior_fill), border = NA, bty = "n", inset = 0.01)
panel_lab("c.")

invisible(dev.off())
