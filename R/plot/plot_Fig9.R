# Figure 9: GMST-BWT coupling sensitivity test

load_run <- function(file){
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  e
}

get_draws <- function(run, parameter){
  draws <- as.matrix(run$inv.out$BUGSoutput$sims.list[[parameter]])
  if (ncol(draws) == length(run$ages)) return(draws)
  if (nrow(draws) == length(run$ages)) return(t(draws))
  stop(parameter, " dimensions do not match ages")
}

qband <- function(draws){
  q <- apply(draws, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  list(q025 = q[1, ], median = q[2, ], q975 = q[3, ])
}

if (!exists("model.output.root", inherits = FALSE)) model.output.root <- "output/model_runs/final_archiveblock_3M"
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"

main <- load_run(file.path(model.output.root, "inv_out_main.rda"))
coupled <- load_run(file.path(model.output.root, "inv_out_coupled.rda"))

if (!isTRUE(all.equal(main$ages, coupled$ages))) stop("The main and coupled runs do not use the same age grid")

ages <- main$ages
keep <- which(ages >= 0 & ages <= 540000)
age.Ma <- ages[keep]/1000

main.q <- lapply(c("d13CO2", "BWT", "GMST"), function(x) qband(get_draws(main, x)))
coupled.q <- lapply(c("d13CO2", "BWT", "GMST"), function(x) qband(get_draws(coupled, x)))

main.fill <- adjustcolor("grey55", 0.30)
main.line <- "black"
coupled.fill <- adjustcolor("#56B4E9", 0.32)
coupled.line <- "#0072B2"

plot_panel <- function(main.band, coupled.band, ylab, panel, show.x = FALSE, show.legend = FALSE){
  ylim <- range(main.band$q025[keep], main.band$q975[keep],
                coupled.band$q025[keep], coupled.band$q975[keep], na.rm = TRUE)
  plot(NA, xlim = c(540, 0), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab,
       xaxs = "i", yaxs = "i", las = 1)
  abline(h = pretty(ylim), col = "grey85", lty = 3)
  polygon(c(age.Ma, rev(age.Ma)),
          c(main.band$q025[keep], rev(main.band$q975[keep])),
          col = main.fill, border = NA)
  polygon(c(age.Ma, rev(age.Ma)),
          c(coupled.band$q025[keep], rev(coupled.band$q975[keep])),
          col = coupled.fill, border = NA)
  lines(age.Ma, main.band$median[keep], col = main.line, lwd = 1.4)
  lines(age.Ma, coupled.band$median[keep], col = coupled.line, lwd = 1.4)
  box()
  usr <- par("usr")
  text(usr[1] + 0.02*(usr[2]-usr[1]), usr[4] - 0.06*diff(usr[3:4]),
       panel, adj = c(0, 1), font = 2)
  if (show.x) axis(1, at = seq(0, 500, 50))
  if (show.legend){
    legend("topright", legend = c("main", "coupled"),
           col = c(main.line, coupled.line), lwd = 2,
           fill = c(main.fill, coupled.fill), border = NA, bty = "n")
  }
}

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure9_coupled_sensitivity.pdf"),
          width = 5.735804, height = 5.2)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1), mar = c(0.8, 4.6, 0.7, 0.8), oma = c(3.3, 0, 0, 0))

plot_panel(main.q[[1]], coupled.q[[1]],
           expression(delta^13*C[CO[2]]~("‰ VPDB")), "a.", show.legend = TRUE)
plot_panel(main.q[[2]], coupled.q[[2]],
           expression(BWT[global]~(degree*C)), "b.")
plot_panel(main.q[[3]], coupled.q[[3]],
           expression(GMST~(degree*C)), "c.", show.x = TRUE)
mtext("Age (Ma)", side = 1, outer = TRUE, line = 2.0)

par(op)
dev.off()
