# Figure S1: GMST-BWT coupling sensitivity test

load_run <- function(file){
  run <- readRDS(file)
  if (!all(c("ages", "d13CO2", "BWT", "GMST") %in% names(run))) {
    stop("Run summary is incomplete: ", file)
  }
  run
}

if (!exists("model.output.root", inherits = FALSE)) {
  source("R/model/d13CO2_RunPaths.R", local = TRUE)
  selected.run <- if (exists("model.run", inherits = FALSE)) model.run else NULL
  model.output.root <- d13CO2_model_run_dir(selected.run, must.exist = TRUE)
}
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"

main <- load_run(file.path(model.output.root, "posterior_summary_main.rds"))
coupled <- load_run(file.path(model.output.root, "posterior_summary_coupled.rds"))

if (!isTRUE(all.equal(main$ages, coupled$ages))) stop("The main and coupled runs do not use the same age grid")

ages <- main$ages
keep <- which(ages >= 0 & ages <= 540000)
age.Ma <- ages[keep]/1000

main.q <- main[c("d13CO2", "BWT", "GMST")]
coupled.q <- coupled[c("d13CO2", "BWT", "GMST")]

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
  lines(age.Ma, main.band$med[keep], col = main.line, lwd = 1.4)
  lines(age.Ma, coupled.band$med[keep], col = coupled.line, lwd = 1.4, lty = ifelse(panel == "A.", 2, 1))
  box()
  usr <- par("usr")
  text(usr[1] + 0.02*(usr[2]-usr[1]), if (panel == "C.") usr[3] + 0.14*diff(usr[3:4]) else usr[4] - 0.06*diff(usr[3:4]),
       panel, adj = c(0, 1), font = 2)
  if (show.x) axis(1, at = seq(0, 500, 50))
  if (show.legend){
    legend("topright", legend = c("main", "coupled"),
           col = c(main.line, coupled.line), lwd = 2, lty = c(1, 2),
           fill = c(main.fill, coupled.fill), border = NA, bty = "n")
  }
}

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "FigureS1_coupled_sensitivity.pdf"),
          width = 5.735804, height = 5.8)

op <- par(no.readonly = TRUE)
layout(matrix(1:4, ncol = 1), heights = c(1, 1, 1, lcm(1.672896)))
par(mar = c(0.8, 4.6, 0.7, 0.8), oma = c(0, 0, 0, 0))

plot_panel(main.q[[1]], coupled.q[[1]],
           expression(delta^13*C[a]~("‰ VPDB")), "A.", show.legend = TRUE)
plot_panel(main.q[[2]], coupled.q[[2]],
           expression(BWT[global]~(degree*C)), "B.")
par(mar = c(0.8, 4.6, 0.7, 0.8))
plot_panel(main.q[[3]], coupled.q[[3]],
           expression(GMST~(degree*C)), "C.")

## GSA timescale strip
par(mar = c(3.3, 4.6, 0.1, 0.8))

gsa <- data.frame(
  from = c(0.00, 2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40),
  to = c(2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40, 540.00),
  lab = c("Q", "Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Є"),
  col = c("#FFF2A6", "#FFD37F", "#FFAB70", "#7FC64E", "#34B2C9", "#812B92", "#F04028", "#67A599", "#CB8C37", "#B3E1B6", "#009270", "#7FA056")
)

plot(NA, xlim = c(540, 0), ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")

for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
}
box()
for (i in seq_len(nrow(gsa))) {
  text((gsa$from[i] + gsa$to[i]) / 2, 0.5, gsa$lab[i], cex = 0.8, adj = c(ifelse(i == 1, 1, 0.5), 0.5), xpd = NA)
}

axis(1, at = seq(0, 500, 50))
mtext("Age (Ma)", side = 1, line = 2.0)

par(op)
dev.off()
