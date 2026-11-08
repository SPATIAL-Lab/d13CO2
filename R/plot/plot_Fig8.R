# Figure 8: paleogeographic-model sensitivity test

load_run <- function(file) {
  run <- readRDS(file)
  if (!all(c("ages", "run.metadata", "d13CO2") %in% names(run))) {
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
runs <- list()
runs$Sc16 <- load_run(file.path(model.output.root, "posterior_summary_main.rds"))
runs$TC17 <- load_run(file.path(model.output.root, "posterior_summary_plate_torsvik2017.rds"))
runs$ME21 <- load_run(file.path(model.output.root, "posterior_summary_plate_merdith2021.rds"))
runs$CAO24 <- load_run(file.path(model.output.root, "posterior_summary_plate_cao2024.rds"))

if (!all(vapply(runs[-1], function(x) isTRUE(all.equal(x$ages, runs$Sc16$ages)), logical(1)))) {
  stop("Plate-model runs do not use the same age grid")
}

qb <- lapply(runs, `[[`, "d13CO2")
ages <- runs$Sc16$ages
idx <- which(ages >= 0 & ages <= 540000)
x_ma <- ages[idx]/1000

ylim_all <- range(unlist(lapply(qb, function(x) c(x$q025[idx], x$q975[idx]))), na.rm = TRUE)

col_Sc16_fill <- adjustcolor("#999999", alpha.f = 0.35); col_Sc16_line <- "#000000"
col_TC17_fill <- adjustcolor("#E69F00", alpha.f = 0.35); col_TC17_line <- "#E69F00"
col_ME21_fill <- adjustcolor("#56B4E9", alpha.f = 0.35); col_ME21_line <- "#0072B2"
col_CAO24_fill <- adjustcolor("#009E73", alpha.f = 0.35); col_CAO24_line <- "#009E73"

fill <- c(Sc16 = col_Sc16_fill, TC17 = col_TC17_fill, ME21 = col_ME21_fill, CAO24 = col_CAO24_fill)
line <- c(Sc16 = col_Sc16_line, TC17 = col_TC17_line, ME21 = col_ME21_line, CAO24 = col_CAO24_line)
ylab_expr <- expression(delta^13*C[a]~"(" * "‰" * ")")

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure8_paleogeography_sensitivity.pdf"), width = 8.544525, height = 3.448176)
layout(matrix(1:2, ncol = 1), heights = c(5, lcm(2.140314381)))
par(mar = c(0.5, 4.4, 0.5, 0.8), mgp = c(2.1, 0.6, 0), xaxs = "i", yaxs = "i")
plot(x_ma, qb$Sc16$med[idx], type = "n", xlab = "", ylab = ylab_expr,
     xlim = rev(range(x_ma)), ylim = ylim_all, xaxt = "n")
grid(nx = NA, ny = NULL)

for (nm in names(qb)) {
  polygon(c(x_ma, rev(x_ma)), c(qb[[nm]]$q025[idx], rev(qb[[nm]]$q975[idx])),
          col = fill[[nm]], border = NA)
}
for (nm in c(names(qb)[-1], "Sc16")) {
  lines(x_ma, qb[[nm]]$med[idx], col = line[[nm]], lwd = 2)
}

legend("bottomright",
       legend = c("Scotese (2016)", "Torsvik & Cocks (2017)",
                  "Merdith et al. (2021)", "Cao et al. (2024)"),
       col = line, lwd = 2, bty = "n", seg.len = 3, cex = 0.8)

## GSA timescale strip
par(mar = c(2.6, 4.4, 0.1, 0.8))

gsa <- data.frame(
  from = c(0.00, 2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40),
  to = c(2.58, 23.03, 66.00, 145.00, 201.30, 251.90, 298.90, 358.90, 419.20, 443.80, 485.40, 540.00),
  lab = c("Q", "Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Є"),
  col = c("#FFF2A6", "#FFD37F", "#FFAB70", "#7FC64E", "#34B2C9", "#812B92", "#F04028", "#67A599", "#CB8C37", "#B3E1B6", "#009270", "#7FA056")
)

plot(NA, xlim = rev(range(x_ma)), ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")

for (i in seq_len(nrow(gsa))) {
  rect(gsa$to[i], 0, gsa$from[i], 1, col = gsa$col[i], border = NA)
}
box()
for (i in seq_len(nrow(gsa))) {
  text((gsa$from[i] + gsa$to[i]) / 2, 0.5, gsa$lab[i], cex = 0.65, adj = c(ifelse(i == 1, 1, 0.5), 0.5), xpd = NA)
}

axis(1)
mtext("Age (Ma)", side = 1, line = 1.7)

invisible(dev.off())
