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
runs$MU22 <- load_run(file.path(model.output.root, "posterior_summary_plate_merdith2021.rds"))
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
col_MU22_fill <- adjustcolor("#56B4E9", alpha.f = 0.35); col_MU22_line <- "#0072B2"
col_CAO24_fill <- adjustcolor("#009E73", alpha.f = 0.35); col_CAO24_line <- "#009E73"

fill <- c(Sc16 = col_Sc16_fill, TC17 = col_TC17_fill, MU22 = col_MU22_fill, CAO24 = col_CAO24_fill)
line <- c(Sc16 = col_Sc16_line, TC17 = col_TC17_line, MU22 = col_MU22_line, CAO24 = col_CAO24_line)
ylab_expr <- expression(delta^13*C[CO[2]]~"(" * "‰" * ")")

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure8_paleogeography_sensitivity.pdf"), width = 5.696350, height = 2.848176)
par(mar = c(3.2, 4.4, 0.5, 0.8), mgp = c(2.1, 0.6, 0), xaxs = "i", yaxs = "i")
plot(x_ma, qb$Sc16$med[idx], type = "n", xlab = "age (Ma)", ylab = ylab_expr,
     xlim = rev(range(x_ma)), ylim = ylim_all)
grid(nx = NA, ny = NULL)

for (nm in names(qb)) {
  polygon(c(x_ma, rev(x_ma)), c(qb[[nm]]$q025[idx], rev(qb[[nm]]$q975[idx])),
          col = fill[[nm]], border = NA)
  lines(x_ma, qb[[nm]]$med[idx], col = line[[nm]], lwd = 2)
}

legend("topright",
       legend = c("Scotese (2016)", "Torsvik & Cocks (2017)",
                  "Müller et al. (2022)", "Cao et al. (2024)"),
       col = line, lwd = 2, bty = "n", seg.len = 3, cex = 0.8)

invisible(dev.off())
