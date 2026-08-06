# Figure 8: paleogeographic-model sensitivity test

load_run <- function(file) {
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  if (!all(c("inv.out", "ages", "run.metadata") %in% ls(e))) stop("Run bundle is incomplete: ", file)
  e
}
get_draws <- function(run, parameter) {
  mat <- as.matrix(run$inv.out$BUGSoutput$sims.list[[parameter]])
  if (ncol(mat) == length(run$ages)) return(mat)
  if (nrow(mat) == length(run$ages)) return(t(mat))
  stop(parameter, " dimensions do not match ages")
}
qband <- function(mat) {
  q <- apply(mat, 2, quantile, probs = c(0.025, 0.975, 0.5), na.rm = TRUE)
  list(q025 = q[1, ], q975 = q[2, ], med = q[3, ])
}

if (!exists("model.output.root", inherits = FALSE)) model.output.root <- "output/model_runs/final_archiveblock_3M"
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"
runs <- list(
  Sc16 = load_run(file.path(model.output.root, "inv_out_main.rda")),
  TC17 = load_run(file.path(model.output.root, "inv_out_plate_torsvik2017.rda")),
  MU22 = load_run(file.path(model.output.root, "inv_out_plate_merdith2021.rda")),
  CAO24 = load_run(file.path(model.output.root, "inv_out_plate_cao2024.rda"))
)

if (!all(vapply(runs[-1], function(x) isTRUE(all.equal(x$ages, runs$Sc16$ages)), logical(1)))) {
  stop("Plate-model runs do not use the same age grid")
}

qb <- lapply(runs, function(x) qband(get_draws(x, "d13CO2")))
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
