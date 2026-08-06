
####################################################################################################
####################################################################################################
# Script to generate Figure 4 in manuscript: NSB prior vs. posterior density panels
#
# Assumes object exists:
# inv.out

if (!exists("model.output.root", inherits = FALSE)) model.output.root <- "output/model_runs/final_archiveblock_3M"
if (!exists("figure.output.root", inherits = FALSE)) figure.output.root <- "output/figures"
load(file.path(model.output.root, "inv_out_main.rda"))

####################################################################################################
####################################################################################################


####################################################################################################
# Hyperparameters from JAGS model
####################################################################################################

bf.nsb.m <- 0
bf.nsb.sd <- 0.25

pf.nsb.m <- 0
pf.nsb.sd <- 0.25

brach.nsb.m <- 0
brach.nsb.sd <- 1

bivalve.nsb.m <- 0
bivalve.nsb.sd <- 1

amm.nsb.m <- 0
amm.nsb.sd <- 1

bel.nsb.m <- 0
bel.nsb.sd <- 1

micrite.nsb.m <- 0
micrite.nsb.sd <- 0.25

bulk.nsb.m <- 0
bulk.nsb.sd <- 0.5

bulk_sr.nsb.m <- 0
bulk_sr.nsb.sd <- 1

bulk_marg.nsb.m <- 0
bulk_marg.nsb.sd <- 0.75

n_prior_draws <- 1e5


####################################################################################################
# Helper functions
####################################################################################################

get_site_draws <- function(inv.out, stem) {
  nm <- paste0(stem, ".nsb_site")
  draws <- inv.out$BUGSoutput$sims.list[[nm]]
  
  if (is.null(draws)) {
    stop("No draws found for ", nm)
  }
  
  as.matrix(draws)
}

rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  out <- numeric(0)
  
  while (length(out) < n) {
    x <- rnorm(n, mean = mean, sd = sd)
    x <- x[x >= lower & x <= upper]
    out <- c(out, x)
  }
  
  out[seq_len(n)]
}

sample_A_terms <- function(n = n_prior_draws) {
  asd <- rtnorm(n, mean = 1, sd = 0.2, lower = 0)
  bpump <- rtnorm(n, mean = 1.2, sd = 0.4, lower = 0)
  remin <- rtnorm(n, mean = 0.6, sd = 0.3, lower = 0)
  
  list(
    Asurf = asd,
    Abf = asd + bpump + remin
  )
}

sample_site_prior <- function(stem, m, sd, A_terms, n = n_prior_draws) {
  eta_site <- rnorm(n, mean = m, sd = sd)
  
  if (stem == "bf") {
    eta_site + A_terms$Abf
  } else {
    eta_site + A_terms$Asurf
  }
}

avg_site_density <- function(draws, n_grid = 512) {
  site_rngs <- apply(draws, 2, range, na.rm = TRUE)
  xlo <- min(site_rngs[1, ], na.rm = TRUE)
  xhi <- max(site_rngs[2, ], na.rm = TRUE)
  
  if (!is.finite(xlo) || !is.finite(xhi) || xlo == xhi) {
    m <- mean(draws, na.rm = TRUE)
    s <- sd(as.numeric(draws), na.rm = TRUE)
    xlo <- m - 4 * s
    xhi <- m + 4 * s
  }
  
  xg <- seq(xlo, xhi, length.out = n_grid)
  yg <- matrix(NA, nrow = n_grid, ncol = ncol(draws))
  
  for (j in seq_len(ncol(draws))) {
    dj <- density(draws[, j], na.rm = TRUE)
    yg[, j] <- approx(dj$x, dj$y, xg, rule = 2)$y
  }
  
  list(x = xg, y = rowMeans(yg, na.rm = TRUE))
}


####################################################################################################
# Family definitions
####################################################################################################

families <- c("bf", "pf", "brach", "bivalve", "amm", "bel", "bulk", "micrite", "bulk_sr", "bulk_marg")

fam_label <- c(
  bf = "a. benthic foraminifera",
  pf = "b. planktic foraminifera",
  brach = "c. brachiopod",
  bivalve = "d. bivalve",
  amm = "e. ammonite",
  bel = "f. belemnite",
  bulk = "g. bulk carbonate (open ocean)",
  micrite = "h. micrite (open ocean)",
  bulk_sr = "i. bulk carbonate (semi-restricted)",
  bulk_marg = "j. bulk carbonate (marginal sea)"
)

fam_mu <- c(
  bf = bf.nsb.m,
  pf = pf.nsb.m,
  brach = brach.nsb.m,
  bivalve = bivalve.nsb.m,
  amm = amm.nsb.m,
  bel = bel.nsb.m,
  bulk = bulk.nsb.m,
  micrite = micrite.nsb.m,
  bulk_sr = bulk_sr.nsb.m,
  bulk_marg = bulk_marg.nsb.m
)

fam_sd <- c(
  bf = bf.nsb.sd,
  pf = pf.nsb.sd,
  brach = brach.nsb.sd,
  bivalve = bivalve.nsb.sd,
  amm = amm.nsb.sd,
  bel = bel.nsb.sd,
  bulk = bulk.nsb.sd,
  micrite = micrite.nsb.sd,
  bulk_sr = bulk_sr.nsb.sd,
  bulk_marg = bulk_marg.nsb.sd
)


####################################################################################################
# Prior draws
####################################################################################################

set.seed(1)

A_terms <- sample_A_terms(n_prior_draws)

prior_draws <- vector("list", length(families))
names(prior_draws) <- families

for (stem in families) {
  prior_draws[[stem]] <- sample_site_prior(
    stem = stem,
    m = fam_mu[[stem]],
    sd = fam_sd[[stem]],
    A_terms = A_terms,
    n = n_prior_draws
  )
}


####################################################################################################
# Plot
####################################################################################################

dir.create(figure.output.root, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure.output.root, "Figure4.pdf"), width = 5.976744, height = 4.524064)

op <- par(no.readonly = TRUE)

layout(matrix(1:10, nrow = 5, byrow = TRUE))
par(xaxs = "i", yaxs = "i", mar = c(1.45, 2.65, 0.35, 0.35),
    oma = c(0, 2.0, 0, 0), mgp = c(1.35, 0.35, 0),
    tcl = -0.16, cex.axis = 0.68, cex.lab = 0.72)

col_prior <- adjustcolor("purple4", 0.45)
col_post <- adjustcolor("coral1", 0.5)
col_pmed <- "purple4"
col_mmed <- "darkred"

xlab_expr <- expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)"))
ylab_str <- "Probability density"

for (i in seq_along(families)) {
  stem <- families[i]
  
  draws_mat <- get_site_draws(inv.out, stem)
  dpost <- avg_site_density(draws_mat)
  
  xprior <- prior_draws[[stem]]
  dprior <- density(xprior)
  
  xlo <- min(dprior$x, dpost$x, na.rm = TRUE)
  xhi <- max(dprior$x, dpost$x, na.rm = TRUE)
  ymax <- 1.05 * max(dprior$y, dpost$y, na.rm = TRUE)
  
  plot(NA, xlim = c(xlo, xhi), ylim = c(0, ymax),
       xlab = "", ylab = "", axes = FALSE)

  axis(1)
  axis(2, labels = i %% 2 == 1, las = 1)
  box()
  if (i %in% c(9, 10)) mtext(xlab_expr, side = 1, line = 1.0, cex = 0.72)
  
  grid(nx = NA, ny = NULL)
  
  polygon(c(dprior$x, rev(dprior$x)),
          c(dprior$y, rep(0, length(dprior$y))),
          col = col_prior, border = NA)
  
  polygon(c(dpost$x, rev(dpost$x)),
          c(dpost$y, rep(0, length(dpost$y))),
          col = col_post, border = NA)
  
  pmed <- median(xprior, na.rm = TRUE)
  pmed_y <- approx(dprior$x, dprior$y, pmed, rule = 2)$y
  segments(pmed, 0, pmed, pmed_y, col = col_pmed, lty = "dotted")
  
  post_med <- median(as.numeric(draws_mat), na.rm = TRUE)
  post_med_y <- approx(dpost$x, dpost$y, post_med, rule = 2)$y
  segments(post_med, 0, post_med, post_med_y, col = col_mmed, lty = "dashed")
  
  if (i == 2) {
    legend("topright", legend = c("prior", "posterior"),
           fill = c(col_prior, col_post), border = NA, bty = "n", cex = 0.62)
  }
  
  text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]),
       par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
       labels = fam_label[[stem]],
       adj = c(0, 1),
       cex = 0.72,
       font = 2)
}

mtext(ylab_str, side = 2, outer = TRUE, line = 0.65, cex = 0.75)

invisible(dev.off())
