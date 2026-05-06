
# Open a device reliably; fallback to PDF if headless
ensure_device <- function(w = 11, h = 9, pdf_file = "figure_nsb.pdf") {
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

## =====================  MULTIPANEL NSB DENSITIES: PRIOR vs POSTERIOR  =====================

## ---- Hyperparameters (from JAGS model) ----
bf.nsb.m <- 0; bf.nsb.sd <- 0.25
pf.nsb.m <- 0; pf.nsb.sd <- 0.25
brach.nsb.m <- 0; brach.nsb.sd <- 1
bivalve.nsb.m <- 0; bivalve.nsb.sd <- 1
amm.nsb.m <- 0; amm.nsb.sd <- 1
bel.nsb.m <- 0; bel.nsb.sd <- 1
micrite.nsb.m <- 0; micrite.nsb.sd <- 0.25
bulk.nsb.m <- 0; bulk.nsb.sd <- 0.5
bulk_sr.nsb.m <- 0; bulk_sr.nsb.sd <- 1
bulk_marg.nsb.m <- 0; bulk_marg.nsb.sd <- 0.75

HYP_SHAPE <- 1e3
HYP_RATE  <- 1e-3

## ---- Helpers ----
get_site_draws <- function(inv.out, stem) {
  S <- inv.out$BUGSoutput
  nm <- paste0(stem, ".nsb_site")
  if (!is.null(S$sims.list[[nm]])) {
    M <- as.matrix(S$sims.list[[nm]])
    if (nrow(M) < ncol(M)) M <- t(M)
    return(M)
  }
  if (!is.null(S$sims.matrix)) {
    sm  <- S$sims.matrix
    pat <- paste0("^", gsub("\\.", "\\\\.", stem), "\\.nsb_site\\[(\\d+)\\]$")
    keep <- grepl(pat, colnames(sm))
    if (!any(keep)) stop("No draws found for ", stem, ".nsb_site")
    idx <- as.integer(sub(pat, "\\1", colnames(sm)[keep]))
    sm  <- sm[, keep, drop = FALSE]
    sm  <- sm[, order(idx), drop = FALSE]
    return(as.matrix(sm))
  }
  stop("inv.out$BUGSoutput has neither sims.list nor sims.matrix for ", stem, ".nsb_site")
}

sample_site_prior <- function(m, sd, n = 1e5, shape = HYP_SHAPE, rate = HYP_RATE) {
  mu  <- rnorm(n, mean = m, sd = sd)
  tau <- rgamma(n, shape = shape, rate = rate)
  rnorm(n, mean = mu, sd = 1/sqrt(tau))
}

avg_site_density <- function(draws, n_grid = 512) {
  site_rngs <- apply(draws, 2, range, na.rm = TRUE)
  xl <- min(site_rngs[1, ], na.rm = TRUE)
  xr <- max(site_rngs[2, ], na.rm = TRUE)
  if (!is.finite(xl) || !is.finite(xr) || xl == xr) {
    m  <- mean(draws, na.rm = TRUE)
    s  <- sd(as.numeric(draws), na.rm = TRUE)
    xl <- m - 4*s; xr <- m + 4*s
  }
  xg <- seq(xl, xr, length.out = n_grid)
  ys <- matrix(NA_real_, nrow = n_grid, ncol = ncol(draws))
  for (j in seq_len(ncol(draws))) {
    dj <- density(draws[, j], na.rm = TRUE)
    ys[, j] <- approx(dj$x, dj$y, xg, rule = 2)$y
  }
  list(x = xg, y = rowMeans(ys, na.rm = TRUE))
}

## ---- Family definitions ----
families <- c("bf","pf","brach","bivalve","amm","bel","bulk","micrite","bulk_sr","bulk_marg")
fam_label <- c(
  bf = "Benthic foram NSB", pf = "Planktic foram NSB",
  brach = "Brachiopod NSB", bivalve = "Bivalve NSB",
  amm = "Ammonite NSB", bel = "Belemnite NSB",
  bulk = "Bulk (open ocean) NSB", micrite = "Micrite (open ocean) NSB",
  bulk_sr = "Bulk (semi-restricted) NSB", bulk_marg = "Bulk (marginal sea) NSB"
)
fam_mu <- c(
  bf = bf.nsb.m, pf = pf.nsb.m, brach = brach.nsb.m, bivalve = bivalve.nsb.m,
  amm = amm.nsb.m, bel = bel.nsb.m, bulk = bulk.nsb.m, micrite = micrite.nsb.m,
  bulk_sr = bulk_sr.nsb.m, bulk_marg = bulk_marg.nsb.m
)
fam_sd <- c(
  bf = bf.nsb.sd, pf = pf.nsb.sd, brach = brach.nsb.sd, bivalve = bivalve.nsb.sd,
  amm = amm.nsb.sd, bel = bel.nsb.sd, bulk = bulk.nsb.sd, micrite = micrite.nsb.sd,
  bulk_sr = bulk_sr.nsb.sd, bulk_marg = bulk_marg.nsb.sd
)

## ---- Graphics ----
ensure_device(12, 10)
layout(matrix(1:10, nrow = 5, byrow = TRUE))  # 5 rows × 2 columns
op <- par(no.readonly = TRUE)
on.exit(safe_par_reset(op), add = TRUE)
par(xaxs = "i", yaxs = "i", mar = c(3.5, 4.8, 2.6, 1.2))

col_prior <- adjustcolor("slateblue", 0.5)
col_post  <- adjustcolor("coral1",   0.5)
col_pmed  <- "navyblue"
col_mmed  <- "darkred"
xlab_expr <- expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)"))
ylab_str  <- "Probability Density"

set.seed(1)
for (stem in families) {
  draws_mat <- get_site_draws(inv.out, stem)
  dpost     <- avg_site_density(draws_mat)
  
  xprior <- sample_site_prior(m = fam_mu[[stem]], sd = fam_sd[[stem]], n = 1e5)
  dprior <- density(xprior)
  
  xmid <- fam_mu[[stem]]
  xhw  <- 2.5 * fam_sd[[stem]]
  xlo  <- min(xmid - xhw, min(dprior$x, dpost$x, na.rm = TRUE))
  xhi  <- max(xmid + xhw, max(dprior$x, dpost$x, na.rm = TRUE))
  ymax <- 1.05 * max(max(dprior$y, na.rm = TRUE), max(dpost$y, na.rm = TRUE))
  
  plot(NA, xlim = c(xlo, xhi), ylim = c(0, ymax), xlab = xlab_expr, ylab = ylab_str)
  grid(nx = NA, ny = NULL)
  
  polygon(c(dprior$x, rev(dprior$x)), c(dprior$y, rep(0, length(dprior$y))),
          col = col_prior, border = NA)
  polygon(c(dpost$x,  rev(dpost$x)),  c(dpost$y,  rep(0, length(dpost$y))),
          col = col_post, border = NA)
  
  pmd   <- median(xprior, na.rm = TRUE)
  pmd_y <- approx(dprior$x, dprior$y, pmd, rule = 2)$y
  segments(pmd, 0, pmd, pmd_y, col = col_pmed, lty = "dotted")
  
  post_md   <- median(as.numeric(draws_mat), na.rm = TRUE)
  post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
  segments(post_md, 0, post_md, post_md_y, col = col_mmed, lty = "dashed")
  
  # Simplified legend
  legend("topright", legend = c("prior", "posterior"),
         fill = c(col_prior, col_post), border = NA, bty = "n", cex = 0.9)
  
  # ---- Panel title inside top-left corner ----
  text(x = par("usr")[1] + 0.02 * diff(par("usr")[1:2]),
       y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
       labels = fam_label[[stem]], adj = c(0, 1),
       cex = 1.0, font = 2)
}

# Save figure (macOS Quartz)
quartz.save("nsb_prior_vs_posterior_panels.pdf", type = "pdf", device = dev.cur())
