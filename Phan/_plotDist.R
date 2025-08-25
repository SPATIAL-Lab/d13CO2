

#  Load objects:
load("Phan/chpc_output/prox.in_PhanDA.rda")
load("Phan/chpc_output/flattened_PhanDA.rda")
load("Phan/chpc_output/inv.out_PhanDA.rda")
load("Phan/chpc_output/BWT_PhanDA.rda")
load("Phan/chpc_output/sites_PhanDA.rda")


#  Assumed objects:
#    prox.in: proxy data.frame
#    flattened: data.frame with columns ai, site.index, ages
#    inv.out: JAGS mcmc output object
#    BWT: data.frame with BWT$age, BWT$BWT, BWT$BWT_2sd
#    sites: data.frame with sites$site (name), sites$site.index

##### Aggregate prior vs posterior across all ai per site #####

if (!exists("toff_sd_uniform_bot")) toff_sd_uniform_bot <- 1

# Categories
surf_cats <- c("bulk", "bulk open water", "bulk open ocean", "micrite open ocean","bulk semi restricted","bulk marginal sea", "bulk marginal sea restricting up section")
bf_cat    <- "bf"

# Ensure GMST PhanDA sd exists on prox.in
if (!("GMST_PhanDA_hi" %in% names(prox.in) && "GMST_PhanDA_lo" %in% names(prox.in)))
  stop("prox.in must have GMST_PhanDA_hi/lo.")
PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                (prox.in$GMST_PhanDA    - prox.in$GMST_PhanDA_lo)) / 2

# Interpolate row-wise (per flattened row) inputs on the exact ages 
gmst_row   <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = flattened$ages, rule = 2)$y
gmstsd_row <- approx(prox.in$age, PhanDA_sd,            xout = flattened$ages, rule = 2)$y
# Choose the temp offset fields that match driver 
toff_row   <- if ("temp_offset_interp" %in% names(flattened)) flattened$temp_offset_interp else flattened$temp_offset_PhanDA_interp
toffsd_row <- if ("temp_offset_sd_interp" %in% names(flattened)) flattened$temp_offset_sd_interp else stop("Need flattened$temp_offset_sd_interp")

# Bottom water (per flattened row)
bwt_row <- approx(BWT$age, BWT$BWT, xout = flattened$ages, rule = 2)$y
bwtsd_row <- approx(BWT$age, BWT$BWT_2sd/2, xout = flattened$ages, rule = 2)$y

# Posterior draws 
sl <- inv.out$BUGSoutput$sims.list
if (is.null(sl$tempC) || is.null(sl$tempC_bot)) {
  stop("tempC / tempC_bot not found in sims.list; include them in `parameters.to.save`.")
}
# Dimensions: [iterations, nrow(flattened)]
stopifnot(ncol(sl$tempC) == nrow(flattened), ncol(sl$tempC_bot) == nrow(flattened))

# Build per-site index sets for SURF (bulk+micrite) and BOT (bf)
site_nums <- sort(unique(flattened$site.index))
site_labs <- sapply(site_nums, function(si) unique(prox.in$site[prox.in$site.index == si])[1])

# For each site, find all flattened rows that correspond to SURF or BOT observations at any ai
idx_by_site_surf <- lapply(site_nums, function(si) {
  ai_surf <- unique(prox.in$ai[prox.in$site.index == si & prox.in$category %in% surf_cats])
  if (length(ai_surf) == 0) return(integer(0))
  idx_site <- which(flattened$site.index == si)
  match(ai_surf, flattened$ai[idx_site]) |> (\(ii) idx_site[ii])() |> na.omit() |> as.integer()
})
names(idx_by_site_surf) <- site_nums

idx_by_site_bot <- lapply(site_nums, function(si) {
  ai_bot <- unique(prox.in$ai[prox.in$site.index == si & prox.in$category == bf_cat])
  if (length(ai_bot) == 0) return(integer(0))
  idx_site <- which(flattened$site.index == si)
  match(ai_bot, flattened$ai[idx_site]) |> (\(ii) idx_site[ii])() |> na.omit() |> as.integer()
})
names(idx_by_site_bot) <- site_nums

# Mixture-of-Normals prior density helper: average component densities
mixnorm_density <- function(x, mu, sd) {
  if (length(mu) == 0) return(rep(0, length(x)))
  # numeric guard
  sd[sd <= 1e-8] <- 1e-8
  rowMeans(sapply(seq_along(mu), function(k) dnorm(x, mean = mu[k], sd = sd[k])))
}

# Colors
prior_col_surf <- grDevices::adjustcolor("coral1",     alpha.f = 0.55)
post_col_surf  <- grDevices::adjustcolor("dodgerblue", alpha.f = 0.55)
prior_col_bot  <- grDevices::adjustcolor("orange2",  alpha.f = 0.55)
post_col_bot   <- grDevices::adjustcolor("royalblue3",  alpha.f = 0.55)

# Layout
per_page <- 18
n_pages  <- ceiling(length(site_nums) / per_page)

for (page in seq_len(n_pages)) {
  dev.new(width = 12, height = 8)
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 6),
      mar   = c(3, 1.5, 1, 0.5),
      oma   = c(0, 0, 0, 0),
      mgp   = c(2, 0.3, 0))
  
  i1 <- (page - 1) * per_page + 1
  i2 <- min(page * per_page, length(site_nums))
  
  for (ii in i1:i2) {
    si  <- site_nums[ii]
    lbl <- site_labs[ii]
    
    idx_surf <- idx_by_site_surf[[as.character(si)]]
    idx_bot  <- idx_by_site_bot [[as.character(si)]]
    
    # SURF priors: mu = GMST + toff; sd = sqrt(gmstsd^2 + toffsd^2)
    mu_surf <- gmst_row[idx_surf] + toff_row[idx_surf]
    sd_surf <- sqrt(pmax(gmstsd_row[idx_surf], 0)^2 + pmax(toffsd_row[idx_surf], 0)^2)
    
    # BOT priors: mu = BWT; sd = sqrt(bwtsd^2 + toff_bot_sd^2)
    mu_bot <- bwt_row[idx_bot]
    sd_bot <- sqrt(pmax(bwtsd_row[idx_bot], 0)^2 + toff_sd_uniform_bot^2)
    
    # SURF posterior draws pooled across all ai for this site
    post_surf_draws <- if (length(idx_surf)) as.vector(sl$tempC[, idx_surf, drop = FALSE]) else numeric(0)
    # BOT posterior draws pooled
    post_bot_draws  <- if (length(idx_bot))  as.vector(sl$tempC_bot[, idx_bot, drop = FALSE]) else numeric(0)
    
    # Build x-limits from all available info
    xmin <- +Inf; xmax <- -Inf
    if (length(mu_surf)) { xmin <- min(xmin, min(mu_surf - 4*sd_surf, na.rm=TRUE)); xmax <- max(xmax, max(mu_surf + 4*sd_surf, na.rm=TRUE)) }
    if (length(mu_bot))  { xmin <- min(xmin, min(mu_bot  - 4*sd_bot,  na.rm=TRUE)); xmax <- max(xmax, max(mu_bot  + 4*sd_bot,  na.rm=TRUE)) }
    if (length(post_surf_draws)) { xmin <- min(xmin, min(post_surf_draws, na.rm=TRUE)); xmax <- max(xmax, max(post_surf_draws, na.rm=TRUE)) }
    if (length(post_bot_draws))  { xmin <- min(xmin, min(post_bot_draws,  na.rm=TRUE)); xmax <- max(xmax, max(post_bot_draws,  na.rm=TRUE)) }
    if (!is.finite(xmin) || !is.finite(xmax)) { xmin <- -5; xmax <- 25 }  # fallback sensible range (°C)
    
    x <- seq(xmin, xmax, length.out = 512)
    
    # Densities
    y_prior_surf <- if (length(mu_surf)) mixnorm_density(x, mu_surf, sd_surf) else rep(0, length(x))
    y_prior_bot  <- if (length(mu_bot))  mixnorm_density(x, mu_bot,  sd_bot)  else rep(0, length(x))
    y_post_surf  <- if (length(post_surf_draws) >= 2) density(post_surf_draws, from = xmin, to = xmax, n = 512)$y else rep(0, length(x))
    y_post_bot   <- if (length(post_bot_draws)  >= 2) density(post_bot_draws,  from = xmin, to = xmax, n = 512)$y else rep(0, length(x))
    
    ymax <- max(y_prior_surf, y_prior_bot, y_post_surf, y_post_bot, 1e-9)
    
    # Panel
    plot(NA, xlim = c(xmin, xmax), ylim = c(0, ymax),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    # SURF
    if (any(y_prior_surf > 0)) polygon(c(x, rev(x)), c(y_prior_surf, rep(0, length(x))), col = prior_col_surf, border = NA)
    if (any(y_post_surf  > 0)) polygon(c(x, rev(x)), c(y_post_surf,  rep(0, length(x))), col = post_col_surf,  border = NA)
    # BOT
    if (any(y_prior_bot  > 0)) polygon(c(x, rev(x)), c(y_prior_bot,  rep(0, length(x))), col = prior_col_bot,  border = NA)
    if (any(y_post_bot   > 0)) polygon(c(x, rev(x)), c(y_post_bot,   rep(0, length(x))), col = post_col_bot,   border = NA)
    
    # Axis + label
    atx <- pretty(c(xmin, xmax), n = 3)
    axis(1, at = atx, labels = atx, cex.axis = 0.6, tck = -0.02, line = 0.7)
    mtext(lbl, side = 1, line = 2, cex = 0.7)
    
    # Legend once per page (first panel)
    if (ii == i1) {
      legend("topright",
             legend = c("Prior (surface)", "Posterior (surface)",
                        "Prior (benthic)", "Posterior (benthic)"),
             fill   = c(prior_col_surf,    post_col_surf,
                        prior_col_bot,     post_col_bot),
             bty    = "n", cex = 0.7, inset = 0.02)
    }
  }
  
  # Fill blanks if last page partly full
  drawn <- i2 - i1 + 1
  if (drawn < per_page) for (j in seq_len(per_page - drawn)) plot.new()
  
  par(old.par)
}
