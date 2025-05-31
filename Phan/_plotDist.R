#  Assumed objects:
#    prox.in       : your proxy data.frame
#    flattened     : data.frame with columns ai, site.index, ages
#    inv.out       : your JAGS fit (with tempC[...] and tempC_bot[...] saved)
#    BWT           : data.frame with BWT$age, BWT$BWT, BWT$BWT_2sd
#    sites         : data.frame with sites$site (name), sites$site.index

# Impute NAs in GMST_PhanDA by nearest‐neighbor 
for (col in c("GMST_PhanDA","GMST_PhanDA_hi","GMST_PhanDA_lo")) {
  ord <- order(prox.in$age)
  df  <- data.frame(age = prox.in$age[ord], val = prox.in[[col]][ord])
  df  <- df[!is.na(df$val), ]
  df  <- df[!duplicated(df$age), ]
  prox.in[[col]] <- approx(x    = df$age,
                           y    = df$val,
                           xout = prox.in$age,
                           method = "constant",
                           rule = 2)$y
}

PhanDA_sd <- ((prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo) +
                (prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA)) / 2

# First ai per site → flattened index
first_ai <- tapply(prox.in$ai, prox.in$site.index, min)
flat_idx <- sapply(names(first_ai), function(si) {
  ai0 <- first_ai[[si]]
  which(flattened$site.index == as.numeric(si) & flattened$ai == ai0)
})
names(flat_idx) <- names(first_ai)

# Extract posteriors from JAGS matrix
post_mat        <- inv.out$BUGSoutput$sims.matrix
post_tempC      <- lapply(flat_idx, function(idx)
  post_mat[, paste0("tempC[", idx, "]")])
post_tempC_bot  <- lapply(flat_idx, function(idx)
  post_mat[, paste0("tempC_bot[", idx, "]")])

# Build prior parms by site
prior_params <- lapply(seq_along(first_ai), function(i) {
  si   <- as.numeric(names(first_ai)[i])
  idxf <- flat_idx[i]
  subset <- (prox.in$site.index == si & prox.in$ai == first_ai[[i]])
  
  # flag bf
  is_bf <- any(prox.in$material[subset] == "bf")
  if (is_bf) {
    age0 <- flattened$ages[idxf]    # <— use flattened$ages!
    mu0  <- approx(BWT$age, BWT$BWT, 
                   xout = age0, rule = 2)$y
    sd0  <- approx(BWT$age, BWT$BWT_2sd/2, 
                   xout = age0, rule = 2)$y
  } else {
    vals <- prox.in[subset, ]
    mu0  <- mean(vals$GMST_PhanDA + vals$temp_offset, na.rm = TRUE)
    sd0  <- mean(sqrt(PhanDA_sd[subset]^2 + vals$temp_offset_sd^2),
                 na.rm = TRUE)
  }
  c(mean = mu0, sd = sd0)
})
names(prior_params) <- names(first_ai)

# Prep for paged plotting
site_nums <- as.numeric(names(first_ai))
site_labs <- sapply(site_nums, function(si)
  unique(prox.in$site[prox.in$site.index == si]))
n_sites   <- length(site_nums)
per_page  <- 18
n_pages   <- ceiling(n_sites / per_page)

prior_col <- grDevices::adjustcolor("coral1",     alpha.f = 0.5)
post_col  <- grDevices::adjustcolor("dodgerblue", alpha.f = 0.5)

# Loop over pages
for (page in seq_len(n_pages)) {
  dev.new(width = 12, height = 8)
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 6),
      mar   = c(3, 1.5, 1, 0.5),
      oma   = c(0, 0, 0, 0),
      mgp   = c(2, 0.3, 0))
  
  i1 <- (page - 1)*per_page + 1
  i2 <- min(page*per_page, n_sites)
  
  for (ii in seq(i1, i2)) {
    si    <- site_nums[ii]
    lbl   <- site_labs[ii]
    params<- prior_params[[ as.character(si) ]]
    mu0   <- params["mean"]
    sd0   <- params["sd"]
    
    # skip if prior bad
    if (!is.finite(mu0) || !is.finite(sd0)) {
      plot.new()
      mtext(lbl, side = 1, line = 2, cex = 0.7)
      next
    }
    
    # prior density
    x1 <- seq(mu0 - 4*sd0, mu0 + 4*sd0, length.out = 512)
    y1 <- dnorm(x1, mu0, sd0)
    
    # pick posterior draws
    vals <- if (any(prox.in$material[prox.in$site.index==si &
                                     prox.in$ai==first_ai[[as.character(si)]]] == "bf")) {
      post_tempC_bot[[ as.character(si) ]]
    } else {
      post_tempC   [[ as.character(si) ]]
    }
    vals <- na.omit(vals)
    
    # posterior density?
    has_dd <- length(vals) >= 2
    if (has_dd) {
      d2   <- density(vals)
      xlim <- range(x1, d2$x)
      ylim <- c(0, max(max(y1), max(d2$y)))
    } else {
      xlim <- range(x1)
      ylim <- c(0, max(y1))
    }
    
    # plot
    plot(NA, xlim = xlim, ylim = ylim,
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    polygon(c(x1, rev(x1)), c(y1, rep(0, length(y1))),
            col = prior_col, border = NA)
    if (has_dd) {
      polygon(c(d2$x, rev(d2$x)), c(d2$y, rep(0, length(d2$y))),
              col = post_col, border = NA)
    }
    
    # x‐axis
    atx <- pretty(xlim, n = 3)
    axis(1, at = atx, labels = atx,
         cex.axis = 0.6, tck = -0.02, line = 0.7)
    
    # site name
    mtext(lbl, side = 1, line = 2, cex = 0.7)
    
    # legend on first panel
    if (ii == i1) {
      legend("topright",
             legend = c("Prior","Posterior"),
             fill   = c(prior_col, post_col),
             bty    = "n", cex = 0.7, inset = 0.02)
    }
  }
  
  # fill blanks
  drawn <- i2 - i1 + 1
  if (drawn < per_page) {
    for (j in seq(drawn+1, per_page)) plot.new()
  }
  
  par(old.par)
}
