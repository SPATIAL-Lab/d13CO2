
# non-secular bias of archival d13C from high-fidelity (true) value (per mille), standard deviation w/ fixed mean = 0 


# Bulk carbonate - benthic forams

bf.nsb.m  <- 0
bf.nsb.sd <- 0.25
shap <- 1e2
rat  <- 2e3
n    <- 1e5

xprior <- rnorm(n, rnorm(n, 0, bf.nsb.sd), rgamma(n, shap, rat))
d      <- density(xprior)
xpost  <- inv.out$BUGSoutput$sims.list$bf.nsb
dpost  <- density(xpost)
ymax <- max(d$y, dpost$y) * 1.05

plot(NA, xlim = c(bf.nsb.sd*-2.5, bf.nsb.sd*2.5), ylim = c(0, ymax),
     xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")

polygon(c(d$x, rev(d$x)),           c(d$y,    rep(0, length(d$y))),
        col = adjustcolor("slateblue", 0.5), border = NA)
polygon(c(dpost$x, rev(dpost$x)),   c(dpost$y, rep(0, length(dpost$y))),
        col = adjustcolor("coral1",  0.5),    border = NA)

pmd <- median(xprior)
pmd_y <- approx(d$x, d$y, pmd, rule = 2)$y
segments(pmd, 0, pmd, pmd_y, col = "navyblue", lty = "dotted")

post_md <- median(xpost)
post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
segments(post_md, 0, post_md, post_md_y, col = "darkred", lty = "dashed")

legend("topright",
       legend = c("foram NSB prior", "foram NSB post"),
       fill   = c(adjustcolor("slateblue", 0.5), adjustcolor("coral1", 0.5)),
       border = NA, bty = "n")



# Bulk carbonate - open ocean

bulk.nsb.m  <- 0
bulk.nsb.sd <- 0.5
shap <- 1e2
rat  <- 2e3
n    <- 1e5

xprior <- rnorm(n, rnorm(n, 0, bulk.nsb.sd), rgamma(n, shap, rat))
d      <- density(xprior)
xpost  <- inv.out$BUGSoutput$sims.list$bulk.nsb
dpost  <- density(xpost)
ymax <- max(d$y, dpost$y) * 1.05

plot(NA, xlim = c(bulk.nsb.sd*-2, bulk.nsb.sd*2), ylim = c(0, ymax),
     xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")

polygon(c(d$x, rev(d$x)),           c(d$y,    rep(0, length(d$y))),
        col = adjustcolor("slateblue", 0.5), border = NA)
polygon(c(dpost$x, rev(dpost$x)),   c(dpost$y, rep(0, length(dpost$y))),
        col = adjustcolor("coral1",  0.5),    border = NA)

pmd <- median(xprior)
pmd_y <- approx(d$x, d$y, pmd, rule = 2)$y
segments(pmd, 0, pmd, pmd_y, col = "navyblue", lty = "dotted")

post_md <- median(xpost)
post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
segments(post_md, 0, post_md, post_md_y, col = "darkred", lty = "dashed")

legend("topright",
       legend = c("bulk open ocean NSB prior", "bulk open ocean NSB post"),
       fill   = c(adjustcolor("slateblue", 0.5), adjustcolor("coral1", 0.5)),
       border = NA, bty = "n")



# Bulk carbonate - micrite

micrite.nsb.m  <- 0
micrite.nsb.sd <- 0.25
shap <- 1e2
rat  <- 2e3
n    <- 1e5

xprior <- rnorm(n, rnorm(n, 0, micrite.nsb.sd), rgamma(n, shap, rat))
d      <- density(xprior)
xpost  <- inv.out$BUGSoutput$sims.list$micrite.nsb
dpost  <- density(xpost)
ymax <- max(d$y, dpost$y) * 1.05

plot(NA, xlim = c(micrite.nsb.sd*-2.5, micrite.nsb.sd*2.5), ylim = c(0, ymax),
     xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")

polygon(c(d$x, rev(d$x)),           c(d$y,    rep(0, length(d$y))),
        col = adjustcolor("slateblue", 0.5), border = NA)
polygon(c(dpost$x, rev(dpost$x)),   c(dpost$y, rep(0, length(dpost$y))),
        col = adjustcolor("coral1",  0.5),    border = NA)

pmd <- median(xprior)
pmd_y <- approx(d$x, d$y, pmd, rule = 2)$y
segments(pmd, 0, pmd, pmd_y, col = "navyblue", lty = "dotted")

post_md <- median(xpost)
post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
segments(post_md, 0, post_md, post_md_y, col = "darkred", lty = "dashed")

legend("topright",
       legend = c("micrite NSB prior", "micrite NSB post"),
       fill   = c(adjustcolor("slateblue", 0.5), adjustcolor("coral1", 0.5)),
       border = NA, bty = "n")



# Bulk carbonate - semi-restricted sea

bulk_sr.nsb.m  <- 0
bulk_sr.nsb.sd <- 1
shap <- 1e2
rat  <- 2e3
n    <- 1e5

xprior <- rnorm(n, rnorm(n, 0, bulk_sr.nsb.sd), rgamma(n, shap, rat))
d      <- density(xprior)
xpost  <- inv.out$BUGSoutput$sims.list$bulk_sr.nsb
dpost  <- density(xpost)
ymax <- max(d$y, dpost$y) * 1.05

plot(NA, xlim = c(bulk_sr.nsb.sd*-2, bulk_sr.nsb.sd*2), ylim = c(0, ymax),
     xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")

polygon(c(d$x, rev(d$x)),           c(d$y,    rep(0, length(d$y))),
        col = adjustcolor("slateblue", 0.5), border = NA)
polygon(c(dpost$x, rev(dpost$x)),   c(dpost$y, rep(0, length(dpost$y))),
        col = adjustcolor("coral1",  0.5),    border = NA)

pmd <- median(xprior)
pmd_y <- approx(d$x, d$y, pmd, rule = 2)$y
segments(pmd, 0, pmd, pmd_y, col = "navyblue", lty = "dotted")

post_md <- median(xpost)
post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
segments(post_md, 0, post_md, post_md_y, col = "darkred", lty = "dashed")

legend("topright",
       legend = c("bulk semi-restricted NSB prior", "bulk semi-restricted NSB post"),
       fill   = c(adjustcolor("slateblue", 0.5), adjustcolor("coral1", 0.5)),
       border = NA, bty = "n")



# Bulk carbonate - marginal sea

bulk_marg.nsb.m  <- 0
bulk_marg.nsb.sd <- 0.75
shap <- 1e2
rat  <- 2e3
n    <- 1e5

xprior <- rnorm(n, rnorm(n, 0, bulk_marg.nsb.sd), rgamma(n, shap, rat))
d      <- density(xprior)
xpost  <- inv.out$BUGSoutput$sims.list$bulk_marg.nsb
dpost  <- density(xpost)
ymax <- max(d$y, dpost$y) * 1.05

plot(NA, xlim = c(bulk_marg.nsb.sd*-2, bulk_marg.nsb.sd*2), ylim = c(0, ymax),
     xlab = expression(paste(epsilon["NSB"], " (", "\u2030", " VPDB)")),
     ylab = "Probability Density")

polygon(c(d$x, rev(d$x)),           c(d$y,    rep(0, length(d$y))),
        col = adjustcolor("slateblue", 0.5), border = NA)
polygon(c(dpost$x, rev(dpost$x)),   c(dpost$y, rep(0, length(dpost$y))),
        col = adjustcolor("coral1",  0.5),    border = NA)

pmd <- median(xprior)
pmd_y <- approx(d$x, d$y, pmd, rule = 2)$y
segments(pmd, 0, pmd, pmd_y, col = "navyblue", lty = "dotted")

post_md <- median(xpost)
post_md_y <- approx(dpost$x, dpost$y, post_md, rule = 2)$y
segments(post_md, 0, post_md, post_md_y, col = "darkred", lty = "dashed")

legend("topright",
       legend = c("bulk marginal NSB prior", "bulk marginal NSB post"),
       fill   = c(adjustcolor("slateblue", 0.5), adjustcolor("coral1", 0.5)),
       border = NA, bty = "n")



