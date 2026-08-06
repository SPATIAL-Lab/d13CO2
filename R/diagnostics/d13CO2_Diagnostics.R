# Save convergence diagnostics for one model run

split.chains <- function(x) {
  n <- nrow(x)
  h <- floor(n/2)
  cbind(x[seq_len(h), , drop = FALSE], x[n-h+seq_len(h), , drop = FALSE])
}

basic.rhat <- function(x) {
  n <- nrow(x)
  W <- mean(apply(x, 2, var))
  B <- n*var(colMeans(x))
  sqrt(((n-1)*W/n+B/n)/W)
}

rank.rhat <- function(x) {
  x <- split.chains(x)
  normalize <- function(z) {
    r <- rank(z, ties.method = "average")
    qnorm((r-3/8)/(length(r)+1/4))
  }
  z <- matrix(normalize(as.vector(x)), nrow = nrow(x), ncol = ncol(x))
  folded <- abs(x-median(x))
  z.folded <- matrix(normalize(as.vector(folded)),
                     nrow = nrow(x), ncol = ncol(x))
  max(basic.rhat(z), basic.rhat(z.folded))
}

autocorrelation.ess <- function(x) {
  chains <- lapply(seq_len(ncol(x)), function(j) coda::mcmc(x[, j]))
  suppressWarnings(as.numeric(coda::effectiveSize(coda::mcmc.list(chains))))
}

variable.names <- dimnames(inv.out$BUGSoutput$sims.array)[[3]]

extract.array <- function(prefix) {
  columns <- grep(paste0("^", prefix, "\\["), variable.names)
  index <- as.integer(gsub(".*\\[|\\]", "", variable.names[columns]))
  inv.out$BUGSoutput$sims.array[, , columns[order(index)], drop = FALSE]
}

archive.names <- run.metadata$active.archives
scalar.names <- intersect(
  c("d13CO2_sigma", "GMST_sigma", "BWT_sigma", "cc_co2_constant1",
    "cc_co2_coeff1", "asd", "bpump", "remin", "f_co3", "f_carbacid",
    "Abf", "Asurf", "d13CO2_level", "BWT_GMST_beta",
    paste0(archive.names, ".nsb_mean"), paste0(archive.names, ".nsb_tau"),
    paste0(archive.names, ".archive_level")),
  variable.names
)

scalar.diagnostics <- do.call(rbind, lapply(scalar.names, function(name) {
  x <- inv.out$BUGSoutput$sims.array[, , name]
  data.frame(parameter = name, median = median(x),
             lower = quantile(x, 0.025), upper = quantile(x, 0.975),
             Rhat = rank.rhat(x), ESS = autocorrelation.ess(x))
}))

keep <- ages >= 0 & ages <= run.metadata$age.max
series <- lapply(c("d13CO2", "GMST", "BWT"), extract.array)
names(series) <- c("d13CO2", "GMST", "BWT")
series <- lapply(series, function(x) x[, , keep, drop = FALSE])

diagnose.series <- function(name) {
  x <- series[[name]]
  center <- apply(x, c(1, 2), mean)

  diagnose <- function(form, n.nodes, get.node) {
    do.call(rbind, lapply(seq_len(n.nodes), function(i) {
      z <- get.node(i)
      data.frame(parameter = name, form = form, time.index = i,
                 Rhat = rank.rhat(z), ESS = autocorrelation.ess(z))
    }))
  }

  rbind(
    diagnose("raw", dim(x)[3], function(i) x[, , i]),
    diagnose("centered", dim(x)[3], function(i) x[, , i]-center),
    diagnose("first_difference", dim(x)[3]-1L,
             function(i) x[, , i+1L]-x[, , i])
  )
}

temporal.diagnostics <- do.call(rbind,
                                lapply(names(series), diagnose.series))

groups <- split(temporal.diagnostics,
                interaction(temporal.diagnostics$parameter,
                            temporal.diagnostics$form, drop = TRUE))
diagnostic.summary <- do.call(rbind, lapply(groups, function(x) {
  data.frame(parameter = x$parameter[1], form = x$form[1],
             max_Rhat = max(x$Rhat), median_Rhat = median(x$Rhat),
             min_ESS = min(x$ESS), median_ESS = median(x$ESS),
             n_Rhat_over_1.03 = sum(x$Rhat > 1.03),
             n_ESS_under_100 = sum(x$ESS < 100),
             n_ESS_under_400 = sum(x$ESS < 400))
}))
rownames(diagnostic.summary) <- NULL

write.csv(scalar.diagnostics,
          file.path(output.root, paste0("scalar_diagnostics_", run.profile, ".csv")),
          row.names = FALSE)
write.csv(temporal.diagnostics,
          file.path(output.root, paste0("diagnostics_", run.profile, ".csv")),
          row.names = FALSE)
write.csv(diagnostic.summary,
          file.path(output.root, paste0("diagnostics_summary_", run.profile, ".csv")),
          row.names = FALSE)
