# Driver for the compact CHPC-equivalent model suite

library(R2jags)

if (!exists("run.profile", inherits = FALSE)) run.profile <- "main"
if (!exists("model.run", inherits = FALSE)) model.run <- NULL
source("R/model/d13CO2_RunPaths.R", local = TRUE)
model.run <- d13CO2_model_run_name(model.run)
output.root <- d13CO2_model_run_dir(model.run, create = TRUE)
if (!exists("n.chains", inherits = FALSE)) n.chains <- 4L
if (!exists("n.iter", inherits = FALSE)) n.iter <- 3e6
if (!exists("n.burnin", inherits = FALSE)) n.burnin <- n.iter/2
if (!exists("n.thin", inherits = FALSE)) n.thin <- 250L
if (!exists("run.seed", inherits = FALSE)) run.seed <- 26072600L
if (!exists("save.run.output", inherits = FALSE)) save.run.output <- TRUE
if (!exists("save.full.output", inherits = FALSE)) save.full.output <- FALSE

model.version <- "2026-08-07-compact-archiveblock-uniform-minus-nsb"
if (exists("model.version.override", inherits = FALSE)) model.version <- model.version.override

run.profiles <- list(
  main = list(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv",
              paleogeography = "PALEOMAP", GMST = "PhanDA",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = FALSE),
  gmst_scotese = list(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv",
              paleogeography = "PALEOMAP", GMST = "Scotese21",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = FALSE),
  plate_torsvik2017 = list(file = "data/processed/PhanCompWithTemp_TorsvikCocks2017.csv",
              paleogeography = "TorsvikCocks2017", GMST = "PhanDA",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = FALSE),
  plate_merdith2021 = list(file = "data/processed/PhanCompWithTemp_MERDITH2021.csv",
              paleogeography = "MERDITH2021", GMST = "PhanDA",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = FALSE),
  plate_cao2024 = list(file = "data/processed/PhanCompWithTemp_CAO2024.csv",
              paleogeography = "CAO2024", GMST = "PhanDA",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = FALSE),
  cenozoic = list(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv",
              paleogeography = "PALEOMAP", GMST = "PhanDA",
              age.max = 66000, step.int = 100, n.spinup = 10, coupled = FALSE),
  coupled = list(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv",
              paleogeography = "PALEOMAP", GMST = "PhanDA",
              age.max = 540000, step.int = 1000, n.spinup = 10, coupled = TRUE)
)
if (!run.profile %in% names(run.profiles)) stop("Unknown run profile: ", run.profile)
run.settings <- run.profiles[[run.profile]]

proxy.file <- run.settings$file
bwt.file <- "data/raw/CenozoicBWT.csv"
if (exists("proxy.file.override", inherits = FALSE)) proxy.file <- proxy.file.override
if (exists("bwt.file.override", inherits = FALSE)) bwt.file <- bwt.file.override
age.min <- 0
age.max <- run.settings$age.max
step.int <- run.settings$step.int
n.spinup <- run.settings$n.spinup
GMST_model <- run.settings$GMST
temp_offset_model <- "Li22"
GMST_sd_Scotese21 <- 5
toff_sd_uniform <- 2
toff_sd_uniform_bot <- 1
d13C.analyt.sd <- 0.1
if (!exists("d13CO2_sigma_upper", inherits = FALSE)) d13CO2_sigma_upper <- 0.1

if (run.profile == "cenozoic") {
  model.file <- "R/model/d13CO2_PSM_cenozoic.R"
} else if (run.settings$coupled) {
  model.file <- "R/model/d13CO2_PSM_coupled.R"
} else {
  model.file <- "R/model/d13CO2_PSM.R"
}
if (exists("model.file.override", inherits = FALSE)) model.file <- model.file.override

nsb.priors <- data.frame(
  archive = c("bf", "pf", "brach", "bivalve", "amm", "bel", "micrite", "bulk", "bulk_sr", "bulk_marg"),
  mean = 0,
  sd = c(0.25, 0.25, 1, 1, 1, 1, 0.25, 0.5, 1, 0.75),
  stringsAsFactors = FALSE
)
d13CO2_level_prior_mean <- -6

# Data preparation
####################################################################################################

prox.raw <- as.data.frame(read.csv(proxy.file))
paleogeo.cols <- paste0(c("plng_", "plat_"), run.settings$paleogeography)
climate.cols <- c("MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                  "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA")
required.cols <- c(names(prox.raw)[1:7], paleogeo.cols, climate.cols)
if (!all(required.cols %in% names(prox.raw))) {
  stop("Missing input columns in ", proxy.file, ": ",
       paste(setdiff(required.cols, names(prox.raw)), collapse = ", "))
}
prox.in <- cbind(prox.raw[, 1:7], prox.raw[, paleogeo.cols], prox.raw[, climate.cols],
                 rep(toff_sd_uniform, nrow(prox.raw)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category",
                    "paleolon", "paleolat", "MAT", "GMST_Scotese21", "GMST_PhanDA",
                    "GMST_PhanDA_hi", "GMST_PhanDA_lo", "temp_offset",
                    "temp_offset_PhanDA", "temp_offset_sd")

prox.in$age <- prox.in$age * 1e3
prox.in <- prox.in[prox.in$age >= age.min & prox.in$age <= age.max, ]
prox.in <- transform(prox.in,
                     ai = n.spinup + as.numeric(1 + floor((max(prox.in$age) - prox.in$age) / step.int)))
prox.in <- prox.in[order(prox.in$age, decreasing = TRUE), ]

ages.short <- seq(max(prox.in$age), min(prox.in$age), by = -step.int) - 0.5*step.int
ages <- seq(n.spinup*step.int + max(prox.in$age), min(prox.in$age), by = -step.int) - 0.5*step.int
ai.all <- seq_along(ages)
age.indices <- data.frame(ai = ai.all, age = ages)
n.steps <- length(ages)
dt.scale <- abs(diff(ages))/1000
age.max.spinup <- age.max + step.int*n.spinup

PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                (prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo)) / 2
if (GMST_model == "PhanDA") {
  GMST.m <- approx(prox.in$age, prox.in$GMST_PhanDA, xout = ages, rule = 2, ties = mean)$y
  GMST.sd <- approx(prox.in$age, PhanDA_sd, xout = ages, rule = 2, ties = mean)$y
} else {
  GMST.m <- approx(prox.in$age, prox.in$GMST_Scotese21, xout = ages, rule = 2, ties = mean)$y
  GMST.sd <- rep(GMST_sd_Scotese21, length(ages))
}

prox.in <- transform(prox.in, site.index = as.numeric(factor(site, ordered = is.ordered(site))))
sites <- data.frame(site = sort(unique(prox.in$site)),
                    site.index = seq_along(sort(unique(prox.in$site))))

flattened <- unique(prox.in[c("ai", "site.index")])
flattened <- flattened[order(flattened$ai, flattened$site.index), ]
flattened$ages <- age.indices$age[match(flattened$ai, age.indices$ai)]
flattened$GMST_PhanDA_interp <- approx(prox.in$age, prox.in$GMST_PhanDA,
                                       xout = flattened$ages, rule = 2)$y
flattened$GMST_Scotese21_interp <- approx(prox.in$age, prox.in$GMST_Scotese21,
                                          xout = flattened$ages, rule = 2)$y
flattened$GMST_PhanDA_sd_interp <- approx(prox.in$age, PhanDA_sd,
                                          xout = flattened$ages, rule = 2)$y
cell.offsets <- aggregate(cbind(temp_offset, temp_offset_sd) ~
                            ai + site.index, data = prox.in, FUN = mean, na.rm = TRUE)
cell.key <- paste(cell.offsets$ai, cell.offsets$site.index, sep = "_")
flattened.key <- paste(flattened$ai, flattened$site.index, sep = "_")
offset.rows <- match(flattened.key, cell.key)
stopifnot(all(is.finite(offset.rows)))
flattened$temp_offset_interp <- cell.offsets$temp_offset[offset.rows]
flattened$temp_offset_sd_interp <- cell.offsets$temp_offset_sd[offset.rows]
flattened$row.index <- seq_len(nrow(flattened))
rownames(flattened) <- NULL

clean.d13C <- prox.in[complete.cases(prox.in$d13C), ]
archive.data.raw <- list(
  bf = clean.d13C[clean.d13C$category == "bf", ],
  pf = clean.d13C[clean.d13C$category == "Planktonic foraminifera", ],
  brach = clean.d13C[clean.d13C$category == "Brachiopod calcite", ],
  bivalve = clean.d13C[clean.d13C$category == "Bivalve", ],
  amm = clean.d13C[clean.d13C$category == "Ammonite", ],
  bel = clean.d13C[clean.d13C$category == "Belemnite", ],
  micrite = clean.d13C[clean.d13C$category == "micrite open ocean", ],
  bulk = clean.d13C[clean.d13C$category %in% c("bulk", "bulk open water", "bulk open ocean"), ],
  bulk_sr = clean.d13C[clean.d13C$category == "bulk semi restricted", ],
  bulk_marg = clean.d13C[clean.d13C$category %in%
                           c("bulk marginal sea", "bulk marginal sea restricting up section"), ]
)
if (run.profile == "cenozoic") {
  archive.data.raw <- archive.data.raw[c("bf", "pf", "brach", "bivalve", "bulk")]
}

active.nsb.priors <- nsb.priors[match(names(archive.data.raw), nsb.priors$archive), ]
n.archive.levels <- nrow(active.nsb.priors) + 1L
level.archive.covariance <- matrix(3^2, nrow = n.archive.levels,
                                   ncol = n.archive.levels)
diag(level.archive.covariance) <- c(3^2, 3^2 + active.nsb.priors$sd^2)
level_archive_precision <- solve(level.archive.covariance)

# Exact sufficient-statistic form of the original raw 0.1 per mille likelihood
aggregate.archive <- function(x) {
  cell.key <- interaction(x$ai, x$site.index, drop = TRUE)
  ans <- lapply(split(x, cell.key), function(z) {
    data.frame(ai = z$ai[1], site.index = z$site.index[1], d13C = mean(z$d13C),
               n.obs = nrow(z), d13C.within.sd = if (nrow(z) > 1) sd(z$d13C) else NA_real_,
               d13C.sd = d13C.analyt.sd/sqrt(nrow(z)))
  })
  ans <- do.call(rbind, ans)
  ans <- ans[order(ans$ai, ans$site.index), ]
  rownames(ans) <- NULL
  ans
}

archive.data <- lapply(archive.data.raw, aggregate.archive)

make.archive.index <- function(x) {
  archive.sites <- sort(unique(x$site.index))
  archive.flat <- unique(x[c("ai", "site.index")])
  archive.flat <- archive.flat[order(archive.flat$ai, archive.flat$site.index), ]
  archive.flat$archive.site.index <- match(archive.flat$site.index, archive.sites)
  archive.key <- paste(archive.flat$ai, archive.flat$site.index, sep = "_")
  global.key <- paste(flattened$ai, flattened$site.index, sep = "_")
  observation.key <- paste(x$ai, x$site.index, sep = "_")
  list(n.sites = length(archive.sites),
       site.index = archive.sites,
       n.flat = nrow(archive.flat),
       ai.flat = archive.flat$ai,
       si.flat = archive.flat$archive.site.index,
       ri.flat = match(archive.key, global.key),
       ri.obs = match(observation.key, archive.key),
       ri.global.obs = match(observation.key, global.key))
}

archive.indexes <- lapply(archive.data, make.archive.index)
stopifnot(all(vapply(archive.indexes, function(x) {
  all(is.finite(x$ri.flat)) && all(is.finite(x$ri.obs)) && all(is.finite(x$ri.global.obs))
}, logical(1))))

make.site.map <- function(idx, archive) {
  data.frame(archive = archive,
             archive.site.index = seq_along(idx$site.index),
             site.index = idx$site.index,
             site = sites$site[match(idx$site.index, sites$site.index)])
}
nsb.site.map <- do.call(rbind, Map(make.site.map, archive.indexes, names(archive.indexes)))
rownames(nsb.site.map) <- NULL

make.cell.map <- function(x, idx, archive) {
  data.frame(archive = archive,
             archive.row.index = seq_len(nrow(x)),
             global.row.index = idx$ri.global.obs,
             ai = x$ai,
             age = age.indices$age[match(x$ai, age.indices$ai)],
             site.index = x$site.index,
             site = sites$site[match(x$site.index, sites$site.index)],
             n.obs = x$n.obs,
             d13C.mean = x$d13C,
             d13C.within.sd = x$d13C.within.sd,
             d13C.sd = x$d13C.sd)
}
d13C.cell.map <- do.call(rbind, Map(make.cell.map, archive.data, archive.indexes,
                                    names(archive.data)))
rownames(d13C.cell.map) <- NULL

BWT.Cen <- as.data.frame(read.csv(bwt.file))
names(BWT.Cen) <- c("age", "BWT", "BWT_2sd")
BWT.Cen$age <- BWT.Cen$age*1e3
BWT.Cen <- BWT.Cen[order(BWT.Cen$age, decreasing = TRUE), ]
BWT.Cen.last <- cbind(age.max.spinup, BWT.Cen[1, 2:3])
names(BWT.Cen.last) <- c("age", "BWT", "BWT_2sd")
BWT <- rbind(BWT.Cen, BWT.Cen.last)
BWT <- BWT[order(BWT$age, decreasing = TRUE), ]
BWT.m <- approx(BWT$age, BWT$BWT, xout = ages, method = "linear", rule = 2)$y
BWT.sd <- approx(BWT$age, BWT$BWT_2sd/2, xout = ages, method = "linear", rule = 2)$y

toff.m <- flattened$temp_offset_interp
toff.sd <- flattened$temp_offset_sd_interp
stopifnot(length(ages) == length(ai.all),
          all(is.finite(GMST.m)), all(is.finite(GMST.sd)),
          all(is.finite(BWT.m)), all(is.finite(BWT.sd)),
          all(is.finite(toff.m)), all(is.finite(toff.sd)))

# JAGS inputs
####################################################################################################

data.pass <- list(
  n.steps = n.steps,
  dt.scale = dt.scale,
  si.flat = flattened$site.index,
  ai.flat = flattened$ai,
  GMST.obs = GMST.m,
  GMST.sd = GMST.sd,
  BWT.obs = BWT.m,
  BWT.sd = BWT.sd,
  toff_sd_uniform_bot = toff_sd_uniform_bot,
  toff.m = toff.m,
  toff.sd = toff.sd,
  d13CO2.l = -12,
  d13CO2.u = 0,
  d13CO2_sigma_upper = d13CO2_sigma_upper,
  d13CO2_level_prior_mean = d13CO2_level_prior_mean,
  level_archive_precision = level_archive_precision
)

stems <- c(bf = "d13Cbf", pf = "d13Cpf", brach = "d13Cbrach",
           bivalve = "d13Cbivalve", amm = "d13Camm", bel = "d13Cbel",
           micrite = "d13Cmicrite", bulk = "d13Cbulk",
           bulk_sr = "d13Cbulk_sr", bulk_marg = "d13Cbulk_marg")

for (archive in names(archive.data)) {
  x <- archive.data[[archive]]
  idx <- archive.indexes[[archive]]
  stem <- stems[[archive]]
  prior.row <- nsb.priors[nsb.priors$archive == archive, ]
  data.pass[[paste0(stem, ".data")]] <- x$d13C
  data.pass[[paste0(stem, ".sd")]] <- x$d13C.sd
  data.pass[[paste0("ri.", stem)]] <- idx$ri.obs
  data.pass[[paste0("n.", stem)]] <- nrow(x)
  data.pass[[paste0("n.sites.", archive)]] <- idx$n.sites
  data.pass[[paste0("n.flat.", archive)]] <- idx$n.flat
  data.pass[[paste0("ai.flat.", archive)]] <- idx$ai.flat
  data.pass[[paste0("si.flat.", archive)]] <- idx$si.flat
  data.pass[[paste0("ri.flat.", archive)]] <- idx$ri.flat
  data.pass[[paste0(archive, ".nsb.m")]] <- prior.row$mean
  data.pass[[paste0(archive, ".nsb.sd")]] <- prior.row$sd
}

parms <- c("d13CO2", "d13CO2_level", "GMST", "BWT", "GMST_sigma", "BWT_sigma", "d13CO2_sigma",
           "cc_co2_constant1", "cc_co2_coeff1", "asd", "bpump", "remin",
           "f_co3", "f_carbacid", "Abf", "Asurf")
parms <- c(parms,
           paste0(names(archive.data), ".nsb_mean"),
           paste0(names(archive.data), ".nsb_tau"),
           paste0(names(archive.data), ".nsb_site"),
           paste0(names(archive.data), ".archive_level"))
if (run.settings$coupled) parms <- c(parms, "BWT_GMST_beta")
if (save.full.output) {
  parms <- c(parms, "tempC", "tempC_bot", "toff", "toff_bot",
             unname(stems[names(archive.data)]))
}
if (exists("extra.parms", inherits = FALSE)) parms <- unique(c(parms, extra.parms))

dir.create(output.root, recursive = TRUE, showWarnings = FALSE)
run.time <- system.time({
  inv.out <- do.call(jags.parallel,
                     list(data = data.pass,
                          model.file = model.file,
                          parameters.to.save = parms,
                          inits = NULL,
                          n.chains = n.chains,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          jags.seed = run.seed,
                          jags.module = c("glm", "dic")))
})

run.metadata <- list(
  model.run = model.run,
  output.directory = output.root,
  run.profile = run.profile,
  run.id = run.profile,
  model.version = model.version,
  proxy.file = proxy.file,
  bwt.file = bwt.file,
  paleogeography = run.settings$paleogeography,
  GMST_model = GMST_model,
  temp_offset_model = temp_offset_model,
  use_BWT_GMST_coupling = run.settings$coupled,
  active.archives = names(archive.data),
  likelihood = "archive-site-time means with 0.1/sqrt(n), exactly proportional to the original raw Gaussian likelihood",
  n.proxy.records = sum(vapply(archive.data.raw, nrow, integer(1))),
  age.min = age.min,
  age.max = age.max,
  step.int = step.int,
  n.spinup = n.spinup,
  n.likelihood.cells = nrow(d13C.cell.map),
  n.nsb.sites = nrow(nsb.site.map),
  n.steps = n.steps,
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  run.seed = run.seed,
  process.sd.priors = c(GMST = "uniform(0,0.5)", BWT = "uniform(0,0.5)",
                        d13CO2 = paste0("uniform(0,", d13CO2_sigma_upper, ")")),
  process.time.reference.kyr = 1000,
  time.evolution = "original one-direction first-order random walk",
  parameterization = "exact d13CO2 level/anomaly and likelihood-aligned archive-level block reparameterization",
  temperature.offset.indexing = "site-time Li22 offsets from the selected plate-model file",
  model.file = model.file,
  run.time = run.time,
  completed = Sys.time()
)

if (save.run.output) {
  save(inv.out, ages, ages.short, age.indices, flattened, prox.in, sites,
       nsb.site.map, d13C.cell.map, GMST.m, GMST.sd, BWT.m, BWT.sd, run.metadata,
       file = file.path(output.root, paste0("inv_out_", run.profile, ".rda")))
  saveRDS(run.metadata, file.path(output.root, paste0("run_metadata_", run.profile, ".rds")))

  posterior_qband <- function(parameter) {
    draws <- as.matrix(inv.out$BUGSoutput$sims.list[[parameter]])
    if (nrow(draws) == length(ages)) draws <- t(draws)
    if (ncol(draws) != length(ages)) stop(parameter, " dimensions do not match ages")
    q <- apply(draws, 2, quantile,
               probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    list(q025 = q[1, ], q25 = q[2, ], med = q[3, ],
         q75 = q[4, ], q975 = q[5, ])
  }
  posterior.summary <- list(
    ages = ages,
    GMST.m = GMST.m,
    GMST.sd = GMST.sd,
    BWT.m = BWT.m,
    BWT.sd = BWT.sd,
    run.metadata = run.metadata,
    d13CO2 = posterior_qband("d13CO2"),
    GMST = posterior_qband("GMST"),
    BWT = posterior_qband("BWT")
  )
  saveRDS(posterior.summary,
          file.path(output.root, paste0("posterior_summary_", run.profile, ".rds")))
}

invisible(run.metadata)
