# Run the definitive model suite

# Choose the folder created under output/model_runs. This can also be set before
# sourcing this script or with the D13CO2_MODEL_RUN environment variable.
if (!exists("model.run", inherits = FALSE)) model.run <- NULL
source("R/model/d13CO2_RunPaths.R", local = TRUE)
model.run <- d13CO2_model_run_name(model.run)

all.profiles <- c("main", "gmst_scotese", "plate_torsvik2017", "plate_merdith2021",
                  "plate_cao2024", "cenozoic", "coupled")
requested.profile <- Sys.getenv("RUN_PROFILE", unset = "")
profiles <- if (nzchar(requested.profile)) requested.profile else all.profiles
if (!all(profiles %in% all.profiles)) {
  stop("Unknown RUN_PROFILE: ", paste(setdiff(profiles, all.profiles), collapse = ", "))
}
output.root <- d13CO2_model_run_dir(model.run, create = TRUE)

for (j in seq_along(profiles)) {
  run.profile <- profiles[j]
  n.chains <- 4L
  n.iter <- 3e6
  n.burnin <- 1.5e6
  n.thin <- 250L
  run.seed <- 26080700L + match(run.profile, all.profiles)
  save.run.output <- TRUE
  save.full.output <- FALSE
  d13CO2_sigma_upper <- 0.1

  source("R/model/d13CO2_Driver.R")
  source("R/diagnostics/d13CO2_Diagnostics.R")
  rm(inv.out)
  invisible(gc())
}

if (identical(profiles, all.profiles)) {
  source("R/diagnostics/d13CO2_CombineDiagnostics.R")
}
