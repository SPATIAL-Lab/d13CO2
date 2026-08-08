# Combine convergence summaries after the suite finishes

if (!exists("model.run", inherits = FALSE)) model.run <- NULL
source("R/model/d13CO2_RunPaths.R", local = TRUE)
model.run <- d13CO2_model_run_name(model.run)
output.root <- d13CO2_model_run_dir(model.run, must.exist = TRUE)

profiles <- c("main", "gmst_scotese", "plate_torsvik2017", "plate_merdith2021",
              "plate_cao2024", "cenozoic", "coupled")

required.files <- unlist(lapply(
  c("diagnostics_summary_", "scalar_diagnostics_"),
  function(prefix) file.path(output.root, paste0(prefix, profiles, ".csv"))
))
missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files)) {
  stop("Cannot combine diagnostics; missing files:\n", paste(missing.files, collapse = "\n"))
}

convergence.summary <- do.call(rbind, lapply(profiles, function(profile) {
  x <- read.csv(file.path(output.root,
                          paste0("diagnostics_summary_", profile, ".csv")))
  cbind(run.profile = profile, x)
}))

scalar.summary <- do.call(rbind, lapply(profiles, function(profile) {
  x <- read.csv(file.path(output.root,
                          paste0("scalar_diagnostics_", profile, ".csv")))
  cbind(run.profile = profile, x)
}))

write.csv(convergence.summary,
          file.path(output.root, "convergence_summary.csv"), row.names = FALSE)
write.csv(scalar.summary,
          file.path(output.root, "scalar_convergence_summary.csv"), row.names = FALSE)
