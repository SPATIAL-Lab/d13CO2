# Combine convergence summaries after the suite finishes

profiles <- c("main", "gmst_scotese", "plate_torsvik2017", "plate_merdith2021",
              "plate_cao2024", "cenozoic", "coupled")
output.root <- "output/model_runs/final_archiveblock_3M"

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
