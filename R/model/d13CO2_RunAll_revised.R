# Run the definitive model suite

profiles <- c("main", "gmst_scotese", "plate_torsvik2017", "plate_merdith2021",
              "plate_cao2024", "cenozoic", "coupled")

for (j in seq_along(profiles)) {
  run.profile <- profiles[j]
  output.root <- "output/model_runs/final_archiveblock_3M"
  n.chains <- 4L
  n.iter <- 3e6
  n.burnin <- 1.5e6
  n.thin <- 250L
  run.seed <- 26080700L + j
  save.run.output <- TRUE
  save.full.output <- FALSE
  d13CO2_sigma_upper <- 0.1

  source("R/model/d13CO2_Driver_revised.R")
  rm(inv.out)
  invisible(gc())
}
