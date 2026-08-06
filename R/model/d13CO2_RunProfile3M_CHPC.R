# Run one model profile on CHPC

work <- Sys.getenv("WORK_DIR", unset = getwd())
setwd(work)

profiles <- c("main", "gmst_scotese", "plate_torsvik2017", "plate_merdith2021",
              "plate_cao2024", "cenozoic", "coupled")
run.profile <- Sys.getenv("RUN_PROFILE")
if (!run.profile %in% profiles) stop("Unknown run profile: ", run.profile)

output.root <- "output/model_runs/final_archiveblock_3M"
n.chains <- 4L
n.iter <- 3e6
n.burnin <- 1.5e6
n.thin <- 250L
run.seed <- 26080700L + match(run.profile, profiles)
save.run.output <- TRUE
save.full.output <- FALSE
d13CO2_sigma_upper <- 0.1

source("R/model/d13CO2_Driver_revised.R")
source("R/diagnostics/d13CO2_Diagnostics.R")
