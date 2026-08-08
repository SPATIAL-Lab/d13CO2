# Run the definitive model suite and generate Figures 2-9

# Choose the folder created under output/model_runs and used by the figures.
if (!exists("model.run", inherits = FALSE)) {
  model.run <- Sys.getenv("D13CO2_MODEL_RUN", unset = "definitive")
}

source("R/model/d13CO2_RunAll.R")
source("R/plot/RunPlots.R")
