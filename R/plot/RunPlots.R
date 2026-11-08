# Generate Figures 2-8 and S1

if (!exists("model.run", inherits = FALSE)) model.run <- NULL
source("R/model/d13CO2_RunPaths.R", local = TRUE)
model.run <- d13CO2_model_run_name(model.run)
model.output.root <- d13CO2_model_run_dir(model.run)
full.output.root <- model.output.root
if (model.run == "definitive" && !dir.exists(model.output.root)) {
  model.output.root <- "results/definitive/posterior_summaries"
  full.output.root <- "results/definitive/full_outputs"
}
if (!exists("figure.output.root", inherits = FALSE)) {
  figure.output.root <- file.path("output", "figures")
}

required.profiles <- c("main", "gmst_scotese", "plate_torsvik2017",
                       "plate_merdith2021", "plate_cao2024", "cenozoic", "coupled")
required.files <- c(
  file.path(full.output.root, c("inv_out_main.rda", "inv_out_cenozoic.rda")),
  file.path(model.output.root, paste0("posterior_summary_", required.profiles, ".rds"))
)
missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files)) {
  stop("Cannot generate all figures from model run '", model.run,
       "'; missing files:\n", paste(missing.files, collapse = "\n"))
}

series.environment <- new.env(parent = globalenv())
series.environment$model.output.root <- model.output.root
series.environment$series.output.root <- "output"
sys.source("R/model/d13CO2_ExportTimeSeries.R", envir = series.environment)

plot.scripts <- file.path("R", "plot", paste0("plot_Fig", 2:9, ".R"))

for (script in plot.scripts) {
  plot.environment <- new.env(parent = globalenv())
  plot.environment$model.output.root <- if (basename(script) %in% paste0("plot_Fig", 2:6, ".R")) full.output.root else model.output.root
  plot.environment$figure.output.root <- figure.output.root
  sys.source(script, envir = plot.environment)
}
