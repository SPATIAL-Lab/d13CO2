# Generate Figures 2-9

model.output.root <- "output/model_runs/final_archiveblock_3M"
figure.output.root <- "output/figures"

plot.scripts <- file.path("R", "plot", paste0("plot_Fig", 2:9, ".R"))

for (script in plot.scripts) {
  plot.environment <- new.env(parent = globalenv())
  plot.environment$model.output.root <- model.output.root
  plot.environment$figure.output.root <- figure.output.root
  sys.source(script, envir = plot.environment)
}
