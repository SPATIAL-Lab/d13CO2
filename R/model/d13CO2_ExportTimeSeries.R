# Export atmospheric carbon-isotope posterior time series

if (!exists("model.output.root", inherits = FALSE)) {
  if (!exists("model.run", inherits = FALSE)) model.run <- NULL
  source("R/model/d13CO2_RunPaths.R", local = TRUE)
  model.output.root <- d13CO2_model_run_dir(model.run, must.exist = TRUE)
}
if (!exists("series.output.root", inherits = FALSE)) series.output.root <- "output"

export_d13CO2 <- function(profile, output.file, prefix) {
  summary.file <- file.path(model.output.root,
                            paste0("posterior_summary_", profile, ".rds"))
  summary <- readRDS(summary.file)
  required <- c("q025", "q25", "med", "q75", "q975")
  if (!all(required %in% names(summary$d13CO2))) {
    stop("Posterior summary does not contain the 50% and 95% intervals: ",
         summary.file)
  }

  ages <- summary$ages
  keep <- ages >= 0 & ages <= summary$run.metadata$age.max
  q <- summary$d13CO2
  result <- data.frame(age_Ma = ages[keep]/1000,
                       q025 = q$q025[keep], q25 = q$q25[keep],
                       q50 = q$med[keep], q75 = q$q75[keep],
                       q975 = q$q975[keep])
  names(result)[-1] <- paste0(prefix, c("_2p5", "_25", "_50", "_75", "_97p5"))
  result <- result[order(result$age_Ma), ]

  dir.create(series.output.root, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(series.output.root, output.file), row.names = FALSE)
  invisible(result)
}

export.specifications <- list(
  main = list(file = "d13Ca_Phan_1Myr.csv", prefix = "d13Ca"),
  cenozoic = list(file = "d13Ca_Cen_100kyr.csv", prefix = "d13CO2")
)
if (!exists("export.profiles", inherits = FALSE)) {
  export.profiles <- names(export.specifications)
}
for (profile in export.profiles) {
  specification <- export.specifications[[profile]]
  export_d13CO2(profile, specification$file, specification$prefix)
}
