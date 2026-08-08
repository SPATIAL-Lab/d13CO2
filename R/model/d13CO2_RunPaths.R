# Shared model-run path handling

d13CO2_model_run_name <- function(model.run = NULL) {
  if (is.null(model.run) || !length(model.run) || !nzchar(model.run)) {
    model.run <- Sys.getenv("D13CO2_MODEL_RUN", unset = "definitive")
  }

  if (length(model.run) != 1L || is.na(model.run) || !nzchar(model.run) ||
      model.run %in% c(".", "..") || basename(model.run) != model.run ||
      grepl("[/\\\\]", model.run)) {
    stop("model.run must be one folder name under output/model_runs: ",
         paste(model.run, collapse = ", "))
  }

  model.run
}

d13CO2_model_run_dir <- function(model.run = NULL,
                                  model.runs.root = file.path("output", "model_runs"),
                                  create = FALSE,
                                  must.exist = FALSE) {
  model.run <- d13CO2_model_run_name(model.run)
  run.dir <- file.path(model.runs.root, model.run)

  if (create) {
    dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (must.exist && !dir.exists(run.dir)) {
    stop("Selected model run does not exist: ", run.dir,
         "\nRun R/model/d13CO2_RunAll.R first or choose another model.run.")
  }

  run.dir
}
