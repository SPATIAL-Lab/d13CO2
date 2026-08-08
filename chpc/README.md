# Definitive model suite

Copy the `R` and `chpc` folders into the repository without changing their structure.

Choose a folder name under `output/model_runs/` and launch the seven profiles from the repository root:

```bash
sbatch --export=ALL,D13CO2_MODEL_RUN=definitive chpc/d13CO2_suite.sh
```

The array runs `main`, `gmst_scotese`, `plate_torsvik2017`,
`plate_merdith2021`, `plate_cao2024`, `cenozoic`, and `coupled` as separate jobs.
Each array task uses `d13CO2_RunAll.R` with `RUN_PROFILE` set to one profile. Model outputs and diagnostic CSVs are written to:

```text
output/model_runs/<D13CO2_MODEL_RUN>/
```

After every array job finishes, combine the diagnostic summaries with:

```bash
D13CO2_MODEL_RUN=definitive Rscript R/diagnostics/d13CO2_CombineDiagnostics.R
```

Generate all figures with:

```bash
D13CO2_MODEL_RUN=definitive Rscript R/plot/RunPlots.R
```
