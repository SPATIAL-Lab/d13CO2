# Definitive model suite

Submit from the cloned repository root. The launcher uses the submission directory by default; `WORK_DIR` can override it. Install the required R packages in your own library (`R_LIBS_USER` may be supplied through the environment). The Slurm account, partition, email, and module settings are CHPC-specific; adjust them for your cluster.

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
