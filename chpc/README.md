# Definitive 3M model suite

Copy the `R` and `chpc` folders into the repository without changing their structure.

Launch the seven profiles from the repository root:

```bash
sbatch chpc/d13CO2_suite_3M.sh
```

The array runs `main`, `gmst_scotese`, `plate_torsvik2017`,
`plate_merdith2021`, `plate_cao2024`, `cenozoic`, and `coupled` as separate jobs.
Model outputs and diagnostic CSVs are written to:

```text
output/model_runs/final_archiveblock_3M/
```

After every array job finishes, combine the diagnostic summaries with:

```bash
Rscript R/diagnostics/d13CO2_CombineDiagnostics.R
```

Generate all figures with:

```bash
Rscript R/plot/RunPlots.R
```
