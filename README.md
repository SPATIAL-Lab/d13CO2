# d13CO2

Proxy-system model for reconstructing the carbon-isotope composition of Phanerozoic atmospheric CO2 from marine carbonate records. The model was developed by Dustin Harper and Gabriel Bowen with input from the Phanerozoic CO2 Proxy Integration Project community.

## Analysis

Run commands from the repository root. The model requires JAGS and the R package `R2jags`. Figure generation also uses `biwavelet` and `readxl`.

The definitive analysis contains seven profiles: the main reconstruction, an alternative GMST prior, three paleogeographic-model sensitivity tests, a high-resolution Cenozoic reconstruction, and a coupled GMST-BWT sensitivity test.

Run the models only:

```r
source("R/model/d13CO2_RunAll_revised.R")
```

Generate Figures 2-9 from saved model outputs:

```r
source("R/plot/RunPlots.R")
```

Run the complete analysis and then generate all figures:

```r
source("RunFullAnalysis.R")
```

Model outputs are written to `output/model_runs/final_archiveblock_3M/`, and figures are written to `output/figures/`.

## Data processing

The current preprocessing scripts are in `R/preprocess/`. Earlier paleogeographic reconstruction files and processing provenance are retained locally under `legacy/`.

The 719 MB climate-simulation NetCDF is configured for Git LFS in `.gitattributes`; initialize Git LFS before adding the data to a remote repository.
