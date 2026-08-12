# d13CO2

Bayesian proxy-system model for reconstructing the carbon-isotope composition of Phanerozoic atmospheric CO2 from marine carbonate records. The model was developed by Dustin Harper and Gabriel Bowen with input from the Phanerozoic CO2 Proxy Integration Project community.

## Citation

Harper, D. T., and Bowen, G. J. (2026). d13CO2: Data, model output, and code for “Reconstructing the carbon isotope composition of the Phanerozoic atmosphere” (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.21910191

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21910191.svg)](https://doi.org/10.5281/zenodo.21910191)

## Software

The definitive analysis used JAGS 4.3.1 with four parallel chains. The R workflow requires `R2jags`, `rjags`, and `coda`; figure generation additionally uses `biwavelet`, `readxl`, and `viridisLite`. The preprocessing scripts use `dplyr`, `openxlsx`, and `ncdf4`. Run all commands from the repository root.

The Li et al. (2022) climate-simulation NetCDF is stored with Git LFS. Install Git LFS and run `git lfs pull` after cloning the repository.

## Definitive analysis

The suite contains seven profiles:

| Profile | Purpose | Paleogeography | GMST | Resolution |
| --- | --- | --- | --- | --- |
| `main` | Primary Phanerozoic reconstruction | PALEOMAP | PhanDA | 1 Myr |
| `gmst_scotese` | GMST sensitivity | PALEOMAP | Scotese et al. (2021) | 1 Myr |
| `plate_torsvik2017` | Paleogeographic sensitivity | Torsvik and Cocks (2017) | PhanDA | 1 Myr |
| `plate_merdith2021` | Paleogeographic sensitivity | Merdith et al. (2021) | PhanDA | 1 Myr |
| `plate_cao2024` | Paleogeographic sensitivity | Cao et al. (2024) | PhanDA | 1 Myr |
| `cenozoic` | Higher-resolution Cenozoic reconstruction | PALEOMAP | PhanDA | 100 kyr |
| `coupled` | Coupled GMST-BWT sensitivity | PALEOMAP | PhanDA | 1 Myr |

All profiles use the Li et al. (2022) simulations to define local temperature offsets. The definitive configuration runs four chains for 3,000,000 iterations, discards the first 1,500,000, and retains every 250th sample.

Run the seven model profiles sequentially:

```r
model.run <- "definitive"
source("R/model/d13CO2_RunAll.R")
```

On CHPC, `chpc/d13CO2_suite.sh` submits the seven profiles as independent array jobs.

Generate Figures 2-9 from a completed run:

```r
model.run <- "definitive"
source("R/plot/RunPlots.R")
```

Run the models and figures in one call:

```r
source("RunFullAnalysis.R")
```

Model outputs are written to `output/model_runs/<model.run>/`; figures and exported time series are written to `output/`. The public results from the definitive CHPC run are in `results/definitive/`, and final figure PDFs are in `results/figures/`.

## Data workflow

The model compilation contains 63,112 marine carbonate measurements assigned to ten archive categories. The data workflow was:

1. **Compile and harmonize records.** Bulk-carbonate data were assembled primarily from Cramer and Jarvis (2020). Component-specific data were assembled from the StabisoDB compilation described by Grossman et al. (2025), with the additional records described in the manuscript. Ages were placed on the GTS2020 timescale (Gradstein et al., 2020).
2. **Standardize sites and archives.** Study locations were assigned consistent site names, modern latitude and longitude, and archive categories. `R/preprocess/SiteAssigner.R` joins the curated study-level location tables and produces the bulk compilation.
3. **Integrate component data.** The harmonized bulk and component records were combined into the final 63,112-row compilation used for paleogeographic reconstruction.
4. **Assign paleocoordinates.** GPlates was used to reconstruct every site with the PALEOMAP, Torsvik and Cocks (2017), Merdith et al. (2021), and Cao et al. (2024) plate models. The exported coordinates for all four models are distributed together.
5. **Calculate local temperature offsets.** Local mean annual temperature and global mean surface temperature were extracted from the 10 Myr CESM simulations of Li et al. (2022) and linearly interpolated to each record's age and paleolocation. Their difference defines the local temperature offset. The primary GMST prior and uncertainty are from the PhanDA reconstruction of Judd et al. (2024); the Scotese et al. (2021) curve is used only in the GMST sensitivity profile. `R/preprocess/Phan_MAT_extractor.R` contains the extraction and interpolation calculation for the Cao et al. (2024) example, and the completed model-ready tables for all four plate models are provided.
6. **Run the inversion.** The driver selects the appropriate processed table, temporal resolution, GMST prior, and proxy-system model for each profile. Cenozoic bottom-water-temperature constraints are supplied separately from Meckler et al. (2022).

The age harmonization, source-data integration, and GPlates reconstruction involved curated steps rather than a single fully automated script. The final compilation, paleocoordinates, processed model inputs, and source-reference table are therefore provided directly so that the inversion begins from a fixed and documented dataset.

## Data files

### Compilation and paleogeography

| File | Role in the workflow |
| --- | --- |
| `data/compilation/paleozoicAssignSites.xlsx` | Paleozoic bulk records and curated study-to-site assignments. |
| `data/compilation/nonpaleozoicAssignSites.xlsx` | Non-Paleozoic bulk records and curated study-to-site assignments. |
| `data/compilation/updated_paleozoicAssignSites.xlsx` | Paleozoic records after modern site coordinates and archive categories were joined. |
| `data/compilation/updated_nonpaleozoicAssignSites.xlsx` | Non-Paleozoic records after modern site coordinates and archive categories were joined. |
| `data/compilation/PhanCompUpdated_bulk.csv` | Combined bulk-carbonate compilation produced by `SiteAssigner.R`. |
| `data/compilation/Phan_component_d13C.xlsx` | Component-specific carbonate compilation. |
| `data/compilation/PhanCompUpdated_WithComponent.csv` | Final 63,112-row bulk-plus-component compilation. |
| `data/compilation/references.csv` | Source-reference lookup for the data compilation. |
| `data/paleogeography/phancomp_paleocoord.csv` | Final compilation with paleolongitude and paleolatitude from all four plate models. |

### Climate and model inputs

| File | Role in the workflow |
| --- | --- |
| `data/raw/High_Resolution_Climate_Simulation_Dataset_540_Myr.nc` | Li et al. (2022) CESM climate fields at 10 Myr intervals. |
| `data/raw/PhanDA_GMST.csv` | PhanDA GMST median and uncertainty used for the primary prior. |
| `data/raw/CenozoicBWT.csv` | Cenozoic bottom-water-temperature constraints. |
| `data/processed/PhanCompWithTemp_PALEOMAP.csv` | Main PALEOMAP model input with Li et al. temperature offsets. |
| `data/processed/PhanCompWithTemp_TorsvikCocks2017.csv` | Torsvik and Cocks plate-model sensitivity input. |
| `data/processed/PhanCompWithTemp_MERDITH2021.csv` | Merdith plate-model sensitivity input. |
| `data/processed/PhanCompWithTemp_CAO2024.csv` | Cao plate-model sensitivity input. |

### Figure source data

`data/raw/ISOORG16.xlsx`, `ISOORG23.xlsx`, `Fosteretal17.xlsx`, `lateCenozoic_plantwax_d13C.xlsx`, `CenCO2PIP_500kyrCO2.csv`, and `CenCO2PIP_500kyrTemp.csv` provide the comparison records used in Figures 5 and 6. `data/processed/Tipple2010_Fig4_benthic_d13CO2.csv` digitizes the benthic-based atmospheric reconstruction shown in Figure 4 of Tipple et al. (2010).

## Definitive results

The repository distributes compact posterior summaries for all seven profiles, run metadata, convergence diagnostics, final figures, and two publication-ready atmospheric carbon-isotope time series:

- `results/definitive/time_series/d13Ca_Phan_1Myr.csv` contains the primary Phanerozoic reconstruction.
- `results/definitive/time_series/d13CO2_Cen_100kyr.csv` contains the higher-resolution Cenozoic reconstruction.

Each CSV reports age in Ma and the 2.5th, 25th, 50th, 75th, and 97.5th posterior percentiles. The 25th-75th columns define the 50% credible interval, the 2.5th-97.5th columns define the 95% credible interval, and the 50th percentile is the posterior median.

The complete `inv_out_main.rda` and `inv_out_cenozoic.rda` objects are stored in `results/definitive/full_outputs/` using Git LFS. These objects support the primary reconstruction, higher-resolution Cenozoic reconstruction, Figures 2-6, and the two CSV exports. Compact posterior summaries in `results/definitive/posterior_summaries/` provide the reported time-series intervals for all seven profiles without requiring the remaining full MCMC objects. Convergence files are in `results/definitive/diagnostics/`, and run settings are in `results/definitive/metadata/`.

To regenerate the figures from the distributed results without rerunning JAGS, place the two full outputs and all seven posterior summaries together in `output/model_runs/definitive/`, then source `R/plot/RunPlots.R` with `model.run <- "definitive"`.

## Primary data sources

- [Cramer and Jarvis (2020)](https://doi.org/10.1016/B978-0-12-824360-2.00011-5)
- [Gradstein et al. (2020), *Geologic Time Scale 2020*](https://doi.org/10.1016/C2018-0-02379-3)
- [Grossman et al. (2025)](https://doi.org/10.1073/pnas.2424291122) and [StabisoDB](https://stabisodb.org/)
- [Judd et al. (2024)](https://doi.org/10.1126/science.adk3705)
- [Li et al. (2022)](https://doi.org/10.1038/s41597-022-01490-4)
- [Merdith et al. (2021)](https://doi.org/10.1016/j.earscirev.2020.103477)
- [Torsvik and Cocks (2017)](https://doi.org/10.1017/9781316225523)
- [Cao et al. (2024)](https://doi.org/10.1016/j.gsf.2024.101922)

## License

See `LICENSE`.
