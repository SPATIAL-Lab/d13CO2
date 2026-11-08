# Schartman et al. (2026): Figure 6 data

Anna K. Schartman, Pratigya J. Polissar, Caroline A. E. Strömberg. The origin and development of the Miocene northwest African savanna. PNAS 123(33), e2527340123. https://doi.org/10.1073/pnas.2527340123

## Provenance and paths

The publisher-supplied article archived at https://pmc.ncbi.nlm.nih.gov/articles/PMC13486532/ explicitly identifies Dataset S3 as containing raw and corrected C25-C35 isotope values and links the four filenames below. Automated supplement retrieval did not produce the spreadsheets. The user supplied these files from Downloads on 2026-11-06; their contents and filenames agree with the article and SI. They were copied unchanged. Independent byte comparison with the publisher download was not possible.

- data/raw/Schartman2026/pnas.2527340123.sd01.xlsx: age models, splice, sample IDs and depths.
- data/raw/Schartman2026/pnas.2527340123.sd02.xlsx: replicate n-alkane concentrations; archived, not used as isotope data.
- data/raw/Schartman2026/pnas.2527340123.sd03.xlsx: isotope measurements and correction fields.
- data/raw/Schartman2026/pnas.2527340123.sapp.pdf: methods, screening and age-model limitations.
- data/processed/Schartman2026_plantwax_raw.csv: all 164 sample rows, including the six excluded from plotting; metadata and raw odd-chain C25-C35 measurements with SEM.

SHA256 (original and archived copies):
- sd01.xlsx: 8f167efd68fd8d3c637321dc1884a494ef5357a6194c75944e0d7f14c1e50f65
- sd02.xlsx: 7a1a054ed3eb86ab8947170dfbf46b6534e60e2e96fc412d050df6b7417bde5f
- sd03.xlsx: 21ea4872efacd9c72e04d56811b33d8214be36d35d4d671dbfdd2406a694adfc
- sapp.pdf: ce0bfdc1d4c345a781cc89bc1532aebf9be9eb62974e26d59b5c5608f17ff141

## Extraction and correction check

S3 Sheet1 rows 2-165 are samples; rows 166-167 are footnotes. Sample IDs are unique. Columns B/C/E supply sample ID/core ID/age (Ma). Raw C31 is AB, SEM is AC; AM is atmosphere-corrected C31 and is NOT used in Figure 6. Raw/SEM pairs for C25, C27, C29, C33, C35 are V/W, X/Y, Z/AA, AD/AE, AF/AG. Source row and DOI are retained in each extracted row. S1 Samples was joined by exact sample ID to retain mbsf/mcd depths; all 164 ages agree with S3 within 1e-10 Ma.

S3 AH contains the authors' Tipple atmospheric estimate and AI the -6.5 per mil preindustrial baseline. Every C31 corrected value satisfies

    corrected = (raw + 1000) * (preindustrial + 1000) / (atmosphere + 1000) - 1000

within 3.2e-13 per mil. Corrected-minus-raw ranges from -1.215593 to -0.181001 per mil. Using AM as raw would therefore introduce a second atmospheric correction. Figure 6 uses AB only and calculates model and Tipple differences independently from that same raw measurement. The existing digitized Tipple Fig. 4 curve is retained for consistency across all records; S3 AH is preserved for audit, not substituted for that common baseline.

Raw C31 exists for all 164 samples (range -31.6633 to -25.1596 per mil VPDB); all have SEM (0.0508836-0.0863534 per mil). These are sample-level reported estimates, not individual instrument injections. SI S4 describes duplicate/triplicate measurements and pooled uncertainty including VPDB scale assignment. No per-sample age uncertainty is supplied in S3 or S1 Samples; the extraction leaves it blank. S1 tie-point errors are not sample errors and were not propagated. S3 labels all depths mcd, but DSDP 368 depths match S1 mbsf, where mcd is blank; both original labels and S1 depths are retained.

## Coverage and screening

- ODP 659: 122 raw samples, 5.476635-28.041905 Ma; 120 plotted.
- DSDP 368: 42 raw samples, 3.954716-23.890837 Ma; 38 plotted.
- CPI >= 3 follows SI S6; six excluded IDs: UC4438, UC4447, UC4957, UC4958, UC4941, UC4966. All remain in the CSV with include_CPI=FALSE.
- No duplicate sample IDs, missing C31 values or exact age/raw pairs matching the existing compilation were found. Existing compilation has only age/d13C/category/reference, so a sample-ID-level cross-study duplicate audit is impossible. The added rows are this study's ODP/DSDP projects, not a re-import of the older Polissar compilation.
- Ages older than 24.4 Ma at ODP 659 are extrapolated beyond the final tie point (SI S2); flagged in the CSV. DSDP 368 age control is weaker and uses broad biozone/depth ranges.

## Figure 6 changes and reproduction

R/plot/plot_Fig6.R retains its base R structure. Added CSV import and CPI selection, one marine Schartman_C31 group (purple open diamonds), and its legend entry. Changed both Delta calculations to atmosphere minus raw wax, reference lines to +24.7/+21.8/+15.1, and the panel-c label to Delta_atm-leaf. Small legend spacing, panel-c letter and bottom axis-label adjustments prevent overlap. Existing input workbook and Tipple curve are unchanged. No model rerun or new statistical analysis.

From the repository root, set `model.output.root <- "results/definitive/full_outputs"` and source `R/plot/plot_Fig6.R`. The script writes `output/figures/Figure6.pdf` and requires R, readxl, and the distributed model/data files (including the full outputs downloaded with Git LFS). Other `model.output.root` and `figure.output.root` choices remain supported.

## Inspection and interpretation

1507 raw-wax rows plotted, including 158 new samples. The existing model interpolation leaves 24 old records without Delta; all new records have model and Tipple values. Tipple covers all plotted ages. Numerical sign checks passed. The PDF uses vector marks and embedded fonts, with no raster images; Illustrator import itself was not tested.

For inspection only, 14-17 Ma contains 27 screened Schartman samples. Their model Delta spans 24.489-25.629 per mil, compared with 24.500-25.324 using Tipple. They cluster near the C3 line and do not reproduce the low-Delta DSDP 94 excursion. This regional comparison does not establish a global MCO physiological response: source vegetation, C3/C4 mixing, wax biosynthetic offsets, geography, chain length and chronology differ. Existing DSDP94 C29 and C31 remain distinct. Delta here is the requested simple atmospheric-minus-wax difference, not ratio-normalized discrimination or a direct bulk-leaf measurement. Reference magnitudes were sign-reversed, not recalibrated. SEM and age uncertainties are preserved where available but not propagated into panel c; the model credible bands remain in panel b. S3's corrected values were retained only for the correction audit.
