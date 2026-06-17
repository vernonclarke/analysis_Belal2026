# Figure 7. ChI-evoked ACh release is elevated in MCI-Park mice

Fluorescence-biosensor imaging of ChI-evoked acetylcholine release using GRAB-ACh3.0 (electrical
stimulation, 300 µA, 1 ms) in the dorsolateral striatum of MCI-Park (DAT-Cre × Ndufs2⁻/⁻) vs
control mice. Responses are quantified as ΔF/F and are larger in MCI-Park

## Scripts
- `Figure 7 analysis.ipynb` / `Figure 7 reanalysis.ipynb` Python notebooks: NWB loading, ROI
  trace extraction, ΔF/F quantification (see `Fluorescence Biosensor Imaging instructions.md`)
- `Figure 7 analysis and stats.R` robust linear mixed-effects model with cluster bootstrap,
  plus a Bayesian model, accounting for the nested ROI/slice

`results.RData` is a generated R workspace file produced by `Figure 7 analysis and stats.R` after running the imaging statistics; it is not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `Figure 7 dfF.csv` per-ROI peak ΔF/F. Columns: `Animal`, `Slice`, `ROI`, `Condition` (Control/MCI-Park), `dff`, `SliceID`
- `Figure 7 lme fixed effect.csv` fixed-effect estimates from the (robust) linear mixed-effects model
- `Figure 7 robust lme bootstrap.csv` cluster-bootstrap distribution / p-value (N_boot = 9,999)
- `Figure 7 bayes summary.csv` posterior summary of the Bayesian model
- `single_example_animal1_slice1_roi1.xlsx`, `single_example_animal1_slice3_roi1.xlsx` example ΔF/F traces (panel b)
- `Figure 7.xlsx` combined outputs

## Panels `svg/`
- example ΔF/F trace SVGs to **panel b**
- ΔF/F box plot to **panel c**

## Statistics
Robust linear mixed-effects model with cluster bootstrap
control vs MCI-Park:
ACh release fixed effect p = 0.04620 (N_boot = 9,999)
