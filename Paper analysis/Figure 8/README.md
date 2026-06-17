# Figure 8. ChI-evoked GABAergic PSCs are attenuated in MCI-Park mice

In the prodromal MCI-Park model (ChAT-Flp × Ndufs2 fl/fl × DAT-Cre), optogenetic ChI activation
evokes smaller GABAergic PSCs in SPNs than in controls (reduced amplitude and charge transfer),
with unchanged decay kinetics. Comparison is unpaired (control vs MCI-Park)

## Scripts (run order)
1. `Figure 8 data processing from downloaded dandiset.R` extracts/averages PSCs to per-condition `xlsx`
2. `Figure 8 ctrl.R` (control) and `Figure 8 test.R` (MCI-Park) perform two-component fits
3. `Figure 8 graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`Control for MCI-Park.RData` and `ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK.RData` are generated R workspace files produced by `Figure 8 ctrl.R` and `Figure 8 test.R` after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `Control for MCI-Park.xlsx` / `... full.xlsx`, `ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK.xlsx` / `... full.xlsx` per-cell averaged traces
- `Figure 8 amplitude.csv`, `Figure 8 charge transfer.csv` peak amplitude / charge transfer
- `Figure 8 amplitude 2 components.csv`, `Figure 8 tau 2 components.csv` fast/slow amplitudes and time constants
- `Figure 8 ctrl single examples.csv`, `Figure 8 test single examples.csv` example traces (panels b/c)
- `Figure 8 statistics.csv` Wilcoxon rank-sum results
- `Figure 8.xlsx` combined workbook

## Panels `svg/`
- `egs_*` (control / MCI-Park, + semilog) to **panels b, c**
- amplitude / charge-transfer box plots (+ semilog) to **panel d**
- scatter (slow vs fast amplitude) to **panel e**
- tau box plots to **panel f**
- per-cell fit-diagnostic SVGs; `dbscan_*` clustering QC

## Statistics
Wilcoxon rank-sum test
control vs MCI-Park:

amplitude W = 0, p = 0.00067
charge transfer W = 48, p = 0.00067
fast tau W = 15, p = 0.56477
slow tau W = 29, p = 0.57276

One MCI-Park cell had a slow-only fit (shown red, unpaired); displayed n adjusted accordingly
