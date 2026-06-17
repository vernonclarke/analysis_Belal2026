# Figure 10. NPY-evoked GABAergic input is intact after 6-OHDA lesion

Control experiment: optogenetic activation (Chronos) of NPY interneurons in NPY-Cre × D1-tdTomato
mice evokes GABAergic PSCs in SPNs that are **unchanged** after 6-OHDA lesion (no change in
amplitude, charge transfer or kinetics), in contrast to the ChI pathway (Fig 9). Both dSPNs and
iSPNs; comparison unpaired (control vs 6-OHDA)

## Scripts (run order)
1. `Figure 10 data processing from downloaded dandiset.R` extracts/averages PSCs to per-condition `xlsx`
2. Per-condition fits:
   - `Figure 10 NPY Cre X dSPN WT.R` / `... 6OHDA.R`
   - `Figure 10 NPY Cre X iSPN WT.R` / `... 6OHDA.R`
3. `Figure 10 graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`NPY Cre X dSPN WT.RData`, `NPY Cre X dSPN 6OHDA.RData`, `NPY Cre X iSPN WT.RData` and `NPY Cre X iSPN 6OHDA.RData` are generated R workspace files produced by the per-condition Figure 10 fitting scripts after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `NPY Cre X dSPN WT[ /6OHDA].xlsx` / `... full.xlsx`, `NPY Cre X iSPN WT[ /6OHDA].xlsx` / `... full.xlsx` per-cell averaged traces
- `Figure 10 NPY_Cre_X_*.csv`  NLS 2 PSC component fit summaries (one per group)
- `Figure 10 amplitude.csv`, `Figure 10 charge transfer.csv` peak amplitude / charge transfer
- `Figure 10 amplitude 2 components.csv`, `Figure 10 tau 2 components.csv` fast/slow amplitudes and time constants
- `Figure 10 iSPN single examples control.csv`, `... 6OHDA.csv` example traces (panels c/d)
- `Figure 10 statistics.csv`, `Figure 10 additional statistics.csv` Wilcoxon rank-sum results
- `Figure 10.xlsx` combined workbook

## Panels `svg/`
- `egs_iSPN.svg`, `egs_iSPN_6OHDA.svg` (+ semilog) to **panels c, d**
- amplitude / charge-transfer box plots (+ semilog) to **panels e, f**
- `scatter_dSPN.svg`, `scatter_iSPN.svg`, `scatter_logxy_*` to **panel g**
- `tau_*` box plots to **panel h**
- `NPY Cre X *_<n>.svg` / `NPY_Cre_X_*_<n>.svg` per-cell fit diagnostics (note: both space- and underscore-named copies are present)
- `dbscan_*.svg` clustering QC

## Statistics
Wilcoxon rank-sum test
control vs 6-OHDA:
dSPN amplitude W = 55, p = 1.00000 (n = 13, 5 animals control; n = 8, 3 animals 6-OHDA)
iSPN amplitude W = 40, p = 0.57862 (n =  6, 2 animals control; n = 9, 3. animals 6-OHDA)
dSPN charge transfer W = 30, p = 0.36376 (n = 13, 5 animals control; n = 8, 3 animals 6-OHDA)
iSPN charge transfer W = 8, p = 0.10230  (n = 6,  2 animals control; n = 9, 3. animals 6-OHDA)
dSPN fast tau W = 34, p = 0.84104
iSPN fast tau W = 26, p = 1.00000
dSPN slow tau W = 33, p = 0.73867
iSPN slow tau W = 23, p = 1.00000
