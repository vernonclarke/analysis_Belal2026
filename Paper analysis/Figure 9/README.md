# Figure 9. ChI-evoked GABAergic input is reduced after 6-OHDA (MFB) lesion

Unilateral 6-OHDA lesion of the medial forebrain bundle in ChAT-Cre × D1-tdTomato mice reduces
the optogenetically ChI-evoked GABAergic PSC in both dSPNs and iSPNs (amplitude and charge
transfer), with unchanged decay kinetics. Comparison is unpaired (control vs 6-OHDA), separately
for each cell type

## Scripts (run order)
1. `Figure 9 data processing from downloaded dandiset.R` extracts/averages PSCs to per-condition `xlsx`
2. Per-condition fits:
   - `Figure 9 ChAT-Cre X tdTomato.R` / `... 6OHDA.R` dSPN control / lesion
   - `Figure 9 ChAT-Cre X De eGFP.R` / `... 6OHDA.R` iSPN control / lesion
3. `Figure 9 graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`ChAT-Cre X tdTomato.RData`, `ChAT-Cre X tdTomato 6OHDA.RData`, `ChAT-Cre X De eGFP.RData` and `ChAT-Cre X De eGFP 6OHDA.RData` are generated R workspace files produced by the per-condition Figure 9 fitting scripts after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `ChAT-Cre X tdTomato[ 6OHDA].xlsx` / `... full.xlsx`, `ChAT-Cre X De eGFP[ 6OHDA].xlsx` / `... full.xlsx` per-cell averaged traces (dSPN/iSPN × control/6-OHDA)
- `Figure 9 ChAT-Cre X_*.csv`  NLS 2 PSC component fit summaries (one per group)
- `Figure 9 amplitude.csv`, `Figure 9 charge transfer.csv` peak amplitude / charge transfer
- `Figure 9 amplitude 2 components.csv`, `Figure 9 tau 2 components.csv` fast/slow amplitudes and time constants
- `Figure 9 dSPN single examples control.csv`, `... 6OHDA.csv` example traces (panels c/d)
- `Figure 9 statistics.csv`, `Figure 9 additional statistics.csv` Wilcoxon rank-sum results
- `Figure 9.xlsx` combined workbook

## Panels `svg/`
- `egs_*` (control / 6-OHDA, + semilog) to **panels c, d**
- amplitude box plots to **panel e**; charge-transfer box plots to **panel f**
- scatter (slow vs fast amplitude; open = dSPN, closed = iSPN) to **panel g**
- tau box plots (dSPN/iSPN, fast/slow) to **panel h**
- per-cell fit-diagnostic SVGs; `dbscan_*` clustering QC

## Statistics
Wilcoxon rank-sum test
control vs 6-OHDA:

dSPN amplitude W = 6, p = 0.00260 (n = 10, 5 animals control; n = 9, 3 animals 6-OHDA)
iSPN amplitude W = 0, p = 0.01299 (n =  6, 4 animals control; n = 5, 4 animals 6-OHDA)
dSPN scatter W = 77,  p = 0.02286 (n = 10, 5 animals control; n = 9, 3 animals 6-OHDA)
iSPN scatter W = 30,  p = 0.01732 (n =  6, 4 animals control; n = 5, 4 animals 6-OHDA)
all tau comparisons p = 1.00000

Cells with only one fitted component are shown red (unpaired)
