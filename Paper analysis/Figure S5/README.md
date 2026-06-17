# Figure S5. Reduced ChI-evoked GABAergic input after 6-OHDA (independent study)

Independent replication of the 6-OHDA result using a different model: ChAT-ChR2-eYFP mice with a
unilateral 6-OHDA lesion of the substantia nigra pars compacta (SNc). Optogenetic ChI activation
evokes smaller GABAergic PSCs in SPNs after lesion (reduced amplitude and charge transfer),
with unchanged decay kinetics. Comparison unpaired (control vs 6-OHDA)

## Scripts (run order)
1. `Figure S5 data processing from downloaded dandiset.R` extracts/averages PSCs to per-condition `xlsx`
2. `Figure S5 ctrl.R` (control) and `Figure S5 6OHDA.R` (lesion) perform two-component fits
3. `Figure S5 graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`ctrl.RData` and `6OHDA.RData` are generated R workspace files produced by `Figure S5 ctrl.R` and `Figure S5 6OHDA.R` after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `ctrl.xlsx` / `ctrl full.xlsx`, `6OHDA.xlsx` / `6OHDA full.xlsx` per-cell averaged traces
- `Figure S5 ctrl.csv`, `Figure S5 OHDA.csv` NLS 2 PSC component fit summaries
- `Figure S5 amplitude.csv`, `Figure S5 charge transfer.csv` peak amplitude / charge transfer
- `Figure S5 amplitude 2 components.csv`, `Figure S5 tau 2 components.csv` fast/slow amplitudes and time constants
- `Figure S5 ctrl single examples.csv`, `Figure S5 OHDA single examples.csv` example traces (panels b/c)
- `Figure S5 statistics.csv` Wilcoxon rank-sum results
- `Figure S5.xlsx` combined workbook

## Panels `svg/`
- `egs_*` (control / 6-OHDA, + semilog) to **panels b, c**
- amplitude / charge-transfer box plots (+ semilog) to **panel d**
- scatter (slow vs fast amplitude) to **panel e**
- tau box plots to **panel f**
- per-cell fit-diagnostic SVGs; `dbscan_*` clustering QC

## Statistics
Wilcoxon rank-sum test
control vs 6-OHDA:

amplitude W = 39, p = 6.3 x 10-7 (n = 21 control, 9 animals; n = 21 6-OHDA, 9 animals)
charge transfer W = 404, p = 4.3 x 10-7 (n = 21 control, 9 animals; n = 21 6-OHDA, 9 animals)
fast tau W = 240, p = 0.31292 (n = 21 control, 9 animals; n = 21 6-OHDA, 9 animals)
slow tau W = 188, p = 0.76828 (n = 21 control, 9 animals; n = 21 6-OHDA, 9 animals)

Cells with only one fitted component are shown red (unpaired)
