# Figure 1. Gabazine reduces ChI-evoked PSCs in SPNs

Bath gabazine (10 µM) blocks the optogenetically (Chronos, whole-field LED, 5 ms) ChI-evoked
postsynaptic current in both direct- (dSPN) and indirect-pathway (iSPN) spiny projection
neurons, confirming the response is largely GABAA-receptor mediated. Comparison is **paired** (control
vs gabazine, within each cell).

## Scripts (run order)
1. `Figure 1 data processing from downloaded dandiset.R` extracts and averages the PSCs from
   the DANDI NWB files to `xlsx/ChAT-Cre X tdTomato.xlsx` (dSPN) and `xlsx/ChAT-Cre X De eGFP.xlsx` (iSPN)
2. `Figure 1 analysis.R` fits each averaged trace, writes summary CSVs and `Figure 1.RData`
3. `Figure 1 graphs.R` plots the SVG panels and runs the paired Wilcoxon tests

`Figure 1.RData` is a generated R workspace file produced by `Figure 1 analysis.R` after fitting the averaged traces; it is not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `ChAT-Cre X tdTomato.xlsx`, `ChAT-Cre X De eGFP.xlsx` per-cell averaged gabazine-sensitive traces (dSPN, iSPN)
- `Figure 1 amplitude.csv`, `Figure 1 charge transfer.csv` per-cell peak amplitude / charge transfer; `condition` = `ctrl` / `GABAzine`
- `Figure 1 dSPN single eg.csv`, `Figure 1 iSPN single eg.csv` the example traces shown in panel d (columns `x, y, yfilter, yfit, yfit1, yfit2`)
- `Figure 1 statistics.csv` Wilcoxon **signed-rank** (paired) results
- `Figure 1.xlsx` complete summary workbook

## Figure panels `svg/`
- `egs_dSPN.svg`, `egs_iSPN.svg` (+ `egs_semilog_*`) to **panel d** (representative dSPN/iSPN PSC, control vs gabazine)
- `amplitude_SPN.svg`, `amplitude_semilog_SPN.svg` to **panel e** (amplitude box plots)
- `charge_transfer_SPN.svg`, `charge_transfer_semilog_SPN.svg` to **panel f** (charge transfer box plots)
- `dbscan_dSPN_svg`, `dbscan_iSPN_svg` clustering QC (not a published panel)

## Statistics
Wilcoxon signed-rank test
control vs gabazine:
dSPN amplitude V = 0, p = 0.00391 (n = 10, 6 animals)
iSPN amplitude V = 0, p = 0.03125 (n = 6, 6 animals)
dSPN charge transfer V = 55, p = 0.00391 (n = 10, 6 animals)
iSPN charge transfer V = 21, p = 0.03125 (n = 6, 6 animals)
