# Figure 2. ChI-evoked PSCs are similar in iSPNs and dSPNs

The gabazine-sensitive ChI-evoked PSC is decomposed into fast- and slow-decaying components
and compared **between** cell types (dSPN vs iSPN, unpaired). Amplitude, charge transfer and
decay kinetics do not differ between pathways; the slow-decaying component dominates in both

## Scripts (run order)
1. `Figure 2 ChAT-Cre X tdTomato.R` (dSPN) and `Figure 2 ChAT-Cre X De eGFP.R` (iSPN) load the
   per-cell averaged traces (shared with Figure 1) and perform two-component fits
2. `Figure 2 graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`ChAT-Cre X tdTomato.RData` and `ChAT-Cre X De eGFP.RData` are generated R workspace files produced by `Figure 2 ChAT-Cre X tdTomato.R` and `Figure 2 ChAT-Cre X De eGFP.R` after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `ChAT-Cre X tdTomato.xlsx`, `ChAT-Cre X De eGFP.xlsx` per-cell averaged traces
- `Figure 2 amplitude.csv`, `Figure 2 charge transfer.csv` peak amplitude / charge transfer per cell
- `Figure 2 amplitude 2 components.csv` fast/slow component amplitudes
- `Figure 2 tau 2 components.csv` fast/slow decay time constants
- `Figure 2 dSPN single examples.csv`, `Figure 2 iSPN single examples.csv` example traces (panels b/c)
- `Figure 2 statistics.csv` Wilcoxon rank-sum results
- `Figure 2.xlsx` combined workbook

## Panels `svg/`
- `egs_dSPN.svg`, `egs_iSPN.svg` (+ `semilog_egs_*`) to **panels b, c** (representative PSCs with fast/slow fit; inset = semi-log)
- `amplitude_SPN.svg` / `*_semilog_*`, `charge_transfer_SPN.svg` / `*_semilog_*` to **panel d**
- `scatter.svg`, `scatter_charge_transfer.svg` to **panel e** (slow- vs fast-decaying amplitude)
- `tau_fast_SPN.svg`, `tau_slow_SPN.svg`, `tau_2component_dSPN.svg` to **panel f** (decay time constants)
- `amplitude_2component_SPN.svg`, `area_2component_SPN.svg` fast/slow component box plots
- `ChAT-Cre X tdTomato_*.svg`, `ChAT-Cre X De eGFP_*.svg` per-cell fit diagnostics
- `dbscan_*.svg` clustering QC

## Statistics
Wilcoxon rank-sum test
dSPN vs iSPN:
amplitude W = 42, p = 0.21978 (n = 10 dSPN, 6 animals; n = 6 iSPN, 6 animals)
charge transfer W = 18, p = 0.21978 (n = 10 dSPN, 6 animals; n = 6 iSPN, 6 animals)
fast tau W = 31, p = 1.00000 (n = 10 dSPN, 6 animals; n = 6 iSPN, 6 animals)
slow tau W = 30, p = 1.00000 (n = 10 dSPN, 6 animals; n = 6 iSPN, 6 animals)
