# Figure 3 (a-f). GABA_A δ-subunit knockdown selectively reduces the slow component

CRISPR knockdown of the GABA_A δ-subunit (vs scrambled control) reduces the total charge
transfer of the ChI-evoked PSC by selectively attenuating the **slow-decaying** component,
without changing peak amplitude or decay kinetics. Comparison is unpaired
(`CRISPR control` vs `CRISPR delta KD`)

> Panels g-l of Figure 3 (NGF-interneuron-evoked currents) are in the **`Figure 3B`** folder

## Scripts (run order)
1. `Figure 3 CRISPR data processing from downloaded dandiset.R` extracts/averages PSCs to per-condition `xlsx`
2. `Figure 3 CRISPR ctrl analysis.R`, `Figure 3 CRISPR KD analysis.R` two-component fits per condition
3. `Figure 3 CRISPR graphs.R` plots SVG panels and runs unpaired Wilcoxon tests

`CRISPR control.RData` and `CRISPR delta KD.RData` are generated R workspace files produced by `Figure 3 CRISPR ctrl analysis.R` and `Figure 3 CRISPR KD analysis.R` after fitting the averaged traces; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `Figure 3 amplitude.csv`, `Figure 3 charge transfer.csv` peak amplitude / total charge transfer
- `Figure 3 amplitude 2 components.csv`, `Figure 3 area 2 components.csv` *(if present)* fast/slow component amplitudes and areas
- `Figure 3 tau 2 components.csv` fast/slow decay time constants
- `Figure 3 * single examples.csv` example traces (panels b/c)
- `Figure 3 statistics.csv` Wilcoxon rank-sum results
- `Figure 3.xlsx` combined workbook

## Panels `svg/`
- `egs_*` (control / KD, + semilog) to **panels b, c**
- `amplitude_*`, `charge_transfer_*` (+ semilog) to **panel d**
- `fast_area_2component_SPN.svg`, `slow_area_2component_SPN.svg` (+ semilog) to **panel e** (selective slow reduction)
- `tau_*` to **panel f**
- `CRISPR control_*.svg`, `CRISPR delta KD_*.svg` per-cell fit diagnostics
- `dbscan_*.svg` clustering QC

## Statistics
Wilcoxon rank-sum test

control vs delta-KD:

amplitude W = 19, p = 0.13307 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
charge transfer W = 60, p = 0.01357 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
fast component W = 23, p = 0.26985 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
slow component W = 62, p = 0.01357 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
fast tau W = 18, p = 0.21761 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
slow tau W = 40, p = 0.66907 (n = 7 control, 6 animals; n = 10 delta-KD, 6 animals)
