# Figure 3 (g-l). NGF-interneuron-evoked GABAergic currents in SPNs

Optogenetic activation (whole-field LED, 5 ms) of NDNF⁺ NGF interneurons evokes a GABAergic
PSC in SPNs that is kinetically similar to the ChI-evoked current: a dominant slow-decaying
component on top of a smaller fast component. Single experimental condition (no drug/lesion
comparison). This is the `Figure 3B` working folder for manuscript panels **3g-3l**

## Scripts (run order)
1. `Figure 3B data processing from downloaded dandiset.R` extracts/averages PSCs to `xlsx/NDNF GABA PSCs.xlsx`
2. `Figure 3B analysis.R` two-component fits to summary CSVs + `.RData`
3. `Figure 3B graphs.R` plots SVG panels (no group comparison)

`NDNF GABA PSCs.RData` is a generated R workspace file produced by `Figure 3B analysis.R` after fitting the averaged traces; it is not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
- `NDNF GABA PSCs.xlsx`, `NDNF GABA PSCs_full.xlsx` per-cell averaged (and un-subtracted) traces
- `Figure 3B NDNF_GABA_PSCs.csv` NLS 2 PSC component fit summaries
- `Figure 3B amplitude.csv`, `Figure 3B charge transfer.csv` peak amplitude / charge transfer
- `Figure 3B amplitude 2 components.csv`, `Figure 3B tau 2 components.csv` fast/slow amplitudes and time constants
- `Figure 3B SPN single examples.csv` example trace (panel i)
- `Figure 3B.xlsx` combined workbook

## Panels `svg/`
- example-trace SVG (+ semilog) to **panel i**
- amplitude / charge-transfer box plots to **panel j**
- scatter (slow vs fast amplitude) to **panel k**
- tau box plots to **panel l**
- per-cell fit-diagnostic SVGs

## Note
In one recording only a slow (no fast) component met the two-component criterion; it is shown
in red with no paired fast value (see manuscript Methods)

## Statistics
No inferential group comparison was run for this folder
single condition:
NDNF GABA PSCs only; no p values reported
