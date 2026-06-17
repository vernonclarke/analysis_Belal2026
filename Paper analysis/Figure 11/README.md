# Figure 11 (+ Figure S6). Reduced ChI input to NGF interneurons after 6-OHDA, and reduced β2-nAChR expression

This folder contains **two** related analyses:

1. **RNAscope in-situ hybridisation** of *Chrnb2* (β2 nAChR subunit) mRNA:
   - in **NDNF⁺ NGF interneurons** to **Figure 11a-b**
   - in **TH⁺ neurons** to **Figure S6** (same pipeline, `TH` file variants)
   Both show reduced β2-nAChR expression after 6-OHDA
2. **Patch-clamp connectivity**: optogenetic ChI activation (ChAT-Cre × NDNF-Flp) evokes
   GABAergic PSCs recorded **from NGF interneurons**; charge transfer is reduced after 6-OHDA
   while peak amplitude is largely preserved to **Figure 11c-f**

## Scripts
- `RNAscope_analysis.R` RNAscope spot-count analysis: robust linear mixed-effects model with
  cluster bootstrap + Bayesian model (see `RNAscope instructions.md`). Produces the NDNF (Fig 11)
  and TH (Fig S6) outputs
- `Figure 11 data processing from downloaded dandiset.R` to `Figure 11 analysis.R` to `Figure 11 graphs.R`
  runs the patch-clamp pipeline for the ChI to NGF recordings (panels c-f)

`Figure 11.RData`, `RNAscope_NDNF.RData` and `RNAscope_TH.RData` are generated R workspace files produced by `Figure 11 analysis.R` and `RNAscope_analysis.R` after the patch-clamp fitting and RNAscope analyses; they are not included in this processed-data deposit

## Processed data spreadsheets `xlsx/` or `csv/`
**RNAscope (NDNF = Fig 11, TH = Fig S6):**
- `Figure 11 NDNF RNAscope NDNF count data.csv`, `Figure 11 TH RNAscope TH count data.csv` per-field spot counts with metadata (`condition` Intact/6-OHDA, `cell_type`, `slice_id`, `field`, `hemisphere`, `Animal`, `Group`, `count`, …)
- `... NDNF lme fixed effect.csv`, `... TH lme fixed effect.csv` mixed-effects fixed-effect estimates
- `... NDNF robust lme bootstrap.csv`, `... TH robust lme bootstrap.csv` cluster-bootstrap p-values (N_boot = 1,999)
- `... NDNF bayes summary.csv`, `... TH bayes summary.csv` Bayesian posterior summaries
- `Figure 11 NDNF RNAscope.xlsx`, `Figure 11 TH RNAscope.xlsx` combined workbooks
- `RNAscope_experimenter.csv`, `RNAscope_experimenter_counts.csv` inter-rater / manual-count comparison data

**Patch clamp (Fig 11c-f):**
- `ChI-NGF control.xlsx`, `ChI-NGF 6OHDA.xlsx` per-cell averaged NGF-interneuron traces
- `Figure 11 amplitude.csv`, `Figure 11 charge transfer.csv` peak amplitude / charge transfer
- `Figure 11 ChI-NGF single examples.csv`, `Figure 11 ChI-NGF 6OHDA single examples.csv` example traces (panels d/e)
- `Figure 11 statistics.csv` Wilcoxon rank-sum results
- `Figure 11.xlsx` combined workbook

## Panels `svg/`
**RNAscope:** `Figure 11_RNAscope_boxplot_NDNF.svg` (R) / `..._NDNF_PYTHON.svg` (Python replica) to **Fig 11b**;
`..._TH.svg` / `..._TH_PYTHON.svg` to **Fig S6b**; `Figure 11_Bayesian_Analysis_NDNF.svg` / `..._TH.svg` Bayesian
posterior plots; `NDNF_user_counts_boxplot.svg`, `TH_user_counts_boxplot.svg` inter-rater QC
**Patch clamp:** `egs_ctrl.svg`, `egs_OHDA.svg` (+ semilog) to **panels d, e**; `amplitude.svg` /
`semilog_amplitude_SPN.svg`, `charge_transfer.svg` / `semilog_charge_transfer_SPN.svg` to **panel f**;
`ChI_NGF_control_<n>.svg`, `ChI_NGF_6OHDA_<n>.svg` per-cell fit diagnostics

## Statistics
Robust linear mixed-effects model with cluster bootstrap
control vs 6-OHDA:
RNAscope NDNF fixed effect p = 0.00400 (N_boot = 1,999)
RNAscope TH   fixed effect p < 0.0005 (N_boot = 1,999)

Wilcoxon rank-sum test
control vs 6-OHDA ChI-NGF PSC:
NGF interneuron amplitude W = 17, p = 0.05556       (n = 7 control, 4 animals; n = 11 6-OHDA, 3 animals)
NGF interneuron charge transfer W = 68, p = 0.00591 (n = 7 control, 4 animals; n = 11 6-OHDA, 3 animals)
 
