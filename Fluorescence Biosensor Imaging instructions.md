# <center>Fluorescence Biosensor Imaging</center>

These instructions describe the fluorescence biosensor imaging workflow used for Figure 7 in this repository.

The workflow starts from NWB optical physiology files, computes baseline-normalized fluorescence traces, averages repeated runs for each ROI, exports the peak response table, and then runs the R statistics/plotting script used for the manuscript figure.

## Relevant Files

- `Paper analysis/Figure 7/Figure 7 analysis.ipynb`
  - Reads fluorescence traces already stored in NWB.
  - Computes dF/F0.
  - Averages repeated runs by animal, slice, and ROI.
  - Exports `results.xlsx` and representative trace files.

- `Paper analysis/Figure 7/Figure 7 reanalysis.ipynb`
  - Re-extracts fluorescence from raw NWB image frames.
  - Uses an interactively drawn ROI.
  - Computes dF/F0 and repeats the same averaging/export workflow.

- `Paper analysis/Figure 7/Figure 7 analysis and stats.R`
  - Reads the exported `results.xlsx`.
  - Fits mixed-effects models.
  - Runs robust bootstrap and Bayesian analyses.
  - Exports the Figure 7 SVG plots.

- `Python functions/master_functions.py`
  - Contains the Python helpers used by the notebooks, including:
    - `load_all_nwb`
    - `load_nwb_fluorescence`
    - `compute_dff`
    - `plot_dff`
    - `boxplot_rtype`

## Data Location

The notebooks expect the downloaded NWB data here:

```text
NWBdata/001832
```

The Figure 7 outputs are stored here:

```text
Paper analysis/Figure 7/xlsx
Paper analysis/Figure 7/svg
```

## Analysis Parameters

The notebooks use:

```python
pmt_background = 144
baseline_window = 1
xlim = [4.5, 7]
```

`compute_dff` subtracts the PMT background, computes the pre-stimulus baseline over the 1 second window before stimulation, and returns:

```text
(F - F0) / F0
```

where `F0` is the mean fluorescence during the pre-stimulus baseline window.

## Workflow 1: Read Stored NWB Fluorescence Traces

Use:

```text
Paper analysis/Figure 7/Figure 7 analysis.ipynb
```

This notebook loads the NWB fluorescence traces using:

```python
df = load_all_nwb(NWB_DIR)
```

It then selects the biosensor fluorescence series:

```python
df_grab = df[df["series_name"].str.contains("GRABCh")].reset_index(drop=True)
```

For each trial, it computes dF/F0:

```python
dff = compute_dff(
    f_raw=f_raw,
    time=t,
    stim_time=stim_time,
    pmt_background=pmt_background,
    baseline_window=baseline_window,
)
```

Repeated runs are averaged within:

```text
group
animal_id
slice_id
roi_id
```

The peak dF/F0 is then calculated for each averaged ROI:

```python
df_fl_average["fmax"] = df_fl_average["fluorescence"].apply(
    lambda f: np.asarray(f, dtype=float).max()
)
```

The exported manuscript-analysis table is:

```text
Paper analysis/Figure 7/xlsx/results.xlsx
```

with columns:

```text
Animal
Slice
ROI
Condition
(F1-F0)/F0
```

The condition labels are exported as:

```text
Control
MCI-Park
```

## Workflow 2: Re-Extract Fluorescence From Image Frames

Use:

```text
Paper analysis/Figure 7/Figure 7 reanalysis.ipynb
```

This notebook loads the raw two-photon image frames from each NWB file, draws a rectangle on the mean image, converts that rectangle to an ellipse mask, and extracts mean fluorescence across the ROI pixels.

The notebook then aligns timestamps to the NWB fluorescence timeline where available, reads stimulus timing from the trials table, computes dF/F0, averages repeated runs, and exports the same kind of peak-response table as Workflow 1.

Use this workflow when the stored NWB ROI fluorescence needs to be regenerated from the raw image data.

## Representative Traces

Representative single-ROI traces are exported from the notebook when `SAVE = True`.

Current outputs include:

```text
Paper analysis/Figure 7/svg/single_example_trace_ii8_animal1_slice3_roi1.svg
Paper analysis/Figure 7/svg/single_example_trace_ii30_animal1_slice1_roi1.svg
```

The matching trace spreadsheets are:

```text
Paper analysis/Figure 7/xlsx/single_example_animal1_slice3_roi1.xlsx
Paper analysis/Figure 7/xlsx/single_example_animal1_slice1_roi1.xlsx
```

Each spreadsheet contains:

- `summary`: animal, slice, ROI, condition, run metadata, stimulus time, and peak dF/F0.
- `trace`: time and dF/F0 values for the plotted trace.

## R Statistics And Figure Export

After `results.xlsx` has been generated, run:

```text
Paper analysis/Figure 7/Figure 7 analysis and stats.R
```

The R script reads:

```text
Paper analysis/Figure 7/xlsx/results.xlsx
```

It renames `(F1-F0)/F0` to `dff`, creates a slice-level identifier from `Animal` and `Slice`, and sets:

```r
d$Condition <- factor(d$Condition, levels = c('Control', 'MCI-Park'))
```

The main mixed-effects model is:

```r
dff ~ Condition + (1 | SliceID)
```

The script also:

- fits a robust mixed model with `rlmer`
- runs a cluster bootstrap resampling whole slices
- runs a Bayesian mixed model with `brms`
- calculates posterior and posterior-predictive probabilities
- exports manuscript SVG plots

The bootstrap defaults to one core:

```r
Sys.getenv('FIGURE7_BOOT_CORES', unset='1')
```

Set `FIGURE7_BOOT_CORES` before running the script if parallel bootstrap execution is needed.

## Figure Outputs

The R script writes SVG files to:

```text
Paper analysis/Figure 7/svg
```

Current outputs include:

```text
amplitude_SPN.svg
pseudoreplication.svg
Bayesian_Analysis.svg
fmax_boxplot.svg
single_example_trace_ii8_animal1_slice3_roi1.svg
single_example_trace_ii30_animal1_slice1_roi1.svg
```

## Minimal Run Order

1. Download or confirm the NWB data are present under `NWBdata/001832`.
2. Run `Paper analysis/Figure 7/Figure 7 analysis.ipynb` to generate `results.xlsx`.
3. If ROI fluorescence needs to be regenerated from raw image frames, use `Paper analysis/Figure 7/Figure 7 reanalysis.ipynb` instead.
4. Run `Paper analysis/Figure 7/Figure 7 analysis and stats.R` to generate the statistical results and SVG plots.
