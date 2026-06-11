<h1 align="center">RNAscope Imaging Analysis</h1>

These instructions describe the RNAscope workflow used for Figure 11 in this repository.

The workflow has two related paths:

- **Experimenter/manual counts**: the experimenter drew ROIs and counted CHRNB2 puncta by eye. These counts are stored in the RNAscope NWB files and are the default data used by the Figure 11 R analysis.
- **Automated/reanalysis counts**: saved ROI JSON files can be re-counted with the Python RNAscope analysis functions and compared with the experimenter counts.

## Relevant Files

- `Python functions/master_RNAscope.py`
  - RNAscope conversion, NWB loading, ROI drawing, dot counting, JSON saving, and count reconstruction helpers.

- `RNAscope notebooks/RNAscope2NWBconversion.ipynb`
  - Converts raw Olympus RNAscope folders to NWB.
  - Stores metadata and experimenter/manual CHRNB2 counts.

- `RNAscope notebooks/RNAscope_analysis.ipynb`
  - Loads one field from an NWB file.
  - Draws ROIs for NDNF+ and TH+ cells.
  - Counts CHRNB2 puncta in the count channel.
  - Optionally saves one `.roi_analysis.json` file.

- `RNAscope notebooks/RNAscope reanalysis.ipynb`
  - Batch reloads saved ROI JSON files.
  - Re-counts the ROIs with new analysis parameters.
  - Compares automated/reanalysis counts with the experimenter counts stored in NWB.

- `RNAscope notebooks/RNAscope_readNWBexample.ipynb`
  - Reads metadata and stored experimenter counts from the NWB files.
  - Plots NDNF+ and TH+ experimenter counts.

- `Paper analysis/Figure 11/RNAscope_analysis.R`
  - Runs the Figure 11 RNAscope statistics and plots.

## Data Locations

Raw Olympus RNAscope data are stored here:

```text
RNAscope data/RAW
```

Converted NWB files are stored here:

```text
RNAscope data/NWB
```

Current sessions:

```text
L1.ST8
L2.ST8
L3.ST6
L4.ST8
```

Each converted session has:

```text
RNAscope data/NWB/<session>/<session>.nwb
```

Downloaded DANDI-style NWB files are stored separately under:

```text
NWBdata/001832/sub-L1-ST8/sub-L1-ST8_ses-20240905T115902.nwb
```

Saved ROI analyses are stored under:

```text
NWBdata/001832/sub-<subject>/analysis
```

Batch reanalysis JSON files are stored under:

```text
NWBdata/001832/sub-<subject>/reanalysis
```

Figure 11 R outputs are stored here:

```text
Paper analysis/Figure 11/xlsx
Paper analysis/Figure 11/svg
```

## Convert Raw RNAscope Data To NWB

Use:

```text
RNAscope notebooks/RNAscope2NWBconversion.ipynb
```

The conversion uses `convert_rnascope_session_to_nwb` from `master_RNAscope.py`.

Minimal session conversion:

```python
from pathlib import Path
import os
import sys

repo = Path(os.environ.get(
    "ANALYSIS_ROOT",
    Path.home() / "Documents" / "Repositories" / "analysis_Belal2026"
))
sys.path.insert(0, str(repo / "Python functions"))

from master_RNAscope import convert_rnascope_session_to_nwb

raw_root = repo / "RNAscope data" / "RAW"
nwb_root = repo / "RNAscope data" / "NWB"

for session in ["L1.ST8", "L2.ST8", "L3.ST6", "L4.ST8"]:
    session_out = nwb_root / session
    session_out.mkdir(parents=True, exist_ok=True)

    nwb_path, fields = convert_rnascope_session_to_nwb(
        session_dir=raw_root / session,
        nwb_path=session_out / f"{session}.nwb",
        timezone="America/Chicago",
        metadata={
            "NWBFile": {
                "session_description": f"{session} RNAscope imaging session",
                "experimenter": ["Zhong, Xie"],
                "lab": "Surmeier Lab",
                "institution": "Northwestern University",
                "keywords": ["RNAscope", "CHRNB2", "NDNF", "TH"],
            },
            "Subject": {
                "subject_id": session,
                "species": "Mus musculus",
                "sex": "U",
            },
            "Custom": {
                "source_repo": "analysis_Belal2026",
                "source_raw_folder": str(raw_root / session),
            },
        },
        overwrite=True,
    )

    print(f"Wrote {nwb_path} with {len(fields)} fields")
```

Per field, the NWB stores:

- raw 16-bit source TIFFs
- `.pty` acquisition metadata
- `.lut` display lookup tables
- exported per-channel display TIFFs, when available
- field metadata such as field name, side, field index, pixel size, and channel labels

## Stored Experimenter Counts

The experimenter/manual CHRNB2 counts are stored in each NWB file at:

```python
nwbfile.processing["rnascope_analysis_metadata"]["experimenter_chrnb2_counts"]
```

The stored count table has this shape:

```text
condition
cell_type
slice_id
field
hemisphere
field_index
replicate
count
session
session_group
```

Read all stored experimenter counts:

```python
from pathlib import Path
import os
import pandas as pd
from pynwb import NWBHDF5IO

repo = Path(os.environ.get(
    "ANALYSIS_ROOT",
    Path.home() / "Documents" / "Repositories" / "analysis_Belal2026"
))
nwb_root = repo / "NWBdata" / "001832"
sessions = ["L1.ST8", "L2.ST8", "L3.ST6", "L4.ST8"]

session_paths = {
    "L1.ST8": {
        "nwb_path": nwb_root / "sub-L1-ST8" / "sub-L1-ST8_ses-20240905T115902.nwb",
    },
    "L2.ST8": {
        "nwb_path": nwb_root / "sub-L2-ST8" / "sub-L2-ST8_ses-20240909T100245.nwb",
    },
    "L3.ST6": {
        "nwb_path": nwb_root / "sub-L3-ST6" / "sub-L3-ST6_ses-20240909T112846.nwb",
    },
    "L4.ST8": {
        "nwb_path": nwb_root / "sub-L4-ST8" / "sub-L4-ST8_ses-20240916T101347.nwb",
    },
}

dfs = []

for session in sessions:
    nwb_path = session_paths[session]["nwb_path"]

    with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
        nwbfile = io.read()
        df = (
            nwbfile.processing["rnascope_analysis_metadata"]["experimenter_chrnb2_counts"]
            .to_dataframe()
            .reset_index(drop=True)
        )

    dfs.append(df)

all_user_counts = (
    pd.concat(dfs, ignore_index=True)
    .sort_values(
        ["session", "condition", "cell_type", "hemisphere", "field_index", "replicate"],
        kind="stable",
    )
    .reset_index(drop=True)
)
```

This is the count table used by default in:

```text
Paper analysis/Figure 11/RNAscope_analysis.R
```

with:

```r
reanalysis = FALSE
RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_experimenter')
```

## Analyse One Field Manually

Use:

```text
RNAscope notebooks/RNAscope_analysis.ipynb
```

Example setup:

```python
from pathlib import Path
import sys

repo = Path.home() / "Documents" / "Repositories" / "analysis_Belal2026"
sys.path.insert(0, str(repo / "Python functions"))

from master_RNAscope import (
    get_rnascope_field,
    RNAscopeAnalysisStartFromFieldData,
    RNAscopeAnalysisFinish,
    save_rnascope_field_analysis,
)

session = "L1.ST8"
field = "L1.ST8_C.L_60x.01"

nwb_root = repo / "NWBdata" / "001832"
session_paths = {
    "L1.ST8": {
        "nwb_path": nwb_root / "sub-L1-ST8" / "sub-L1-ST8_ses-20240905T115902.nwb",
    },
    "L2.ST8": {
        "nwb_path": nwb_root / "sub-L2-ST8" / "sub-L2-ST8_ses-20240909T100245.nwb",
    },
    "L3.ST6": {
        "nwb_path": nwb_root / "sub-L3-ST6" / "sub-L3-ST6_ses-20240909T112846.nwb",
    },
    "L4.ST8": {
        "nwb_path": nwb_root / "sub-L4-ST8" / "sub-L4-ST8_ses-20240916T101347.nwb",
    },
}

nwb_path = session_paths[session]["nwb_path"]
```

Current channel setup:

```python
display_mode = "rendered_from_raw"
roi_specs = (
    ("NDNF+", "s_C003"),
    ("TH+", "s_C002"),
)
count_channel = "s_C004"
blind = True
```

Load the field:

```python
field_data = get_rnascope_field(
    nwb_path,
    field=field,
    display_mode=display_mode,
)
```

Draw ROIs:

```python
%matplotlib widget

state = RNAscopeAnalysisStartFromFieldData(
    field_data,
    roi_specs=roi_specs,
    count_channel=count_channel,
    blind=blind,
)
```

Draw polygons in the widget figures. Press `q`, `escape`, or `enter` inside each figure when finished.

## Count CHRNB2 Puncta

Current one-field analysis parameters:

```python
analysis_params = {
    "detection_method": "DoG",
    "sigma_small": 1.0,
    "sigma_large": 2.8,
    "threshold_percentile": 99.9,
    "peak_footprint": 4,
    "maxima_tolerance": 170,
    "show_detected": True,
    "show_verify": True,
}
```

Run counting:

```python
results, state = RNAscopeAnalysisFinish(
    state,
    **analysis_params,
    verify_image="count",
)

display(results)
```

`verify_image="count"` displays the count channel during verification. `verify_image="group"` displays the ROI-group channel.

The result table includes:

```text
field
group
image used
roi
pixel count
count
```

Save one field analysis:

```python
SAVE = True

if SAVE:
    json_path = save_rnascope_field_analysis(
        state,
        results,
        analysis_params=analysis_params,
    )
    print(json_path)
```

This writes:

```text
NWBdata/001832/sub-<subject>/analysis/<field>.roi_analysis.json
```

The JSON stores ROI vertices, count settings, dot counts, pixel counts, channel assignments, display mode, and the saved analysis parameters.

## Detection Methods

`detection_method="DoG"` subtracts a coarsely blurred image from a finely blurred image to enhance puncta before peak detection.

`dog_mode="tolerance"` uses a tolerance-style local-maxima procedure on the DoG image. Higher `maxima_tolerance` generally gives fewer detected puncta.

`dog_mode="legacy"` uses `threshold_percentile` and `peak_footprint` to keep positive DoG local maxima.

`detection_method="fiji"` uses a Fiji/ImageJ Find Maxima-style local-maxima procedure on the raw 16-bit count image. This is Fiji-like, but the exact count can still depend on image scaling, channel choice, ROI definition, and tolerance.

For both methods, increasing `maxima_tolerance` reduces over-counting by merging or rejecting less prominent nearby maxima.

## Batch Reanalysis

Use:

```text
RNAscope notebooks/RNAscope reanalysis.ipynb
```

This notebook:

1. Reads saved ROI JSON files from `NWBdata/001832/sub-<subject>/analysis`.
2. Reconstructs the analysis state from each JSON and NWB file.
3. Re-runs `RNAscopeAnalysisFinish` with new parameters.
4. Saves reanalysis JSON files to `NWBdata/001832/sub-<subject>/reanalysis`.
5. Reconstructs a count table from the reanalysis JSON files.
6. Compares reanalysis counts with the experimenter counts stored in NWB.

Current batch settings:

```python
run_reanalysis = True
save_reanalysis = True
source_analysis_folder = "analysis"
reanalysis_folder = "reanalysis"
display_mode = "rendered_from_raw"

detection_method = "DoG"
dog_mode = "tolerance"

new_params = {
    "detection_method": detection_method,
    "dog_mode": dog_mode,
    "sigma_small": 1.0,
    "sigma_large": 2.8,
    "threshold_percentile": 99,
    "peak_footprint": 4,
    "maxima_tolerance": 170,
    "show_detected": False,
    "show_verify": False,
}
```

Reconstruct count rows from saved JSON files:

```python
from master_RNAscope import counts_from_roi_jsons

reanalysis_counts = counts_from_roi_jsons(reanalysis_json_paths)
```

The reconstructed count table includes:

```text
condition
cell_type
slice_id
field
hemisphere
field_index
replicate
count
session
session_group
count_channel
analysis_json
```

## Figure 11 R Analysis

Run:

```text
Paper analysis/Figure 11/RNAscope_analysis.R
```

Choose:

```r
cell_type <- 'NDNF'
```

or:

```r
cell_type <- 'TH'
```

The R script loads experimenter counts by default:

```r
reanalysis = FALSE
RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_experimenter')
```

To use automated reanalysis counts instead:

```r
reanalysis = TRUE
RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_reanalysed')
```

The script maps:

```text
UL -> Control
L  -> 6OHDA
```

Python count tables store `condition` as `Intact` or `Lesioned`; the Figure 11 R script relabels these groups from `hemisphere` for plotting and statistics.

The main mixed-effects model is:

```r
var ~ Group + (1 | Animal/field)
```

The script also runs robust bootstrap and Bayesian analyses, then exports CSV/XLSX summaries and SVG plots to:

```text
Paper analysis/Figure 11/xlsx
Paper analysis/Figure 11/svg
```

## Function Reference

| Function | Purpose |
|---|---|
| `convert_rnascope_session_to_nwb()` | Convert one raw Olympus RNAscope session folder to NWB |
| `apply_session_metadata_to_nwb()` | Add or update NWBFile, Subject, and custom metadata |
| `get_rnascope_field()` | Load one field from an NWB file or raw folder |
| `load_field_from_nwb()` | Load one field directly from NWB |
| `load_field_from_folder()` | Load one field directly from a raw folder |
| `build_analysis_state_from_field_data()` | Build an analysis state dict without opening figures |
| `RNAscopeAnalysisStartFromSource()` | Load one field and open ROI drawing widgets |
| `RNAscopeAnalysisStartFromFieldData()` | Open ROI drawing widgets from pre-loaded field data |
| `RNAscopeAnalysisFinish()` | Count puncta in the drawn ROIs |
| `measure_rois()` | Count puncta programmatically without opening figures |
| `save_rnascope_field_analysis()` | Save ROIs and counts to JSON |
| `load_rnascope_field_analysis()` | Load one saved ROI/count JSON |
| `load_session_analysis_jsons()` | Load all saved ROI/count JSON files from one session |
| `apply_saved_analysis_to_state()` | Restore saved ROIs into an analysis state |
| `reconstruct_state_from_saved_analysis()` | Reload images and saved ROIs from an analysis JSON |
| `counts_from_roi_jsons()` | Reconstruct a tidy count table from saved ROI JSON files |
