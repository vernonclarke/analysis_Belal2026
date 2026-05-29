# RNAscope Analysis Pipeline

This guide covers the full workflow in `master_RNAscope.py`: converting raw Olympus data to NWB,
loading fields, drawing ROIs, counting dots, saving results, and re-analysing saved analyses.

---

## Overview

```
Raw Olympus folder
        │
        ▼
convert_rnascope_session_to_nwb()        one-time conversion
        │
        ▼
    .nwb file
        │
        ▼
RNAscopeAnalysisStartFromSource()        load + draw ROIs
        │
        ▼
RNAscopeAnalysisFinish()                 count dots in ROIs
        │
        ▼
save_rnascope_field_analysis()           save ROIs + counts to JSON
        │
        ▼
reconstruct_state_from_saved_analysis()  reload and re-count later
```

---

## Step 1 - Convert a session folder to NWB

```python
from master_RNAscope import convert_rnascope_session_to_nwb

nwb_path, fields = convert_rnascope_session_to_nwb(
    session_dir="path/to/session_folder",
    nwb_path="path/to/output.nwb",        # optional; defaults to session_folder/session_folder.nwb
    timezone="America/Chicago",           # timezone for capture timestamps
    metadata=None,                        # optional dict - see apply_session_metadata_to_nwb
    session_description=None,             # freetext
    experiment_description=None,          # freetext
    subject_id=None,                      # e.g. "mouse_001"
    subject_sex="U",                      # "M", "F", or "U"
    subject_species="Mus musculus",
    include_exported_channel_tiffs=True,  # store per-channel display TIFFs in NWB
    overwrite=True,                       # overwrite existing NWB file
)
```

**Returns:** `(nwb_path, fields)` - path to the written NWB file and list of field dicts found.

What gets stored per field in the NWB:
- Raw 16-bit source TIFFs (`s_C00x`)
- `.pty` acquisition metadata per channel
- `.lut` display lookup tables
- Exported per-channel display TIFFs (if `include_exported_channel_tiffs=True`)

---

## Step 2 - Add / update session metadata

```python
from master_RNAscope import apply_session_metadata_to_nwb

apply_session_metadata_to_nwb(
    nwb_path="path/to/output.nwb",
    metadata={
        "NWBFile": {
            "session_description": "L1 ST8 RNAscope",
            "experimenter": ["Smith, Jane"],
            "keywords": ["RNAscope", "ChRNB2"],
        },
        "Subject": {
            "subject_id": "mouse_001",
            "species": "Mus musculus",
            "sex": "M",
            "age": "P60D",
        },
        "Custom": {
            "cohort": "L1",
            "treatment": "control",
        },
    },
)
```

---

## Step 3 - Start an analysis (draw ROIs)

### From an NWB file or folder (most common)

```python
%matplotlib widget

from master_RNAscope import RNAscopeAnalysisStartFromSource

state = RNAscopeAnalysisStartFromSource(
    source="path/to/output.nwb",            # or a session folder path
    field="fieldname",                      # e.g. "L1_ST8_s001"
    source_type=None,                       # "nwb" or "folder"; auto-detected if None
    display_mode="exported_channel_tiffs",  # or "rendered_from_raw"
    roi_specs=None,                         # see below
    roi_groups=("GB", "RB"),                # names for each ROI group
    roi_channels=("s_C002", "s_C003"),      # channel shown when drawing each group's ROIs
    count_channel="s_C004",                 # channel used for dot counting
    blind=False,                            # if True, hides the count channel while drawing
)
```

**`roi_specs`** - alternative to `roi_groups` / `roi_channels`. Pass a list of dicts:
```python
roi_specs=[
    {"group": "GB", "channel": "s_C002"},
    {"group": "RB", "channel": "s_C003"},
]
```

**`display_mode`**
- `"exported_channel_tiffs"` - uses the per-channel display TIFFs exported by Olympus (recommended)
- `"rendered_from_raw"` - reconstructs the display image from raw data + LUT metadata

After the figures appear, draw polygons in each figure window, then press **q / Escape / Enter** to confirm each ROI set.

### From pre-loaded field data

```python
from master_RNAscope import get_rnascope_field, RNAscopeAnalysisStartFromFieldData

field_data = get_rnascope_field(source="path/to/output.nwb", field="L1_ST8_s001")

state = RNAscopeAnalysisStartFromFieldData(
    field_data=field_data,
    roi_groups=("GB", "RB"),
    roi_channels=("s_C002", "s_C003"),
    count_channel="s_C004",
    blind=False,
)
```

---

## Step 4 - Finish analysis and count dots

```python
results, state = RNAscopeAnalysisFinish(
    state,
    sigma_small=1.0,            # Gaussian sigma for fine blur (DoG only)
    sigma_large=2.8,            # Gaussian sigma for coarse blur (DoG only)
    detection_method="DoG",     # "DoG" or "fiji"
    DoG_mode="tolerance",       # "tolerance" or "legacy" - only used when detection_method="DoG"
    threshold_percentile=99,    # percentile cutoff - only used when DoG_mode="legacy"
    peak_footprint=4,           # spatial footprint - only used when DoG_mode="legacy"
    maxima_tolerance=200,       # peak separation tolerance (see below)
    show_verify=True,           # show verification figure after counting
    show_detected=False,        # overlay detected dot positions on verification figure
    verify_circle_radius=2,     # radius of overlay circles (pixels)
    verify_circle_color="white",
)
```

**Returns:** `(results, state)` where `results` is a `pd.DataFrame`:

| Column | Description |
|---|---|
| `field` | Field name |
| `group` | ROI group name (e.g. "GB") |
| `image used` | Channel label used for ROI display |
| `roi` | ROI index (1-based) |
| `pixel count` | Number of pixels inside the ROI |
| `count` | Number of detected dots inside the ROI |

### Detection methods

**`detection_method="DoG"`** - subtracts a coarsely-blurred image from a finely-blurred image
to highlight spots, then finds peaks in the result.
- `DoG_mode="tolerance"` (default) - flood-fill from each peak claiming all connected pixels
  within `maxima_tolerance` intensity units (DoG scale). Higher tolerance → fewer counts.
- `DoG_mode="legacy"` - keeps only peaks in the top `threshold_percentile` of positive DoG
  values that are local maxima within a `peak_footprint` pixel window.

**`detection_method="fiji"`** - operates directly on the raw 16-bit image with no preprocessing,
matching Fiji/ImageJ's Find Maxima algorithm exactly. `maxima_tolerance` is in raw pixel
intensity units (0–65535), the same value you would type into Fiji.
Higher tolerance → fewer counts.

---

## Step 5 - Save results

```python
from master_RNAscope import save_rnascope_field_analysis

out_path = save_rnascope_field_analysis(
    state,
    results,
    analysis_params={                   # optional - records the detection settings used
        "detection_method": "DoG",
        "DoG_mode": "tolerance",
        "maxima_tolerance": 200,
    },
    out_dir=None,                       # defaults to <session_folder>/analysis/
)
# Saves to: <out_dir>/<field>.roi_analysis.json
```

The JSON stores: ROI vertices, dot counts, pixel counts, channel assignments, display mode,
and `analysis_params`.

---

## Step 6 - Reload and re-analyse

### Full reconstruct (images + ROIs)

```python
from master_RNAscope import reconstruct_state_from_saved_analysis, RNAscopeAnalysisFinish

state, analysis_payload = reconstruct_state_from_saved_analysis(
    source="path/to/output.nwb",
    analysis_path="analysis/field.roi_analysis.json",
    source_type=None,        # auto-detected
    display_mode=None,       # defaults to whatever was used originally
)

# Re-run counting with different parameters - ROIs are already loaded from the JSON
results, state = RNAscopeAnalysisFinish(
    state,
    detection_method="fiji",
    maxima_tolerance=1000,
    show_verify=True,
    show_detected=True,
)
```

### Load JSON only

```python
from master_RNAscope import load_rnascope_field_analysis

payload = load_rnascope_field_analysis("analysis/field.roi_analysis.json")
# plain dict with keys: session, field, groups, roi_specs, analysis_params, ...
```

### Apply saved ROIs to an existing state

```python
from master_RNAscope import apply_saved_analysis_to_state

state = apply_saved_analysis_to_state(state, payload)
# restores roi_sets, count_channel, roi_specs from the payload into state
```

---

## Reference - all public functions

| Function | Purpose |
|---|---|
| `convert_rnascope_session_to_nwb()` | Convert raw Olympus folder → NWB |
| `apply_session_metadata_to_nwb()` | Add/edit NWBFile, Subject, or custom metadata in an existing NWB |
| `get_rnascope_field()` | Load a single field from NWB or folder |
| `load_field_from_nwb()` | Load directly from NWB |
| `load_field_from_folder()` | Load directly from folder |
| `build_analysis_state_from_field_data()` | Build analysis state dict without opening figures |
| `RNAscopeAnalysisStartFromSource()` | Load field + open ROI drawing widgets |
| `RNAscopeAnalysisStartFromFieldData()` | Open ROI drawing widgets from pre-loaded field data |
| `RNAscopeAnalysisFinish()` | Count dots in drawn ROIs |
| `measure_rois()` | Count dots programmatically (no figures) |
| `save_rnascope_field_analysis()` | Save ROIs + counts to JSON |
| `load_rnascope_field_analysis()` | Load a saved JSON |
| `apply_saved_analysis_to_state()` | Restore ROIs from JSON into a state dict |
| `reconstruct_state_from_saved_analysis()` | Full reload: images + ROIs from saved JSON |

---

## Typical notebook workflow

```python
%matplotlib widget
from master_RNAscope import (
    convert_rnascope_session_to_nwb,
    RNAscopeAnalysisStartFromSource,
    RNAscopeAnalysisFinish,
    save_rnascope_field_analysis,
)

# 1. Convert (once per session)
nwb_path, _ = convert_rnascope_session_to_nwb("data/L1_ST8/")

# 2. Analyse one field
state = RNAscopeAnalysisStartFromSource(nwb_path, field="L1_ST8_s001")
# ... draw ROIs in the figures, press Enter when done ...

# 3. Count
results, state = RNAscopeAnalysisFinish(state, show_detected=True)
print(results)

# 4. Save
save_rnascope_field_analysis(
    state, results,
    analysis_params={
        "detection_method": "DoG",
        "DoG_mode": "tolerance",
        "maxima_tolerance": 200,
    },
)
```

