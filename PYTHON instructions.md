<h1 align="center">Python Analysis Instructions</h1>

These instructions describe the Python environment used by this repository for:

- fluorescence biosensor imaging analysis for Figure 7
- RNAscope NWB conversion, metadata inspection, ROI analysis, and reanalysis for Figure 11

## Index

- [Setup Instructions](#setup-instructions)
- [Run Jupyter Notebooks](#run-jupyter-notebooks)
- [Download DANDI Data](#download-dandi-data)
- [Python Workflow Files](#python-workflow-files)
- [Export Environment](#export-environment)

## Setup Instructions

Run from the repository root:

```bash
cd "$HOME/Documents/Repositories/analysis_Belal2026"

conda create -n image_analysis python=3.12 -y
conda activate image_analysis

pip install matplotlib seaborn python-dotenv remfile
pip install notebook jupyterlab
pip install plotly
pip install -U kaleido
pip install ipykernel
pip install ipympl
pip install pytz
pip install dandi
pip install pynwb hdmf
pip install neuroconv xmltodict
pip install ndx-optogenetics
pip install openpyxl
pip install scipy

python -m ipykernel install --user --name image_analysis --display-name "Python (image_analysis)"
```

Use the `Python (image_analysis)` kernel in the notebooks.

## Run Jupyter Notebooks

```bash
cd "$HOME/Documents/Repositories/analysis_Belal2026"
conda activate image_analysis
jupyter notebook
```

The relevant notebooks are:

```text
Paper analysis/Figure 7/Figure 7 analysis.ipynb
Paper analysis/Figure 7/Figure 7 reanalysis.ipynb
RNAscope notebooks/RNAscope2NWBconversion.ipynb
RNAscope notebooks/RNAscope_analysis.ipynb
RNAscope notebooks/RNAscope reanalysis.ipynb
RNAscope notebooks/RNAscope_readNWBexample.ipynb
```

## Download DANDI Data

Figure 7 uses NWB data from DANDI dataset `001832`.

```bash
cd "$HOME/Documents/Repositories/analysis_Belal2026"
mkdir -p NWBdata

dandi download "DANDI:001832/draft" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBdata --existing error --format pyout --path-type exact

dandi download "DANDI:001832/0.XXXXX" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBdata --existing REFRESH --format PYOUT --path-type EXACT
```

Different versions of the DANDI command line tool use different casing for option values.
Use the uppercase or lowercase command according to what `dandi download --help` shows on that machine.

Use `--existing ERROR` or `--existing error` for a first clean download.
Use `--existing REFRESH` or `--existing refresh` when the files already exist and need to be checked or updated.
If the dandiset has a published version, replace `0.XXXXX` with the published version identifier.

The fluorescence biosensor notebooks expect the downloaded data under:

```text
NWBdata/001832
```

## Python Workflow Files

Core helper files:

```text
Python functions/master_functions.py
Python functions/master_RNAscope.py
```

Figure 7 fluorescence biosensor workflow:

```text
Fluorescence Biosensor Imaging instructions.md
Paper analysis/Figure 7/Figure 7 analysis.ipynb
Paper analysis/Figure 7/Figure 7 reanalysis.ipynb
```

Figure 11 RNAscope workflow:

```text
RNAscope instructions.md
RNAscope notebooks/RNAscope2NWBconversion.ipynb
RNAscope notebooks/RNAscope_analysis.ipynb
RNAscope notebooks/RNAscope reanalysis.ipynb
RNAscope notebooks/RNAscope_readNWBexample.ipynb
```

RNAscope data are stored locally under:

```text
RNAscope data/RAW
RNAscope data/NWB
```

## Export Environment

After the environment is working:

```bash
conda env export --from-history > environment.yml
conda env export > environment.lock.yml
pip freeze > requirements.txt
```
