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

Install Anaconda or Miniconda before running these commands.

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
pip install tzdata
pip install dandi
pip install pynwb hdmf
pip install neuroconv xmltodict
pip install ndx-optogenetics
pip install openpyxl
pip install scipy
pip install pillow

python -m ipykernel install --user --name image_analysis --display-name "Python (image_analysis)"
```

On Windows, use Anaconda Prompt or Miniconda Prompt. If using PowerShell, install Anaconda or Miniconda first, run `conda init powershell`, then reopen PowerShell before using `conda`.

On Windows PowerShell:

```powershell
cd "C:\Users\<USERNAME>\Documents\Repositories\analysis_Belal2026"

conda create -n image_analysis python=3.12 -y
conda activate image_analysis

pip install matplotlib seaborn python-dotenv remfile
pip install notebook jupyterlab
pip install plotly
pip install -U kaleido
pip install ipykernel
pip install ipympl
pip install pytz
pip install tzdata
pip install dandi
pip install pynwb hdmf
pip install neuroconv xmltodict
pip install ndx-optogenetics
pip install openpyxl
pip install scipy
pip install pillow

$env:ANALYSIS_ROOT = "C:\Users\<USERNAME>\Documents\Repositories\analysis_Belal2026"
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

The fluorescence biosensor notebooks expect the downloaded data under:

```text
NWBdata/001832
```

## Download DANDI Data

First complete [Setup Instructions](#setup-instructions). The `dandi` package is installed in the `image_analysis` environment there.

From the repository root:

```bash
conda activate image_analysis
mkdir -p NWBdata
python -m dandi download "DANDI:001832/draft" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT
```

On Windows PowerShell:

If `conda` is not recognized, open Anaconda Prompt or Miniconda Prompt instead. To use PowerShell, install Anaconda or Miniconda first, run `conda init powershell`, then reopen PowerShell.

```powershell
conda activate image_analysis
New-Item -ItemType Directory -Force -Path NWBdata
python -m dandi download "DANDI:001832/draft" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT
```

If `python -m dandi download --help` shows lowercase option values on that machine, use lowercase values instead, such as `--existing error --format pyout --path-type exact`.

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
