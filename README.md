<h1 align="center">Analysis code</h1>

The code presented here was used to analyse the datasets in the following manuscript:

**Cholinergic interneuron control of intrastriatal GABAergic circuits targeting spiny projection neurons is disrupted in Parkinson’s disease models**

Belal, M.<sup>1</sup>, Perez-Rosello, T.<sup>1</sup>, Guven E. B.<sup>4</sup>, Kocaturk S.<sup>5</sup>, Xie, Z.<sup>1</sup>, Li, J.<sup>3</sup>, Dauer, W.<sup>2</sup>, Tkatch, T.<sup>1</sup>, Assous, M.<sup>5</sup>, Tepper, J.M.<sup>4</sup>, Clarke, V.R.J.<sup>1</sup>, Surmeier, D.J.<sup>1</sup>

Affiliations

<sup>1</sup> Department of Neuroscience, Feinberg School of Medicine, Northwestern University, Chicago, IL USA  
<sup>2</sup> Peter O’Donnell Jr. Brain Institute, Departments of Neurology and Neuroscience, University of Texas Southwestern Medical Center, Dallas, TX, USA  
<sup>3</sup> Department of Internal Medicine, University of Michigan Medical School. Ann Arbor, MI, USA  
<sup>4</sup> Molecular and Behavioral Neuroscience, Rutgers University, Newark, NJ USA  
<sup>5</sup> School of Biosciences, Cardiff University, Cardiff, UK

## Table of Contents

- [Initial Setup](#initial-setup)
- [Quick Start Guide](#quick-start-guide)
- [Output Structure](#output-structure)
- [macOS XQuartz Troubleshooting](#macos-xquartz-troubleshooting)
- [Fluorescence Biosensor Imaging instructions](#fluorescence-biosensor-imaging-instructions)
- [RNAscope instructions](#rnascope-instructions)
- [Non-linear Curve Fitting instructions](#non-linear-curve-fitting-instructions)
- [Technical Documentation](#technical-documentation)



First complete the Python environment setup in [PYTHON instructions.md](PYTHON%20instructions.md#setup-instructions). The `dandi` package is installed in the `image_analysis` environment there.

```bash

cd "$HOME/Documents/Repositories/analysis_Belal2026"

conda activate image_analysis

mkdir -p NWBdata

dandi download "DANDI:001832/draft" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBdata --existing error --format pyout --path-type exact

dandi download "DANDI:001832/0.XXXXX" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBdata --existing REFRESH --format PYOUT --path-type EXACT

```

On Windows, use Anaconda Prompt or Miniconda Prompt. If using PowerShell, install Anaconda or Miniconda first, accept the Conda/Anaconda Terms of Service if prompted, run `conda init powershell`, then reopen PowerShell before using `conda`. If `conda` is still not recognized, use Anaconda Prompt/Miniconda Prompt or add the Miniconda `condabin` folder to PATH. If PowerShell reports that scripts are disabled, run `Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned`.

On Windows PowerShell:

```powershell
cd "C:\Users\<USERNAME>\Documents\Repositories\analysis_Belal2026"

conda activate image_analysis

New-Item -ItemType Directory -Force -Path NWBdata

dandi download "DANDI:001832/draft" -o NWBdata --existing error --format pyout --path-type exact

dandi download "DANDI:001832/draft" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT
```

Different versions of the DANDI command line tool use different casing for option values.
If Windows cannot find `dandi` after activating `image_analysis`, run `& "$env:CONDA_PREFIX\Scripts\dandi.exe" download "DANDI:001832/draft" -o NWBdata --existing error --format pyout --path-type exact`.
Use the lowercase or uppercase command according to what `dandi download --help` shows on that machine.
For example, if dandi download --help shows --existing [ERROR|SKIP|REFRESH], use uppercase values such as --existing ERROR; if it shows lowercase values, use lowercase instead.

Use `--existing ERROR` or `--existing error` for a first clean download.
Use `--existing REFRESH` or `--existing refresh` when the files already exist and need to be checked or updated.
If the dandiset has a published version, replace `0.XXXXX` with the published version identifier.

## Initial Setup

Repository components:

- `nNLS functions.R` (core fitting functions)
- `examples/` (example input files)

Requirements:

- [R](https://cran.r-project.org/) (tested on R 4.4.x to 4.5.x)
- Packages: `robustbase`, `minpack.lm`, `Rcpp`, `signal`, `openxlsx`, `dbscan`
- C++ toolchain for `Rcpp`, required because `setup.R` compiles C++ code at load time:
- macOS: install Xcode Command Line Tools with `xcode-select --install`
- Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and make sure it is available on `PATH`; for R 4.4.x use Rtools 4.4, and for R 4.5.x or R 4.6.x use [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)
- Linux (Debian/Ubuntu): install build tools with `sudo apt-get install build-essential`
- To check the compiler from R, run `Sys.which("make")`; it should return a path, not an empty string.

## Quick Start Guide

### 1. Open R and run once

```r
rm(list = ls(all = TRUE))

load_required_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages)
  invisible(lapply(packages, library, character.only = TRUE))
}

required_packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "openxlsx", "dbscan")
load_required_packages(required_packages)
```

### 2. Source repository functions

```r
root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))
source(file.path(root_dir, 'R functions', 'setup.R'))
```

Set `ANALYSIS_ROOT` before sourcing `setup.R` if the repository is checked out somewhere else.
On Windows, set `ANALYSIS_ROOT` explicitly before sourcing `setup.R`; Rtools must be installed and available on `PATH` because `setup.R` compiles Rcpp code at load time.

### 3. Optional: R scripts that read downloaded DANDI/NWB files

The `reticulate` and `NWBenv` setup is only needed for R scripts that read downloaded DANDI/NWB files, such as `Paper analysis/*/*data processing from downloaded dandiset.R`. It is not needed just to source the core R functions.

For R scripts that use `reticulate` to read downloaded DANDI/NWB files, R may also need the Conda path written to `~/.Renviron`. If `file.edit("~/.Renviron")` does not open, create it from R with:

```r
writeLines(
  'RETICULATE_CONDA="C:/Users/<USERNAME>/miniconda3/condabin/conda.bat"',
  con = path.expand("~/.Renviron")
)
```

Replace `<USERNAME>` with the real Windows username, restart R, then check:

```r
Sys.getenv("RETICULATE_CONDA")
file.exists(Sys.getenv("RETICULATE_CONDA"))
```

The `NWBenv` Conda environment used by these R scripts must also have `pip`, `numpy`, and `pynwb` installed. The scripts create `NWBenv` if missing and install only missing packages:

```r
env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
  reticulate::conda_install(env_name, packages = c('pip', 'pynwb', 'numpy'), channel = 'conda-forge')
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()
```

### 4. Single-trace fitting example

```r
dx <- 0.1
stimulation_time <- 150
baseline <- 150
xmax <- 1000
x <- seq(dx, xmax, dx)

N <- 3
IEI <- 50

params <- c(-150, -250, -300, 1, 30, 4)
params_adj <- params
params_adj[N + 3] <- params_adj[N + 3] + stimulation_time

set.seed(42)
y <- product1N(params = params_adj, x = x, N = N, IEI = IEI) + rnorm(length(x), sd = 10)

analyse_PSC(
  response = y,
  dt = dx,
  n = 30,
  N = N,
  IEI = IEI,
  stimulation_time = stimulation_time,
  baseline = baseline,
  func = product1N,
  return.output = FALSE
)
```

## macOS XQuartz Troubleshooting

If GUI windows fail on macOS:

```bash
open -a XQuartz
export DISPLAY=:0
xhost +localhost
```

If needed permanently, add `export DISPLAY=:0` to your shell startup file.

## Technical Documentation

### [Fluorescence Biosensor Imaging instructions](Fluorescence%20Biosensor%20Imaging%20instructions.md)

This code performs fluorescence biosensor imaging analysis in Python/R, including NWB-based data loading, ROI trace analysis, response quantification, statistics, and figure generation.

### [RNAscope instructions](RNAscope%20instructions.md)

This code performs RNAscope analysis in Python/R, including NWB conversion/loading, spot-count analysis, cell-level summaries, statistics, and figure generation.

### [Non-linear Curve Fitting instructions](Non-linear%20Curve%20Fitting%20instructions.md)

This code performs non-linear fitting for postsynaptic current/potential analysis in R, including single-trace and batch workflows, UI-based fitting, and export helpers.

### [Curve Fitting Equations](Curve%20Fitting%20Equations.md)

Equations for fitting models.
