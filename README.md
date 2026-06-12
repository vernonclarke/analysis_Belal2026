<h1 align="center">Analysis code</h1>

The code presented here was used to analyse the datasets in the following manuscript:

**Cholinergic interneuron control of intrastriatal GABAergic circuits targeting spiny projection neurons is disrupted in Parkinson’s disease models**

Belal, M. <a href="https://orcid.org/0000-0001-8778-0617"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>,
Perez-Rosello, T. <a href="https://orcid.org/0009-0007-8952-2276"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>,
Guven E. B. <a href="https://orcid.org/0000-0002-9634-0485"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>4</sup>,
Kocaturk S.<sup>5</sup>,
Xie, Z. <a href="https://orcid.org/0000-0002-8348-4455"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>,
Tkatch, T. <a href="https://orcid.org/0000-0001-6626-7435"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>,
Li, J.<sup>3</sup>,
Dauer, W. <a href="https://orcid.org/0000-0003-1775-7504"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>2</sup>,
Assous, M. <a href="https://orcid.org/0000-0001-6039-816X"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>5</sup>,
Tepper, J. M. <a href="https://orcid.org/0000-0002-8643-4082"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>4</sup>,
Clarke, V. R. J. <a href="https://orcid.org/0000-0002-6154-6555"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>,
Surmeier, D. J. <a href="https://orcid.org/0000-0002-6376-5225"><img src="examples/orcid_16x16.png" width="16" height="16" alt="ORCID iD"></a> <sup>1,6</sup>

**Affiliations**

<sup>1</sup> Department of Neuroscience, Feinberg School of Medicine, Northwestern University, Chicago, IL USA  
<sup>2</sup> Peter O’Donnell Jr. Brain Institute, Departments of Neurology and Neuroscience, University of Texas Southwestern Medical Center, Dallas, TX, USA  
<sup>3</sup> Department of Internal Medicine, University of Michigan Medical School. Ann Arbor, MI, USA  
<sup>4</sup> Molecular and Behavioral Neuroscience, Rutgers University, Newark, NJ USA  
<sup>5</sup> School of Biosciences, Cardiff University, Cardiff, UK  
<sup>6</sup> Aligning Science Across Parkinson's (ASAP) Collaborative Research Network, Chevy Chase, MD 20815
## Funding

This research was funded by grants to DJS from:

Aligning Science Across Parkinson’s [ASAP020551] through the Michael J. Fox Foundation for Parkinson’s Research (MJFF); Aligning Science
Across Parkinson’s Collaborative Research Network, Chevy Chase, MD, 20815; https://parkinsonsroadmap.org.  

Freedom Together Foundation [MR-2021-2960], 875 Third Avenue, 29th Floor, New York, NY 10022; https://www.freedomtogether.org.    

National Institute of Neurological Disorders and Stroke [R37 NS034696], P.O. Box 5801. Bethesda, MD 20824; https://www.ninds.nih.gov.  

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

# Use one of the following, depending on your DANDI CLI version:

dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT

# or, if your DANDI CLI expects lowercase option values:
dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing error --format pyout --path-type exact

# or, if the files already exist and you want to refresh/check them:
dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing REFRESH --format PYOUT --path-type EXACT

# or, lowercase refresh:
dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing refresh --format pyout --path-type exact

```

On Windows, use `Anaconda Prompt` or `Miniconda Prompt`.  

If using `PowerShell`, install `Anaconda` or `Miniconda` first, accept the Conda/Anaconda Terms of Service if prompted, run `conda init powershell`, then reopen `PowerShell` before using `conda`.  

If `conda` is still not recognized or no prompt shortcut is available, see [PYTHON instructions.md](PYTHON%20instructions.md#setup-instructions) for the direct `conda.bat` setup commands for `Command Prompt` and `PowerShell`.  

If `PowerShell` reports that scripts are disabled, run `Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned`.

On Windows `PowerShell`:

```powershell
cd "C:\Users\<USERNAME>\Documents\Repositories\analysis_Belal2026"

conda activate image_analysis

New-Item -ItemType Directory -Force -Path NWBdata

# Use one of the following, depending on your DANDI CLI version:
dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing error --format pyout --path-type exact

# or, if your DANDI CLI expects uppercase option values:
dandi download "DANDI:001832/0.260611.2102" -o NWBdata --existing ERROR --format PYOUT --path-type EXACT
```

Different versions of the DANDI command line tool use different casing for option values.  

If Windows cannot find `dandi` after activating `image_analysis`, run `& "$env:CONDA_PREFIX\Scripts\dandi.exe" download "DANDI:001832/0.260611.2102" -o NWBdata --existing error --format pyout --path-type exact`. 

Use the lowercase or uppercase command according to what `dandi download --help` shows on that machine.  

For example, if dandi download --help shows --existing [ERROR|SKIP|REFRESH], use uppercase values such as --existing ERROR; if it shows lowercase values, use lowercase instead.  

Use `--existing ERROR` or `--existing error` for a first clean download.

Use `--existing REFRESH` or `--existing refresh` when the files already exist and need to be checked or updated.

## Initial Setup

Repository components:

- `nNLS functions.R` (core fitting functions)
- `examples/` (example input files)

Requirements:

- [R](https://cran.r-project.org/) (tested on R 4.4.x to 4.6.x)
- Packages loaded by `setup.R`: `robustbase`, `minpack.lm`, `Rcpp`, `signal`, `openxlsx`, `dbscan`, `shiny`, `shinybusy`, `readxl`
- C++ toolchain for `Rcpp`, required because `setup.R` compiles C++ code at load time:
- macOS: install Xcode Command Line Tools with `xcode-select --install`
- Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and make sure it is available on `PATH`; for R 4.4.x use Rtools 4.4, and for R 4.5.x or R 4.6.x use [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)
- Linux (Debian/Ubuntu): install build tools with `sudo apt-get install build-essential`
- To check the compiler from R, run `Sys.which("make")`; it should return a path, not an empty string.

## Quick Start Guide

### 1. Open R and source repository functions

```r
rm(list = ls(all = TRUE))

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))
source(file.path(root_dir, 'R functions', 'setup.R'))
```

Set `ANALYSIS_ROOT` before sourcing `setup.R` if the repository is checked out somewhere else.

On Windows, set `ANALYSIS_ROOT` explicitly before sourcing `setup.R`; Rtools must be installed and available on `PATH` because `setup.R` compiles Rcpp code at load time.

Sourcing `setup.R` installs and loads the R packages needed for both the core fitting functions and the Shiny UI functions.

On Windows, you do not need to launch R from a terminal. You can open **R x64** from the Start menu or open **RStudio**, then paste the R code above into the R console.

If you want R to work from `Command Prompt` like it does on macOS `Terminal`, first find the installed `R.exe` path:

```cmd
where /r "C:\Program Files\R" R.exe
```

This should print a path similar to:

```cmd
C:\Program Files\R\R-x.y.z\bin\x64\R.exe
```

Test that path directly:

```cmd
"C:\Program Files\R\R-x.y.z\bin\x64\R.exe" --no-save
```

In `PowerShell`, use `&` before the quoted path:

```powershell
& "C:\Program Files\R\R-x.y.z\bin\x64\R.exe" --no-save
```

Replace `R-x.y.z` with the R version folder found on your computer. To make `R.exe --no-save` work from a new `Command Prompt`, add that folder to your user PATH:

```cmd
setx PATH "%PATH%;C:\Program Files\R\R-x.y.z\bin\x64"
```

Close `Command Prompt` and reopen it, then run:

```cmd
R.exe --no-save
```

In `PowerShell`, also use:

```powershell
R.exe --no-save
```

Use `R.exe`, not `R`, because `R` can conflict with `PowerShell`'s command-history alias.

### 2. Windows note for optional R scripts that read downloaded DANDI/NWB files

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

### 3. Single-trace fitting example

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


## Contact

The analysis code was written by Vernon Clarke.

The provided code was executed on a `MacBook M2 pro 32GB` and a `Mac mini M4 Pro 64 GB`. I have tried to ensure that the code works on other operating systems but it's inevitable that some errors and bugs exist. 

If any bug fixes are necessary (most likely related to providing help on other operating systems), it will be provided as an update on the parent [`GitHub`](https://github.com/vernonclarke/analysis_Belal2026).

For queries related to this repository, please [open an issue](https://github.com/vernonclarke/analysis_Belal2026/issues) or [email](mailto:WOPR2@proton.me) directly 

