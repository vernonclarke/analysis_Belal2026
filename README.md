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

## Initial Setup

Repository components:

- `nNLS functions.R` (core fitting functions)
- `examples/` (example input files)
- `images/` (plot assets)

Requirements:

- R (tested on R 4.4.x to 4.5.x)
- Packages: `robustbase`, `minpack.lm`, `Rcpp`, `signal`, `openxlsx`, `dbscan`
- C++ toolchain for `Rcpp`:
- macOS: Xcode Command Line Tools
- Windows: Rtools
- Linux (Debian/Ubuntu): build-essential

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
UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')
source(file.path(root_dir, 'R functions', 'setup.R'))
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

### 4. Dataset quick load example

```r
csv_path <- file.path(root_dir, "examples", "data.csv")
traces <- as.matrix(read.csv(csv_path, check.names = FALSE))

dim(traces)
head(traces[, 1, drop = FALSE])
```

## Output Structure

Typical per-trace output fields include:

- `output`
- `fits`
- `fits.se`
- `gof`
- `AIC`
- `BIC`
- `model.message`
- `sign`
- `traces`
- `fit_results` (when available)

Example:

```r
names(out_list[[1]])
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

### [Curve Fitting Equations.md](Curve-Fitting-Equations.md) 

Equations for fitting models.
