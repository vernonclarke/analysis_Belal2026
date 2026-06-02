# Rfits: Non-Linear Curve Fitting for Postsynaptic Current Analysis

Rfits provides non-linear fitting tools for postsynaptic current/potential analysis in R, including single-trace and batch workflows, UI-based fitting, and export helpers.

## Table of Contents

- [Initial Setup](#initial-setup)
- [Quick Start Guide](#quick-start-guide)
- [Output Structure](#output-structure)
- [macOS XQuartz Troubleshooting](#macos-xquartz-troubleshooting)
- [Technical Documentation](#technical-documentation)



```bash

cd "$HOME/Documents/Repositories/analysis_Belal2026"

mkdir -p NWBdata

dandi download "DANDI:001832/draft" -o NWBData --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBData --existing error --format pyout --path-type exact

dandi download "DANDI:001832/0.XXXXX" -o NWBData --existing ERROR --format PYOUT --path-type EXACT

dandi download "DANDI:001832/draft" -o NWBData --existing REFRESH --format PYOUT --path-type EXACT

```



## Initial Setup

Repository components:

- `nNLS functions.R` (core fitting functions)
- `analyseResponse.R` (modular Shiny workflow)
- `examples/` (example input files)
- `images/` (plot assets)

Requirements:

- R (tested on R 4.4.x to 4.5.x)
- Packages: `robustbase`, `minpack.lm`, `Rcpp`, `signal`, `writexl`
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

required_packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "writexl")
load_required_packages(required_packages)
```

### 2. Source repository functions

```r
UserName <- "YourUserName"  # replace with your macOS username
root_dir <- file.path("/Users", UserName, "Documents", "Repositories", "Rfits")
source(file.path(root_dir, "nNLS functions.R"))
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

Use `TECHNICAL_GUIDE.md` for equations, model details, fitting strategy, and AIC/BIC conventions.
