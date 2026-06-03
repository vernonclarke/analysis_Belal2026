script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

read_file <- function(path) {
  readLines(file.path(repo, path), warn = FALSE)
}

contains <- function(lines, pattern) {
  any(grepl(pattern, lines, fixed = TRUE))
}

setup <- read_file("R functions/setup.R")
readme <- read_file("README.md")
nonlinear <- read_file("Non-linear Curve Fitting instructions.md")
python_doc <- read_file("PYTHON instructions.md")
rnascope <- read_file("RNAscope instructions.md")

stopifnot(contains(setup, "ANALYSIS_ROOT"))
stopifnot(contains(setup, "path.expand('~')"))
stopifnot(!contains(setup, "Sys.getenv('USER')"))
stopifnot(!contains(setup, "file.path('/Users'"))

stopifnot(!contains(readme, "`images/`"))
stopifnot(contains(readme, "ANALYSIS_ROOT"))
stopifnot(contains(readme, "Rtools must be installed and available on `PATH`"))
stopifnot(!contains(readme, "shinybusy"))

stopifnot(contains(python_doc, "pip install pillow"))
stopifnot(contains(python_doc, "pip install tzdata"))

stopifnot(!contains(nonlinear, "writexl"))
stopifnot(!contains(nonlinear, "write_xlsx"))
stopifnot(contains(nonlinear, "openxlsx::write.xlsx"))
stopifnot(contains(nonlinear, "readxl"))
stopifnot(contains(nonlinear, "Rtools must be installed and available on `PATH`"))
stopifnot(!contains(nonlinear, "Sys.getenv('USER')"))
stopifnot(!contains(nonlinear, "file.path('/Users'"))
stopifnot(!contains(nonlinear, "'roma'"))
stopifnot(!contains(nonlinear, "Sheet 1"))

stopifnot(contains(rnascope, "ANALYSIS_ROOT"))
stopifnot(contains(rnascope, "NWBdata/001832/sub-L1-ST8/sub-L1-ST8_ses-20240905T115902.nwb"))

cat("portability/doc checks passed\n")
