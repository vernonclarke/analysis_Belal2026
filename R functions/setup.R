# ==============================================
# SHARED ANALYSIS SETUP
# ==============================================

load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}

required.packages <- c('robustbase', 'minpack.lm', 'openxlsx',
                       'Rcpp', 'signal', 'dbscan')
load_required_packages(required.packages)

if (!exists('root_dir')) {
  root_dir <- file.path(path.expand('~'),
                        'Documents', 'Repositories', 'analysis_Belal2026')
}

repo_root <- Sys.getenv('ANALYSIS_ROOT', unset=root_dir)
repo_root <- normalizePath(repo_root, mustWork=TRUE)

source(file.path(repo_root, 'R functions', 'nNLS functions.R'))

make_paths <- function(identifier) {
  analysis_path <- paste0(repo_root, '/Paper analysis/', identifier)
  xlsx_path <- paste0(analysis_path, '/xlsx')
  svg_path <- paste0(analysis_path, '/svg')

  if (!dir.exists(analysis_path)) stop('Analysis folder not found: ', analysis_path)
  if (!dir.exists(xlsx_path)) dir.create(xlsx_path, recursive=TRUE)
  if (!dir.exists(svg_path)) dir.create(svg_path, recursive=TRUE)

  list(analysis_path=analysis_path,
       xlsx_path=xlsx_path,
       svg_path=svg_path)
}
