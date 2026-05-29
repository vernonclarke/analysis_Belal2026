# ==============================================
# EXTRACT FIGURE 3B SUMMARY DATA FROM DANDI NWB FILES
# Extracting data for fits to XLSX files
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}

required.packages <- c(
  'dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp',
  'reticulate', 'robustbase', 'signal', 'readABF', 'yaml'
)
load_required_packages(required.packages)

env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
  reticulate::conda_install(env_name, packages = c('pynwb', 'numpy'), pip = TRUE)
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

# paths
repo_root <- normalizePath(
  '~/Documents/Repositories/analysis_Belal2026',
  mustWork = TRUE
)

identifier <- 'Figure 3B'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)
xlsx_path <- file.path(analysis_path, 'xlsx')
dandi_root <- file.path(repo_root, 'NWBdata', '001832')

if (!dir.exists(dandi_root)) {
  stop('Downloaded DANDI folder not found: ', dandi_root)
}

if (!dir.exists(analysis_path)) {
  stop('Analysis folder not found at: ', analysis_path)
}

# source helper functions
abf2nwb_functions_path <- file.path(repo_root, 'R functions', 'ABF2NWB_functions.R')
if (!file.exists(abf2nwb_functions_path)) {
  stop('ABF2NWB_functions.R not found at: ', abf2nwb_functions_path)
}
source(abf2nwb_functions_path)

# Output folder
if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

# locate NWB files from mapping in dataset_description_path
dataset_description_path <- file.path(dandi_root, 'dataset_description', 'dataset_description.json')
if (!file.exists(dataset_description_path)) {
  stop('dataset_description.json not found at: ', dataset_description_path)
}

dataset_description <- jsonlite::fromJSON(dataset_description_path, simplifyVector = FALSE)

if (is.null(dataset_description$FigureMappings[['Figure 3']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig3b <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 3']], as.data.frame)
)

# FigureMappings[['Figure 3']] contains Figure 3 and Figure 3B, so keep only Figure 3B.
fig3b <- subset(
  fig3b,
  grepl('^Figure 3B/', original_path) & grepl('_icephys\\.nwb$', dandi_path)
)

fig3b$nwb_filepath <- file.path(dandi_root, fig3b$dandi_path)
fig3b$subject <- sub('/.*$', '', fig3b$dandi_path)

missing_files <- fig3b$nwb_filepath[!file.exists(fig3b$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 3B NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig3b <- fig3b[file.exists(fig3b$nwb_filepath), ]
if (nrow(fig3b) == 0) {
  stop('No local Figure 3B icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- fig3b$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
fig3b <- fig3b[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig3b$subject)

cat('\nFigure 3B NWB files found:', nrow(fig3b), '\n')
print(fig3b[, c('subject', 'dandi_path')])

# extract averages
summary <- lapply(seq_len(nrow(fig3b)), function(ii) {
  nwb_filepath <- fig3b$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary) <- cell_ids

expt_name <- 'NDNF GABA PSCs'

if (save) {
  list2excel(
    summary,
    paste0(expt_name, '_full.xlsx'),
    wd = xlsx_path,
    center_align = TRUE
  )
}

# SAVE CONTROL - GABAZINE DIFFERENCE
summary2 <- sapply(summary, function(cell_data) {
  if ('control' %in% colnames(cell_data) && 'GABAzine' %in% colnames(cell_data)) {
    matrix(cell_data[, 'control'] - cell_data[, 'GABAzine'], ncol = 1)
  } else {
    NULL
  }
})

colnames(summary2) <- cell_ids

summary2[1:10, ]
#          24d02000 24d02006  24d03000   24d03002   24d05002    24d05007    26211001    26211007
#  [1,]  2.55715930 1.481018 -2.659241 -1.9159667  1.5502013  0.12875366  0.10556640  0.99190262
#  [2,]  1.51956170 1.786194 -1.804748 -0.9882324 -0.7080993 -0.23745726 -0.35830076  1.39880364
#  [3,] -0.18942260 1.877746 -2.079407 -2.1601073 -1.3794860  2.44808948  0.61826169  1.11397293
#  [4,] -0.03683472 1.328430 -2.354065 -1.0614746 -0.6775818 -1.64126579  1.05771479  0.21879068
#  [5,]  0.20730590 2.152405 -3.269592  0.1592285 -0.8301696 -2.89248643 -0.30947264 -0.39156085
#  [6,] -0.64718625 1.206360 -3.391662 -0.9394043  1.8248595 -2.46524036 -0.40712889 -0.26949055
#  [7,] -1.28805536 2.854309 -3.117004 -1.8671386  1.3365783 -3.07559189 -0.08974609  0.05603027
#  [8,] -0.86080929 3.006897 -2.476135 -1.3056152  1.8553771 -0.75625607 -0.28505858  0.95121252
#  [9,] -0.28097533 2.549133 -1.163879 -0.9394043  3.1371153  0.52548215  0.44736326  0.46293129
# [10,]  0.20730590 1.236877 -1.621643 -1.4765136  2.9234923  0.09823608  0.88681636 -0.71708167

save <- TRUE
if (save) {
  wb <- createWorkbook()
  addWorksheet(wb, 'Sheet1')
  writeData(wb, sheet = 'Sheet1', x = summary2, startRow = 1, startCol = 1, rowNames = FALSE)

  center_style <- createStyle(halign = 'center', valign = 'center')
  addStyle(
    wb,
    sheet = 'Sheet1',
    style = center_style,
    rows = 1:(nrow(summary2) + 1),
    cols = 1:ncol(summary2),
    gridExpand = TRUE
  )

  saveWorkbook(
    wb,
    file = file.path(xlsx_path, paste0(expt_name, '.xlsx')),
    overwrite = TRUE
  )
}