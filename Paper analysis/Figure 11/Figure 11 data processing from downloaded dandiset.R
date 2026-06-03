# ==============================================
# EXTRACT FIGURE 11 SUMMARY DATA FROM DANDI NWB FILES
# Extracting data for fits to XLSX files
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

save <- TRUE

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))
load_required_packages(c('jsonlite', 'reticulate', 'readABF', 'yaml'))

env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
  reticulate::conda_install(env_name, packages = c('pynwb', 'numpy'), pip = TRUE)
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

identifier <- 'Figure 11'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

dandi_root <- file.path(repo_root, 'NWBdata', '001832')

if (!dir.exists(dandi_root)) {
  stop('Downloaded DANDI folder not found: ', dandi_root)
}

abf2nwb_functions_path <- file.path(repo_root, 'R functions', 'ABF2NWB_functions.R')
if (!file.exists(abf2nwb_functions_path)) {
  stop('ABF2NWB_functions.R not found at: ', abf2nwb_functions_path)
}
source(abf2nwb_functions_path)

if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

dataset_description_path <- file.path(dandi_root, 'dataset_description', 'dataset_description.json')
if (!file.exists(dataset_description_path)) {
  stop('dataset_description.json not found at: ', dataset_description_path)
}

dataset_description <- jsonlite::fromJSON(dataset_description_path, simplifyVector = FALSE)

if (is.null(dataset_description$FigureMappings[['Figure 11']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig11 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 11']], as.data.frame)
)

fig11 <- subset(fig11, grepl('_icephys\\.nwb$', dandi_path))

fig11$nwb_filepath <- file.path(dandi_root, fig11$dandi_path)
fig11$subject <- sub('/.*$', '', fig11$dandi_path)

missing_files <- fig11$nwb_filepath[!file.exists(fig11$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 11 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig11 <- fig11[file.exists(fig11$nwb_filepath), ]
if (nrow(fig11) == 0) {
  stop('No local Figure 11 icephys NWB files found under: ', dandi_root)
}

sort_key <- fig11$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)

fig11 <- fig11[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig11$subject)

cat('\nFigure 11 NWB files found:', nrow(fig11), '\n')
print(fig11[, c('subject', 'dandi_path')])

expt_id <- c(
  'ChI-NGF control',
  'ChI-NGF 6OHDA'
)

summary_all <- lapply(seq_len(nrow(fig11)), function(ii) {
  nwb_filepath <- fig11$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig11))) {
  nwb_filepath <- fig11$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  result <- py_run_string("
from pynwb import NWBHDF5IO

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
session_desc = str(nwb.session_description)
io.close()
", convert = TRUE)

  session_desc <- result$session_desc
  is_6ohda <- grepl('6-OHDA', session_desc, ignore.case = TRUE)

  if (is_6ohda) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  }
}

cat('\nChI-NGF control:', length(summary[[expt_id[1]]]), 'subjects\n')
cat('ChI-NGF 6OHDA:', length(summary[[expt_id[2]]]), 'subjects\n')


summary2 <- lapply(seq_along(summary), function(iii) {
  out <- sapply(seq_along(summary[[iii]]), function(ii) {
    matrix(summary[[iii]][[ii]][, 'control'], ncol = 1)
  })

  colnames(out) <- names(summary[[iii]])
  out
})

names(summary2) <- expt_id

summary2 <- summary2[order(names(summary2))]

names(summary2)

summary2[[1]][1:10,]
#         25310000    25310001   25310005    25310006    25310007   25310012   25311000    25311002  25311007    25312000     25312010
#  [1,] -1.9458739 -0.79577633 -0.3421936 -1.59694817 -0.14765624  0.3149414  2.2543212 -1.33498529 -3.186767  0.90334468  1.496166921
#  [2,] -2.2388427 -0.42956541  0.3291931 -2.01198721  0.46269529  0.6323242  0.2767822 -0.26076659 -2.271240  0.02443848  1.129956001
#  [3,] -1.3599365 -0.73474118 -0.7694397 -1.40163568  0.36503905  0.2172851 -0.4800537 -0.13869628 -2.393310 -0.24411620  0.226635731
#  [4,] -1.5796630 -0.49060056  0.2376404 -0.30300292  0.07207031 -0.1245117 -0.2603271  0.03220215 -3.125732  0.19533690 -0.261645495
#  [5,] -1.2622802  0.18078612  0.5733337 -0.03444824  0.56035154 -0.4663086  0.5941650 -0.48049314 -2.973144  0.85451656 -0.261645495
#  [6,] -1.0913818 -0.55163572  0.1766052 -0.25417479  0.48710935 -0.4907226  0.4720947  0.12985839 -1.966064  0.75686032  0.006909179
#  [7,] -1.0425537  0.27233885  1.2447204 -1.01101069 -0.17207030  0.2661133 -1.1392333  1.03317866 -2.362793  0.68361813 -0.066333005
#  [8,] -0.1392334 -0.39904783  0.6648864 -1.59694817  0.21855468  0.6323242  0.2035400 -0.04104004 -1.752441  0.31740721  0.788159142
#  [9,] -0.4810303 -0.49060056  0.3597107 -0.62038571  0.46269529  0.5590820  0.1547119 -0.38283690 -1.538818 -0.48825681  2.130932516
# [10,] -1.2622802  0.05871582 -0.8915100  0.16086425  0.12089843 -0.1977539 -0.3335693 -0.13869628 -1.660889 -0.87888179  1.984448148

summary2[[2]][1:10,]
#          24d12009   24d12013   25304011    25306003    25306008    25409011   25409015
#  [1,] -0.48649900  0.6825439  1.1930847 -1.39470208 -0.01547851 1.034252881 -0.7939941
#  [2,] -0.36442869 -0.5381592  0.2165222 -0.56462400 -0.06430664 0.008862304 -0.2080566
#  [3,]  0.17268066 -0.2696045 -0.2717590 -0.49138181 -0.84555660 0.008862304  0.4999511
#  [4,] -0.53532712 -0.8311279  0.2470398 -0.05192871 -0.45493162 0.423901347  0.3290527
#  [5,] -0.55974118  0.1210205  0.9794616 -0.63786618 -0.77231442 0.423901347 -0.3057129
#  [6,] -0.51091306 -0.3916748  1.4677429 -0.46696775 -0.57700193 0.448315408 -1.6729003
#  [7,]  0.07502441 -0.2207764  1.3151550  0.21662597 -0.21079101 0.375073224 -1.1113769
#  [8,]  0.49006345  0.1942627  0.9794616  0.33869627 -0.60141599 0.204174795 -1.3066894
#  [9,]  0.09943847 -0.5381592  0.1860046  0.77814938 -0.67465817 0.985424758  0.3046387
# [10,] -0.70622555 -1.1240966  0.6437683  1.58381340 -0.50375974 1.351635678  0.4267090

if (save) {
  for (iii in seq_along(summary2)) {
    out <- summary2[[iii]]

    wb <- createWorkbook()
    addWorksheet(wb, 'Sheet1')
    writeData(wb, sheet = 'Sheet1', x = out, startRow = 1, startCol = 1, rowNames = FALSE)

    center_style <- createStyle(halign = 'center', valign = 'center')
    addStyle(
      wb,
      sheet = 'Sheet1',
      style = center_style,
      rows = 1:(nrow(out) + 1),
      cols = 1:ncol(out),
      gridExpand = TRUE
    )

    saveWorkbook(
      wb,
      file = file.path(xlsx_path, paste0(names(summary2)[iii], '.xlsx')),
      overwrite = TRUE
    )
  }
}
