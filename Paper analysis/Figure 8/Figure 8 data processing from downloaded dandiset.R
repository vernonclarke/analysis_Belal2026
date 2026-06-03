# ==============================================
# EXTRACT FIGURE 8 SUMMARY DATA FROM DANDI NWB FILES
# Extracting data for fits to XLSX files
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))
load_required_packages(c('jsonlite', 'reticulate', 'readABF', 'yaml'))

env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
}
reticulate::conda_install(env_name, packages = c('pip', 'pynwb', 'numpy'), channel = 'conda-forge')
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

identifier <- 'Figure 8'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

dandi_root <- file.path(repo_root, 'NWBdata', '001832')

if (!dir.exists(dandi_root)) {
  stop('Downloaded DANDI folder not found: ', dandi_root)
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

if (is.null(dataset_description$FigureMappings[['Figure 8']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig8 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 8']], as.data.frame)
)

fig8 <- subset(fig8, grepl('_icephys\\.nwb$', dandi_path))

fig8$nwb_filepath <- file.path(dandi_root, fig8$dandi_path)
fig8$subject <- sub('/.*$', '', fig8$dandi_path)

missing_files <- fig8$nwb_filepath[!file.exists(fig8$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 8 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig8 <- fig8[file.exists(fig8$nwb_filepath), ]
if (nrow(fig8) == 0) {
  stop('No local Figure 8 icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- fig8$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
fig8 <- fig8[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig8$subject)

cat('\nFigure 8 NWB files found:', nrow(fig8), '\n')
print(fig8[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c('Control for MCI-Park', 'ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK')

summary_all <- lapply(seq_len(nrow(fig8)), function(ii) {
  nwb_filepath <- fig8$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# group by genotype metadata
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig8))) {
  nwb_filepath <- fig8$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  genotype <- py_run_string("
from pynwb import NWBHDF5IO

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
genotype = str(nwb.subject.genotype)
io.close()
", convert = TRUE)$genotype

  if (grepl('Ndufs2', genotype)) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  }
}

cat('\nControl group:', length(summary[[expt_id[1]]]), 'subjects\n')
cat('MCI-PARK group:', length(summary[[expt_id[2]]]), 'subjects\n')

# SAVE DATA IN SEPARATE SHEETS
invisible(
  sapply(seq_along(names(summary)), function(ii) {
    list2excel(
      summary[[ii]],
      paste0(names(summary)[[ii]], ' full.xlsx'),
      wd = xlsx_path,
      center_align = TRUE
    )
  })
)

# SAVE AP5+NBQX+CGP55845A - GABAzine TRACES IN SINGLE SPREADSHEET
summary2 <- lapply(seq_along(summary), function(iii) {
  out <- sapply(seq_along(summary[[iii]]), function(ii) {
    matrix(
      summary[[iii]][[ii]][, 'AP5+NBQX+CGP55845A'] - summary[[iii]][[ii]][, 'GABAzine'],
      ncol = 1
    )
  })

  colnames(out) <- names(summary[[iii]])
  out
})

names(summary2) <- expt_id

# Sort alphabetically
summary2 <- summary2[order(names(summary2))]

names(summary2)

summary2[[1]][1:10,]
#          24112003   24112011   24125001 24125007    24125011    24209016
#  [1,]  0.49273679  0.1373779 -0.7027588 1.715963  0.95855708  0.56657712
#  [2,]  0.01462809 -0.7903564  0.7620849 1.858378  0.13458251  1.17692866
#  [3,] -0.39227293  0.3815185 -0.5196533 2.204244  0.40924070 -0.33674315
#  [4,] -1.35866286  1.5289794 -1.4962157 2.102519 -1.26922601 -1.16682123
#  [5,] -1.50107822  1.6266357 -2.8389891 2.407694 -0.81146236 -1.21564936
#  [6,] -0.85003658 -0.1311768 -3.3883055 1.838033 -0.04852295 -0.11701660
#  [7,]  0.56394447 -0.9368408 -0.2755127 1.207336  0.95855708 -0.33674315
#  [8,]  0.60463457  0.5035888 -0.5806884 1.044576 -0.75042721  0.05388183
#  [9,] -0.20916747  1.2115966 -1.7403564 1.736308 -0.17059325  0.61540524
# [10,] -1.12469477  0.5035888  0.1517334 1.003886  0.13458251  0.73747555

summary2[[2]][1:10,]
#         24115002    24115010   24306011   24312006   24312007  24313006    24313010   24313012
#  [1,]  1.1040527 -0.02106933  1.7751769  0.0723470  1.5287109 -4.120044  0.01062012 -2.0571532
#  [2,] -0.2387207 -0.48493650  2.0193175  1.6592610  0.9671875 -2.801684 -1.30163568 -0.6411377
#  [3,] -2.6557128  1.19963373  0.5239563  0.1130371  1.2113281 -1.825122 -1.51525872  1.0922607
#  [4,] -1.8744628  1.34611810 -0.3610534  0.1130371  1.5775390 -2.288989  0.04113769 -0.4214111
#  [5,] -2.1430175 -0.99763179 -0.6357116  0.2757975  1.2357421 -4.242114  0.86511226 -0.8364502
#  [6,] -2.4604003  0.46721189 -1.0934753  0.8454589  2.9203124 -6.048755  2.14685048 -1.0805908
#  [7,] -2.5580565  1.39494622 -0.6357116  2.4730630  1.9437499 -5.365161  2.78771959  0.3354248
#  [8,] -0.4096191  0.44279783  0.3103332  2.3509927 -1.1080078 -4.193286  1.65856926  1.8979247
#  [9,]  1.4946777  0.56486814  0.1577454  0.5606282  0.4300781 -3.143481  0.37683104  0.7992920
# [10,]  2.2515136  1.07756343 -0.4220886 -1.1483561  1.6263671 -3.021411 -0.29455565 -0.6167236

save <- TRUE

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
