# ==============================================
# EXTRACT FIGURE 9 SUMMARY DATA FROM DANDI NWB FILES
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
installed_packages <- reticulate::py_list_packages(envname = env_name, type = 'conda')$package
if (length(setdiff(c('pip', 'pynwb', 'numpy'), installed_packages))) {
  reticulate::conda_install(env_name, packages = c('pip', 'pynwb', 'numpy'), channel = 'conda-forge')
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

identifier <- 'Figure 9'
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

if (is.null(dataset_description$FigureMappings[['Figure 9']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig9 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 9']], as.data.frame)
)

fig9 <- subset(fig9, grepl('_icephys\\.nwb$', dandi_path))

fig9$nwb_filepath <- file.path(dandi_root, fig9$dandi_path)
fig9$subject <- sub('/.*$', '', fig9$dandi_path)

missing_files <- fig9$nwb_filepath[!file.exists(fig9$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 9 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig9 <- fig9[file.exists(fig9$nwb_filepath), ]
if (nrow(fig9) == 0) {
  stop('No local Figure 9 icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- fig9$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
fig9 <- fig9[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig9$subject)

cat('\nFigure 9 NWB files found:', nrow(fig9), '\n')
print(fig9[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c(
  'ChAT-Cre X De eGFP',
  'ChAT-Cre X De eGFP 6OHDA',
  'ChAT-Cre X tdTomato',
  'ChAT-Cre X tdTomato 6OHDA'
)

summary_all <- lapply(seq_len(nrow(fig9)), function(ii) {
  nwb_filepath <- fig9$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# group by genotype and session metadata
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig9))) {
  nwb_filepath <- fig9$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  result <- py_run_string("
from pynwb import NWBHDF5IO

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
genotype = str(nwb.subject.genotype)
session_desc = str(nwb.session_description)
io.close()
", convert = TRUE)

  genotype <- result$genotype
  session_desc <- result$session_desc

  is_iSPN <- grepl('Drd2', genotype)
  is_dSPN <- grepl('Drd1a', genotype)
  is_6ohda <- grepl('6-OHDA', session_desc)

  if (is_iSPN && !is_6ohda) {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_iSPN && is_6ohda) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_dSPN && !is_6ohda) {
    summary[[expt_id[3]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_dSPN && is_6ohda) {
    summary[[expt_id[4]]][[cell_id]] <- summary_all[[ii]]
  }
}

cat('\nChAT-Cre X De eGFP (iSPN control):', length(summary[[expt_id[1]]]), 'subjects\n')
cat('ChAT-Cre X De eGFP 6OHDA (iSPN lesion):', length(summary[[expt_id[2]]]), 'subjects\n')
cat('ChAT-Cre X tdTomato (dSPN control):', length(summary[[expt_id[3]]]), 'subjects\n')
cat('ChAT-Cre X tdTomato 6OHDA (dSPN lesion):', length(summary[[expt_id[4]]]), 'subjects\n')

# Save data in separate sheets
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

# Save AP5+NBQX+CGP55845A traces in single spreadsheet
summary2 <- lapply(seq_along(summary), function(iii) {
  out <- sapply(seq_along(summary[[iii]]), function(ii) {
    matrix(summary[[iii]][[ii]][, 'AP5+NBQX+CGP55845A'], ncol = 1)
  })

  colnames(out) <- names(summary[[iii]])
  out
})

names(summary2) <- expt_id

# Sort alphabetically
summary2 <- summary2[order(names(summary2))]

names(summary2)

summary2[[1]][1:10,]
#         24610000    24610005    24610006    24610007     24611003   24618001
#  [1,]  1.3520914 -0.26652831  1.80830883 -1.33229974  0.922949175 -1.5257975
#  [2,]  2.1658935  0.07526855  0.42484536 -1.60085442  1.460058524 -0.9561360
#  [3,]  0.8638102 -0.19328612 -1.48758945 -1.13698725  0.166113273 -0.5899251
#  [4,]  1.5962320  0.29499510  0.14001464 -0.16042480 -0.004785156 -0.5899251
#  [5,]  1.9624429 -0.09562988  1.23864740  0.03488769  0.629980439 -0.3457845
#  [6,]  0.3755290 -0.68156735 -0.06343587 -0.77077633  0.434667948 -0.3457845
#  [7,] -0.4789632  0.12409667 -0.67378740 -1.96706534  0.996191359 -0.4678548
#  [8,]  0.5789795  0.34382323 -0.02274577 -2.28444813  1.069433543 -0.1423340
#  [9,]  3.0610757  0.78327633  0.62829587 -0.96608882  0.703222623 -0.4271647
# [10,]  3.9969481  0.61237790 -0.02274577 -0.57546384  0.581152316 -1.5257975

summary2[[2]][1:10,]
#         23n13001   23n13007   23n14002    23n16015   23d04007
#  [1,] -1.5887450 -0.5537719 -4.1129148 -0.11442057  0.3717041
#  [2,] -2.5042723 -0.4927368 -2.3428954  1.59456373 -0.2996826
#  [3,] -2.0159911 -0.3706665 -0.8780517  0.86214189  0.7989502
#  [4,] -1.9549560 -0.1875610 -1.9156493 -1.82340486  1.2261962
#  [5,] -1.8939208 -0.5537719 -3.0142821 -1.70133455  2.2027587
#  [6,] -2.4432372  0.6669311 -2.9532469 -1.25374343  1.8365478
#  [7,] -2.3822020  0.7279663 -3.8687742 -1.00960282  1.4093017
#  [8,] -1.8939208  1.2772827 -2.4039305 -0.48063149  0.6768799
#  [9,] -0.5511474  0.9721069 -2.6480712  0.04833984 -0.1165771
# [10,] -1.1004638 -0.2485962 -3.6856688 -1.00960282 -0.2386474

summary2[[3]][1:10,]
#        24610004    24611000    24611006    24613001    24613002    24613003    24710004    24710006    24711000    24711005
#  [1,] 1.3216857 -0.48500974  0.26745604  0.08459472  1.41210931 -0.08306884  0.86849971 -0.65002438 -0.25374348  0.42807615
#  [2,] 1.8404845 -0.33852537 -0.14758300  0.35314940  0.28906249 -0.38824461  0.13607787 -0.24312336 -0.41650389 -0.13344726
#  [3,] 0.2840881 -0.26528319  0.78015133  0.10900878 -0.29687499 -1.30377191 -0.56582639 -0.08036295  0.64143877  0.89194332
#  [4,] 1.3522033  0.36948240  0.48718259  0.52404783 -0.49218748 -0.75445553 -1.38980096 -0.56864418  1.86214184  1.59995110
#  [5,] 0.3756409  0.17416991 -0.88000484 -0.55017087 -0.56542966  0.03900146 -1.51187127 -1.05692541  1.69938143  0.89194332
#  [6,] 0.7418518  0.66245114  0.73132321 -0.98962398  0.50878904  0.25262450 -0.50479124 -0.32450357 -0.05029297  0.03745117
#  [7,] 0.1315002  0.66245114 -0.04992676 -0.74548336  0.92382808  0.52728269 -0.04702759  0.44860838  0.35660806  1.23374018
#  [8,] 0.5587463  0.63803708 -0.02551269  0.52404783 -0.15039062  0.46624754 -0.65737912  0.48929848  0.60074867  1.99057608
#  [9,] 1.4742736 -0.04555664  0.02331543  0.71936032 -0.02832031  0.92401119 -1.05410762  1.22172032  0.31591795  0.30600584
# [10,] 0.9554748  0.12534179 -0.26965331  0.40197752  1.11914057  0.64935300 -2.03067007  1.54724114  0.84488928 -0.62172849

summary2[[4]][1:10,]
#         23n08005   23n08007   23n08010    23n08020   23n09001   23n09013   23n16010    23n22029     23n22031
#  [1,] -0.4714762 -0.8011474  0.1453369  0.06722005  1.5424194  0.1277262 -0.2062988 -0.24798583  0.440592427
#  [2,]  0.3830159  0.3381754 -0.6359131  0.02652995  1.2372436  0.2497965 -1.3049316 -0.55316159  0.542317683
#  [3,]  1.1968180 -0.3128662 -0.3917724 -0.09554036  0.2606811  0.2497965  1.1364746  0.17926025  0.359212223
#  [4,]  0.8306071 -1.3301188  1.0486572 -0.42106118 -0.4107055  0.6160075  1.5026855  0.11822509  0.684733040
#  [5,]  0.5864664 -0.4756266  0.4871338 -1.31624343  2.0307006  2.0808511  1.9909667  0.05718994  0.094726558
#  [6,]  0.5864664 -1.2080484 -0.3673584 -1.47900384  3.2514037  0.1277262 -0.2062988 -1.89593497 -0.291829413
#  [7,]  0.7085368 -1.2080484 -1.3683349  0.14860025  0.5048218 -0.4419352 -1.4270019 -1.28558344  0.176106762
#  [8,]  0.9119873 -1.4928792 -1.0265380  0.92171220  1.2982787 -0.3605550 -0.6335449  1.09478755  0.562662734
#  [9,] -0.7156168 -0.7197672 -0.4161865  0.39274087  1.0541381 -0.8081461  1.2585449  1.52203362 -0.006998698
# [10,] -0.6749267 -0.2721761 -0.2452881 -0.33968097 -0.4107055 -0.6046956  2.5402831 -0.67523190 -0.047688800

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
