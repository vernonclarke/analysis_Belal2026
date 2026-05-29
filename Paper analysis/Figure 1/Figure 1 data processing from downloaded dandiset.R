##################################### EXTRACT FIGURE 1 SUMMARY DATA FROM DANDI NWB FILES #####################################

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

identifier <- 'Figure 1'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)
xlsx_path <- file.path(analysis_path, 'xlsx')
figure2_xlsx_path <- file.path(repo_root, 'Paper analysis', 'Figure 2', 'xlsx')
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
if (is.null(dataset_description$FigureMappings[['Figure 1']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig1 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 1']], as.data.frame)
)
fig1 <- subset(fig1, grepl('_icephys\\.nwb$', dandi_path))
fig1$nwb_filepath <- file.path(dandi_root, fig1$dandi_path)
fig1$subject <- sub('/.*$', '', fig1$dandi_path)

missing_files <- fig1$nwb_filepath[!file.exists(fig1$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 1 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig1 <- fig1[file.exists(fig1$nwb_filepath), ]
if (nrow(fig1) == 0) {
  stop('No local Figure 1 icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- fig1$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
fig1 <- fig1[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig1$subject)

cat('\nFigure 1 NWB files found:', nrow(fig1), '\n')
print(fig1[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c('ChAT-Cre X De eGFP', 'ChAT-Cre X tdTomato')

summary_all <- lapply(seq_len(nrow(fig1)), function(ii) {
  nwb_filepath <- fig1$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# group by genotype
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig1))) {
  nwb_filepath <- fig1$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  genotype <- py_run_string('
from pynwb import NWBHDF5IO
io = NWBHDF5IO(nwb_filepath, "r")
nwb = io.read()
genotype = str(nwb.subject.genotype)
io.close()
', convert = TRUE)$genotype

  if (grepl('Drd2-eGFP', genotype)) {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  } else if (grepl('Drd1a-tdTomato', genotype)) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else {
    warning('Unrecognized genotype for ', cell_id, ': ', genotype)
  }
}

cat('\niSPN group:', length(summary[[expt_id[1]]]), 'subjects\n')
cat('dSPN group:', length(summary[[expt_id[2]]]), 'subjects\n')

summary <- summary[order(names(summary))]

# GABAzine traces
summary2 <- lapply(seq_along(summary), function(iii) {
  out <- lapply(seq_along(summary[[iii]]), function(ii) {
    matrix(
      summary[[iii]][[ii]][, 'control'] - summary[[iii]][[ii]][, 'GABAzine'],
      ncol = 1
    )
  })
  names(out) <- names(summary[[iii]])
  out
})
names(summary2) <- names(summary)

out <- lapply(seq_along(names(summary2)), function(ii) {
  valid_matrices <- Filter(function(x) !is.null(x) && is.matrix(x), summary2[[ii]])

  max_rows <- max(sapply(valid_matrices, nrow))

  padded <- lapply(valid_matrices, function(m) {
    if (nrow(m) < max_rows) {
      rbind(m, matrix(NA, max_rows - nrow(m), ncol(m)))
    } else {
      m
    }
  })

  out1 <- do.call(cbind, padded)
  colnames(out1) <- names(valid_matrices)
  out1
})

names(out) <- names(summary2)
names(summary2)

out[[1]][1:10, ]
#        22316001   22322006   22324018  22329005   24o21000   24o22005
#  [1,] -2.667175 -3.4880980 -0.4484456 -2.482568 -1.1648559  0.6208191
#  [2,] -7.122741 -3.4880980  1.5046793 -3.337060 -4.2166136  1.5973815
#  [3,] -3.948913 -2.5725707  2.5626219 -3.190576 -6.9631955 -0.9660949
#  [4,]  3.009094 -1.7180785  2.2777912 -3.752099 -4.5217893 -2.5835265
#  [5,] -1.385437 -2.5115355  1.9929605 -2.629053 -1.1648559 -0.2947082
#  [6,] -8.526550 -0.1311645  2.3998615 -1.896631 -0.5545044  1.5058288
#  [7,]  2.581848 -0.7415161  3.7019448 -1.676904 -2.0803832  0.2240906
#  [8,] -3.521667 -0.9856567  4.0274656 -2.751123 -1.1648559 -0.2641907
#  [9,] -5.902038 -1.4129028  2.9695230 -3.068506 -1.1648559  0.8954773
# [10,]  2.459778 -1.5960082  2.4405516 -3.630029 -4.2166136 -0.6914367

out[[2]][1:10, ]
#       22405005    22405015   22406005   22406013  22411000   22411004  22412005    24o21006    24o21010   24o22003
#  [1,] 5.641072  0.94885249  1.5003356 -1.8062743 -4.265137  1.9136352  8.364705 -0.40327960 -1.55700676 -0.5395914
#  [2,] 7.472127 -2.26566558  1.9886168 -1.9283446 -3.135986  1.6084594  7.632283 -0.52534991 -0.78389482 -2.1265054
#  [3,] 9.099731 -1.24841303  1.7749938 -1.7248941 -4.661865  0.6318969  8.608846 -0.52534991  0.23335774 -2.2078856
#  [4,] 9.303182  0.54195147  0.7679138 -1.1145426 -3.074951 -0.8939819  8.039184 -1.66467277  1.16923009 -3.3065184
#  [5,] 7.879028  0.94885249  2.1106872 -1.2366129 -3.013916 -1.4432983  7.795044 -0.89156083 -0.49906410 -3.4285887
#  [6,] 7.268676 -0.71944170  1.9886168 -0.9517822 -3.532715  0.2656860  9.951619 -0.68811032 -1.55700676 -2.1671955
#  [7,] 8.285929  0.05367024  0.7679138 -0.7076416 -2.525635 -0.7108764 10.114379 -0.03706868 -1.15010574 -2.3706460
#  [8,] 7.879028 -1.61462395 -1.2157287 -0.4228109 -3.105469  0.1436157 10.114379  0.41052244 -0.53975421 -0.3361409
#  [9,] 7.268676 -1.24841303 -0.6053772  0.1061605 -4.295654  0.5098266  9.666788  0.16638183 -0.09216308 -1.2313232
# [10,] 6.251424 -2.67256660  0.4017029  0.8385823 -4.143066 -0.2836304  8.324015 -0.32189940 -0.37699380 -2.6147867


# save full extracted data to XLSX
save <- TRUE
invisible(
  sapply(seq_along(names(summary)), function(ii) {
    list2excel(
      summary[[ii]],
      paste0(names(summary)[[ii]], '.xlsx'),
      wd = xlsx_path,
      center_align = TRUE
    )
  })
)

# save control - GABAzine traces for Figure 2
if (save) {
  if (!dir.exists(figure2_xlsx_path)) {
    dir.create(figure2_xlsx_path, recursive = TRUE)
  }

  invisible(
    sapply(seq_along(names(summary2)), function(ii) {
      list2excel(
        summary2[[ii]],
        paste0(names(summary2)[[ii]], '.xlsx'),
        wd = figure2_xlsx_path,
        center_align = TRUE
      )
    })
  )
}
