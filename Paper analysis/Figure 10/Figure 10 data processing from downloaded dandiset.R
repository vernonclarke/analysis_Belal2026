# ==============================================
# EXTRACT FIGURE 10 SUMMARY DATA FROM DANDI NWB FILES
# Extracting data for fits to XLSX files
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')
load_required_packages(c('jsonlite', 'reticulate', 'readABF', 'yaml'))

env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
  reticulate::conda_install(env_name, packages = c('pynwb', 'numpy'), pip = TRUE)
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

identifier <- 'Figure 10'
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

if (is.null(dataset_description$FigureMappings[['Figure 10']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig10 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 10']], as.data.frame)
)

fig10 <- subset(fig10, grepl('_icephys\\.nwb$', dandi_path))

fig10$nwb_filepath <- file.path(dandi_root, fig10$dandi_path)
fig10$subject <- sub('/.*$', '', fig10$dandi_path)

missing_files <- fig10$nwb_filepath[!file.exists(fig10$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 10 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig10 <- fig10[file.exists(fig10$nwb_filepath), ]
if (nrow(fig10) == 0) {
  stop('No local Figure 10 icephys NWB files found under: ', dandi_root)
}

# Sort by actual date embedded in subject id
# Format: sub-mYY[m]DD### where [m] is 1-9 or o/n/d for Oct/Nov/Dec
sort_key <- sapply(fig10$subject, function(f) {
  f_stripped <- gsub('^sub-m', '', f)
  yr <- substr(f_stripped, 1, 2)
  mo_char <- substr(f_stripped, 3, 3)
  mo <- switch(
    mo_char,
    'o' = '10',
    'n' = '11',
    'd' = '12',
    sprintf('%02d', as.integer(mo_char))
  )
  day <- substr(f_stripped, 4, 5)
  rest <- substr(f_stripped, 6, nchar(f_stripped))
  paste0(yr, mo, day, rest)
})

fig10 <- fig10[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig10$subject)

cat('\nFigure 10 NWB files found:', nrow(fig10), '\n')
print(fig10[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c(
  'NPY Cre X dSPN WT',
  'NPY Cre X dSPN 6OHDA',
  'NPY Cre X iSPN WT',
  'NPY Cre X iSPN 6OHDA'
)

summary_all <- lapply(seq_len(nrow(fig10)), function(ii) {
  nwb_filepath <- fig10$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# group by genotype and session metadata
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig10))) {
  nwb_filepath <- fig10$nwb_filepath[ii]
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

  if (is_dSPN && !is_6ohda) {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_dSPN && is_6ohda) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_iSPN && !is_6ohda) {
    summary[[expt_id[3]]][[cell_id]] <- summary_all[[ii]]
  } else if (is_iSPN && is_6ohda) {
    summary[[expt_id[4]]][[cell_id]] <- summary_all[[ii]]
  }
}

cat('\nNPY Cre X dSPN WT (dSPN control):', length(summary[[expt_id[1]]]), 'subjects\n')
cat('NPY Cre X dSPN 6OHDA (dSPN lesion):', length(summary[[expt_id[2]]]), 'subjects\n')
cat('NPY Cre X iSPN WT (iSPN control):', length(summary[[expt_id[3]]]), 'subjects\n')
cat('NPY Cre X iSPN 6OHDA (iSPN lesion):', length(summary[[expt_id[4]]]), 'subjects\n')

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

# SAVE CONTROL TRACES IN SINGLE SPREADSHEET
summary2 <- lapply(seq_along(summary), function(iii) {
  out <- sapply(seq_along(summary[[iii]]), function(ii) {
    matrix(summary[[iii]][[ii]][, 'control'], ncol = 1)
  })

  colnames(out) <- names(summary[[iii]])
  out
})

names(summary2) <- expt_id

# Sort alphabetically
summary2 <- summary2[order(names(summary2))]

names(summary2)

summary2[[1]][1:10,]
#          23406002    23406004   23406005   23406007   23406008   23406010   23o30006 23o31021
#  [1,]  0.29165648 -1.32666009 -0.1302490 -0.1090698 -0.6156158 -0.1953430 -0.5953613 2.859090
#  [2,]  0.74942013 -1.40804030 -0.7100830 -0.3532104 -1.2564849  0.1098328 -0.5221191 2.859090
#  [3,]  0.41372679 -1.02148433 -0.2828369  1.0505981 -1.0886383  0.7201843 -0.8639160 3.265991
#  [4,] -0.07455444 -0.22802733 -0.6795654  0.8064575 -0.9970855  0.4760437 -0.7662597 2.859090
#  [5,]  0.16958617 -0.06526692 -0.6795654 -0.6787313 -1.0123443 -0.6225891  0.4056152 3.469442
#  [6,]  0.29165648 -0.87906897 -1.1373290 -0.6380411 -0.5850982 -1.8127746  0.6497558 2.899780
#  [7,] -0.04403686 -1.08251948 -1.3509521  0.8268025 -0.3409576 -2.0569152  0.9427246 2.818400
#  [8,] -0.22714232 -0.65527341 -1.0152587  0.9895629 -0.9360504 -1.1719055  1.0159667 2.859090
#  [9,] -0.28817748 -0.12630208 -0.5880127 -0.2311401 -1.2717437 -0.2563782 -0.4244629 2.981160
# [10,] -0.47128294 -0.77734371 -0.1912842  0.2367960 -1.3022613 -1.2634582 -0.5465332 3.143921

summary2[[2]][1:10,]
#          23315002    23315010    23315012   23315013   23316006    23316014   23330003   23330008    23d18010      23d18011   23d18013    23d20005   23d20008
#  [1,] -0.79193111  1.08734126 -1.36856073  2.7578734 -0.8133951  0.24200438  0.7878418 -1.5223632  0.32849120 -0.8242797460 -0.6147868 -0.44318846 -0.4968261
#  [2,] -1.50400790 -0.05706787 -1.13967890  1.3845825 -1.2202962  0.42510984  0.3605957 -1.3514648  1.24401850 -0.4885864026  0.1583252 -0.17463378 -1.6442870
#  [3,] -1.09710688 -0.05706787 -0.68191525 -1.0568237 -1.2202962 -0.00213623 -0.1276855 -1.2782226  0.02331543 -0.3054809425  0.6872965 -0.78498531 -2.0837401
#  [4,] -0.89365637  0.47698972 -0.14785766 -0.9042358 -1.7289224 -0.24627684 -0.4023437 -0.7655273 -0.40393064  0.0607299776  0.9721272  0.04509277  0.3820801
#  [5,]  0.22532144  1.08734126  1.07284541 -1.5145873 -2.0340982 -0.79559322 -0.1887207 -0.3871094  0.87780758 -0.0003051758  0.2803955 -0.07697754  1.4562988
#  [6,] -0.07985433  0.62957761  0.08102417 -2.1249389 -1.8306477 -0.97869868 -0.1124268 -0.7899414  2.03747549 -0.8242797460 -0.3706461 -0.32111815  1.1633300
#  [7,]  0.02187093  1.62139885 -0.14785766 -2.1249389 -0.9151204 -0.73455807  0.0859375 -1.4369140 -0.46496580 -1.1294555128  0.4024658  1.02165522  0.6018066
#  [8,] -0.07985433  2.76580797  0.38619993 -0.2938843 -1.8306477  0.54718015 -0.1887207 -1.8763671 -0.40393064 -1.5872191629 -0.2485758  1.04606929  0.3332519
#  [9,] -0.07985433  2.07916250  0.30990599  0.6216430 -1.6271972  1.03546138 -0.2955322 -1.4857421 -0.03771972 -0.3359985192 -0.8182373 -0.17463378 -0.5456543
# [10,] -0.07985433  1.16363520  0.76766964  1.3845825 -1.2202962  0.66925046  0.1164551 -1.2293945 -0.28186034  0.5795287811 -0.1265055  0.04509277 -0.4235840

summary2[[3]][1:10,]
#          23412000   23412004    23412005   23412007   23412010    23413000    23413001  23413005   23n02013
#  [1,]  0.40589598  0.4155477  0.40637205  1.0885772  1.8497314  1.38387038 -0.42947386 1.2524617 -1.6004638
#  [2,] -0.03355713  0.4969279  0.04016113  0.5392608  1.2393798  0.79386389  0.30294798 1.3541869 -3.5047606
#  [3,] -0.36314696  0.9038289 -0.32604979  0.4782257  0.7816162  0.06144205  0.57760617 1.6593627 -3.7489012
#  [4,]  0.18616942  1.0869344  0.38602700  0.5697784  0.7205810 -0.65063473  0.54708860 1.5372924 -1.4295654
#  [5,]  0.57679441  1.5141804  0.12154134  0.5697784  1.1173095 -0.91512040  0.18087768 1.2931518 -0.9168701
#  [6,]  0.69886471  1.5141804 -0.40742999  1.2259063  0.7510986 -0.67097979 -0.42947386 1.7814330 -2.7235106
#  [7,]  0.57679441  1.3921101 -0.38708494  0.9359893 -0.2864990  0.91593420 -0.33792113 0.7031453 -3.4071043
#  [8,]  0.40589598  0.8224487 -0.69226071  0.3561554  0.2017822  2.58422839 -0.58206174 0.2352091 -2.1375731
#  [9,]  0.49134519 -0.2354940 -0.30570474 -0.1779022  0.3848877  1.91284171  0.05880737 0.6421102 -2.1375731
# [10,]  0.58900144 -1.3544718 -0.36673989 -0.2847137 -0.2254639  1.89249665 -0.36843870 0.6014201 -2.3084716
 
summary2[[4]][1:10,]
#          23331001    23331004   23331006   23d20000   23d20003   23d20010
#  [1,]  0.61978146 -0.05711364  1.1205749 -0.2655640  0.6435791 -0.7521565
#  [2,]  0.54348752 -0.24021910  0.8459167  0.1616821  1.4248290 -0.3045654
#  [3,]  0.23831176  0.11073303 -0.4358215 -0.4486694  0.7168213  1.6078694
#  [4,] -0.08212280 -0.37754820  0.1440124 -0.2350464 -0.2841553  0.3057861
#  [5,] -0.34152220 -0.33177183  0.4491882  0.4363403 -0.1132568  0.4278564
#  [6,] -0.28048705 -0.10289001 -0.5884094 -0.1434936  0.1552978  1.2823486
#  [7,] -0.05160522 -0.27073668 -1.2597961 -0.2655640  0.4238525  1.0382080
#  [8,]  0.17727660 -0.31651305 -0.7715149  0.1311645  0.4726806 -0.5080159
#  [9,]  0.10098266 -0.53013608 -0.3442688  0.9246215 -0.3329834 -0.5487060
# [10,] -0.17367553 -0.22496032  0.1440124  0.1006470 -0.2109131 -0.4266357

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
