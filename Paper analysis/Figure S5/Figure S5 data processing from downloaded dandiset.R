
# ==============================================
# EXTRACT FIGURE S5 SUMMARY DATA FROM DANDI NWB FILES
# Loads NWB files from downloaded DANDI dataset, sorts ctrl vs 6OHDA,
# saves averaged traces as XLSX
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))
load_required_packages(c('jsonlite', 'reticulate', 'readABF', 'yaml'))

env_name <- 'NWBenv'
if (!env_name %in% reticulate::conda_list()$name) {
  reticulate::conda_create(env_name, python_version = '3.11')
  reticulate::conda_install(env_name, packages = c('pynwb', 'numpy'), pip = TRUE)
}
reticulate::use_condaenv(env_name, required = TRUE)
reticulate::py_config()

identifier <- 'Figure S5'
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

if (is.null(dataset_description$FigureMappings[['Figure S5']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

figS5 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure S5']], as.data.frame)
)

figS5 <- subset(figS5, grepl('_icephys\\.nwb$', dandi_path))

figS5$nwb_filepath <- file.path(dandi_root, figS5$dandi_path)
figS5$subject <- sub('/.*$', '', figS5$dandi_path)

missing_files <- figS5$nwb_filepath[!file.exists(figS5$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure S5 NWB files:\n', paste(missing_files, collapse = '\n'))
}

figS5 <- figS5[file.exists(figS5$nwb_filepath), ]
if (nrow(figS5) == 0) {
  stop('No local Figure S5 icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- figS5$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
figS5 <- figS5[order(sort_key), ]

# Match original AXG-style column names:
# sub-megcell030723-001 -> egcell030723_001
# sub-mskcell010423-001 -> skcell010423_001
# sub-mskd221230-005    -> skd221230_005
cell_ids <- figS5$subject
cell_ids <- gsub('^sub-m', '', cell_ids)
cell_ids <- gsub('-', '_', cell_ids)

cat('\nFigure S5 NWB files found:', nrow(figS5), '\n')
print(figS5[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c('ctrl', '6OHDA')

summary_all <- lapply(seq_len(nrow(figS5)), function(ii) {
  nwb_filepath <- figS5$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# classify ctrl vs 6OHDA from session_description
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(figS5))) {
  nwb_filepath <- figS5$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  result <- py_run_string("
from pynwb import NWBHDF5IO

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
session_desc = str(nwb.session_description)
io.close()
", convert = TRUE)

  is_ctrl <- grepl('control hemisphere', result$session_desc, ignore.case = TRUE)
  grp <- if (is_ctrl) 'ctrl' else '6OHDA'

  summary[[grp]][[cell_id]] <- summary_all[[ii]]
}

cat('\nctrl  :', length(summary[['ctrl']]), 'subjects\n')
cat('6OHDA :', length(summary[['6OHDA']]), 'subjects\n')

# save per-group full xlsx
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

# extract first trace column from each cell
summary2 <- lapply(seq_along(summary), function(iii) {
  folders <- summary[[iii]]

  if (is.null(folders) || length(folders) == 0) {
    return(NULL)
  }

  vals <- lapply(folders, function(dat) {
    dat[, 1]
  })

  n_pts <- sapply(vals, length)

  if (length(unique(n_pts)) > 1) {
    warning(names(summary)[iii], ': unequal lengths -- truncating to shortest')
    min_pts <- min(n_pts)
    vals <- lapply(vals, `[`, seq_len(min_pts))
  }

  out <- do.call(cbind, vals)
  colnames(out) <- names(folders)
  out
})

names(summary2) <- names(summary)

ctrl_order <- c(
  'skcell010423_001', 'skcell010423_006', 'skcell010423_007',
  'egcell030723_001', 'egcell030723_009', 'egcell030923_006', 'egcell030923_007',
  'egcell041223_003', 'egcell041823_002', 'egcell050223_006', 'egcell122822_002',
  'egcell123022_003', 'egcell123022_004', 'egcell123022_007', 'egcell123022_011',
  'egcell123022_013', 'egcell123022_014', 'skd221230_005', 'skd221230_008',
  'skcell122822_010', 'skcell122822_013', 'skcell122822_004'
)

ohda_order <- c(
  'skcell010423_000', 'skcell010423_009', 'egcell030723_003',
  'egcell030723_008', 'egcell030723_013', 'egcell030923_008', 'egcell041223_004',
  'egcell041823_003', 'egcell041823_005', 'egcell050223_007', 'egcell122822_001',
  'egcell122822_004', 'egcell123022_001', 'egcell123022_006', 'egcell123022_015',
  'skd221230_010', 'skd221230_004', 'skd221230_006', 'skd221230_007',
  'skcell122822_000', 'skcell122822_001', 'skcell122822_006', 'skcell122822_009'
)

missing_ctrl <- setdiff(ctrl_order, colnames(summary2[['ctrl']]))
missing_ohda <- setdiff(ohda_order, colnames(summary2[['6OHDA']]))

if (length(missing_ctrl) > 0) {
  stop('Missing ctrl columns: ', paste(missing_ctrl, collapse = ', '))
}

if (length(missing_ohda) > 0) {
  stop('Missing 6OHDA columns: ', paste(missing_ohda, collapse = ', '))
}

summary2[['ctrl']] <- summary2[['ctrl']][, ctrl_order, drop = FALSE]
summary2[['6OHDA']] <- summary2[['6OHDA']][, ohda_order, drop = FALSE]

print(names(summary2))

old_width <- getOption('width')
options(width = 500)

tmp <- summary2[['ctrl']][1:10, , drop = FALSE]
colnames(tmp) <- abbreviate(colnames(tmp), minlength = 8)
print(tmp)

#      s010423_001  s010423_006 s010423_007 e030723_001 e030723_009 e030923_006 e030923_007 e041223_  e041823_ e050223_  e122822_ e123022_003 e123022_004 e123022_007 e123022_011 e123022_013 e123022_014 s221230_005 s221230_008 s122822_010 s122822_013 s122822_00
#  [1,]   0.4693255  0.003964844   1.4076380   1.3584663  2.73106514  -2.5525428  0.66753983  -3.0579  1.055843 1.863288 -3.339204  1.31708626  0.82329556  -1.8999766 -0.33418732  -2.3831674   -2.851821   -3.147856   1.4941929  -1.0806274    2.582153  1.5278015
#  [2,]   0.1423514  0.199277344   1.7467222   0.9323294  0.23106493  -2.1548154 -1.18028580  -2.2829 -1.037907 3.063288 -3.367614  1.00458644 -0.62556819  -1.3907173  0.09081271  -0.7354406   -2.308343   -3.812062   1.2871094  -0.9789022    2.625750  1.3370667
#  [3,]  -2.1028704 -0.008242187   1.8823559   0.2789206 -0.02461628   0.6292758 -0.09332919  -1.1579  1.212093 2.388288 -3.026704  0.47149796 -0.45511395  -2.1546062 -0.40918705  -3.4058947   -1.656169   -4.601929   1.0800258  -0.7245890    2.832833  1.4515076
#  [4,]  -2.9312047 -0.496523437   1.8710531   1.4436932 -1.30302516  -0.5070880  0.83058272  -1.7329  2.805843 1.963288 -2.771023 -0.09835401 -0.45511426  -0.4416436 -1.05918738  -1.6729400   -3.259430   -4.996862   0.9928327  -0.3482056    2.985421  2.4624023
#  [5,]  -1.9066860 -1.070253906   0.4242938   0.3925568 -0.05302569  -1.9275427  0.28710562  -1.2579  2.180843 2.413288 -2.430113  0.89429272  0.11306792  -2.7101618 -2.15918748  -3.8320311   -2.145299   -4.287777   1.5813860  -0.4194132    2.909127  2.1572266
#  [6,]  -1.7976946 -0.838320313  -0.6833813   1.7846025 -1.44507094   0.2599574 -0.55528636  -1.4579  1.680843 2.163288 -2.458523 -0.22703111 -0.39829576  -1.2055321 -2.08418747  -3.4058947   -1.547473   -2.896535   1.5050921  -0.3176880    2.909127  1.1844788
#  [7,]  -2.7132220 -0.630800781  -0.4008111   2.6936935  1.25379197  -0.3082244 -0.58246027  -2.2079  1.649593 3.638288 -5.242614  0.19576243  0.08465851  -1.9231247 -1.05918738  -3.3206674   -3.368125   -2.420819   1.4614955  -0.9890747    3.007220  1.4896545
#  [8,]  -2.0374756 -1.009218750   0.2208433   0.9607386 -0.84848030  -1.1320883  0.55884477  -1.9829  2.993342 3.138288 -2.941477  0.78399798 -0.68238669  -1.9231252 -0.88418737  -3.1786220   -3.368125   -2.950389   2.1808384  -1.0704549    3.029018  0.9174500
#  [9,]  -1.8412912 -1.021425781   1.2041875   1.5005117 -0.62120787   0.1463210 -0.47376340  -1.4829  2.493343 3.013288 -2.856250 -0.11673718 -1.19375006  -0.9971992  0.81581277  -0.7922581   -1.764864   -3.973626   1.7230748  -0.3380330    1.873710  0.4978333
# [10,]  -0.8603690 -1.106875000   1.4076380   3.2334659 -1.98484372  -0.3650426 -1.12593797  -1.7829  3.399593 2.513288 -2.714205 -0.24541366 -0.88124972  -1.8536808 -1.55918743  -4.3149857   -3.938777   -3.740256   0.9928327  -0.5007935    2.538557  1.3561401


tmp <- summary2[['6OHDA']][1:10, , drop = FALSE]
colnames(tmp) <- abbreviate(colnames(tmp), minlength = 8)
print(tmp)

#       s010423_000 s010423_009 e030723_003 e030723_008 e030723_01   e030923_   e041223_ e041823_003 e041823_005   e050223_ e122822_001 e122822_004 e123022_001 e123022_006 e123022_01 s221230_01 s221230_004 s221230_006 s221230_007 s122822_000 s122822_001 s122822_006 s122822_009
#  [1,]  -0.5033140    2.470772 -0.85825470    4.481100  -9.111440  1.6456535  0.4536148    1.624038 -0.12518090 -0.1504577   2.6384374   -4.869520  0.47190348   -4.737946  -2.397738 -1.2268066  -1.2727475    6.326751  0.37821452  -0.3135986    4.393389   -1.684262    6.432359
#  [2,]  -0.2659550    3.510277 -0.59783801    2.280014  -9.356976  0.4524719 -0.9451922    4.532692 -0.20670264 -0.1876598   1.6562945   -4.557020 -0.06786866   -5.303422  -2.427500 -1.9219293  -2.3209600    6.501137  0.41890462  -1.6399395    4.709464   -2.922916    6.678847
#  [3,]   0.6947835    3.662865  0.44382932    5.024579  -7.861439 -0.5702554  0.2155228    3.763462 -0.32445638 -0.5224813   1.2098655   -3.667597 -2.02809673   -4.708184  -1.058453 -1.1589898  -3.3824410    6.658085  1.49719238  -2.8723802    5.101833   -2.527983    5.775057
#  [4,]   1.1695014    4.273216  0.02716147    1.790884 -10.517690  0.1115628  0.7512345    3.547115 -0.01648554 -0.4034336   1.0908179   -3.763751  0.92644897   -3.606994  -3.588215 -0.2265082  -3.7804963    7.242278  1.95495605  -2.9545429    4.415187   -1.432940    5.505094
#  [5,]   1.5764025    4.912178  0.36570374    3.584361 -10.897154  2.7251992 -0.2606684    2.537500 -0.11612333 -0.7159338   1.1801036   -3.499328  0.50031257   -5.749851  -1.713214 -1.9388835  -2.5730617    7.957261  2.09737142  -3.8583327    3.009199   -1.594504    5.728107
#  [6,]   1.9041838    4.988472 -0.07700463    4.345231  -9.758761  0.5945171  0.6917114    2.585577 -0.13423907 -0.3141480   1.9836751   -3.307020  0.47190348   -6.166517  -3.498929 -1.4302572  -0.9543032    7.747998  1.89392090  -4.4686843    2.562334   -2.258710    4.929955
#  [7,]   1.3616491    4.072945 -1.43117141    2.660448  -8.620368  0.3956534  0.3643315    2.681731 -0.06177540  0.3034115   0.2872468   -3.162788 -0.38036932   -3.696280  -3.469167 -0.7012261  -0.1714610    8.053174 -0.04903158  -3.4944693    3.150888   -1.774019    4.953430
#  [8,]   1.2938323    2.661507 -0.62387968    3.665883  -8.040011  1.1342899 -0.7368575    3.354808 -0.53278993  0.1843637   0.5551039   -3.523365  0.47190285   -4.559375  -3.677500 -2.4136013   0.4388905    8.340911 -0.22196452  -2.8254301    3.968323   -2.115098    4.530879
#  [9,]   0.7286919    1.822273  0.10528647    2.035448  -9.915011  1.7024720 -0.3797166    4.388462  0.07409428  0.3182924   0.8229614   -2.129135 -1.14741547   -4.946280  -2.070357 -2.8544108   0.2531314    7.582331  0.83597819  -3.8700702    5.036438   -2.007389    3.720990
# [10,]   1.1921070    1.898567  0.52195317    1.410448  -8.352511  0.8786080 -0.2011427    3.667308 -0.17047116 -0.2323028   1.0908182   -2.393558 -1.03377783   -5.481994  -1.177500 -1.8710666   0.1204462    8.593771  1.12080892  -4.1869835    5.472403   -2.258710    3.767940

options(width = old_width)


save <- TRUE

if (save) {
  for (iii in seq_along(summary2)) {
    out <- summary2[[iii]]

    if (is.null(out)) {
      next
    }

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
