# ==============================================
# EXTRACT FIGURE 3 SUMMARY DATA FROM DANDI NWB FILES
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

identifier <- 'Figure 3'
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
if (is.null(dataset_description$FigureMappings[['Figure 3']])) {
  stop('Figure mapping not found in: ', dataset_description_path)
}

fig3 <- do.call(
  rbind,
  lapply(dataset_description$FigureMappings[['Figure 3']], as.data.frame)
)

# FigureMappings[['Figure 3']] also includes Figure 3B entries, so keep only Figure 3.
fig3 <- subset(
  fig3,
  grepl('^Figure 3/', original_path) & grepl('_icephys\\.nwb$', dandi_path)
)

fig3$nwb_filepath <- file.path(dandi_root, fig3$dandi_path)
fig3$subject <- sub('/.*$', '', fig3$dandi_path)

missing_files <- fig3$nwb_filepath[!file.exists(fig3$nwb_filepath)]
if (length(missing_files) > 0) {
  warning('Missing Figure 3 NWB files:\n', paste(missing_files, collapse = '\n'))
}

fig3 <- fig3[file.exists(fig3$nwb_filepath), ]
if (nrow(fig3) == 0) {
  stop('No local Figure 3 icephys NWB files found under: ', dandi_root)
}

# Sort like Mac Finder: letters before numbers
sort_key <- fig3$subject
sort_key <- gsub('o', '01', sort_key)
sort_key <- gsub('n', '02', sort_key)
sort_key <- gsub('d', '03', sort_key)
fig3 <- fig3[order(sort_key), ]

cell_ids <- gsub('^sub-m', '', fig3$subject)

cat('\nFigure 3 NWB files found:', nrow(fig3), '\n')
print(fig3[, c('subject', 'dandi_path')])

# extract averages
expt_id <- c('CRISPR control', 'CRISPR delta KD')

summary_all <- lapply(seq_len(nrow(fig3)), function(ii) {
  nwb_filepath <- fig3$nwb_filepath[ii]
  cat('Loading:', nwb_filepath, '\n')
  load_nwb_averages(nwb_filepath)
})

names(summary_all) <- cell_ids

# group by virus metadata
summary <- setNames(vector('list', length(expt_id)), expt_id)

for (ii in seq_len(nrow(fig3))) {
  nwb_filepath <- fig3$nwb_filepath[ii]
  cell_id <- cell_ids[ii]

  py$nwb_filepath <- nwb_filepath

  virus_info <- py_run_string("
from pynwb import NWBHDF5IO
import re

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
description = str(nwb.subject.description)
virus_match = re.search(r'Virus: ([^|]+)', description)
virus = virus_match.group(1).strip() if virus_match else 'None'
io.close()
", convert = TRUE)$virus

  if (grepl('saCas', virus_info)) {
    summary[[expt_id[2]]][[cell_id]] <- summary_all[[ii]]
  } else {
    summary[[expt_id[1]]][[cell_id]] <- summary_all[[ii]]
  }
}

cat('\nCRISPR control:', length(summary[[expt_id[1]]]), 'subjects\n')
cat('CRISPR delta KD:', length(summary[[expt_id[2]]]), 'subjects\n')

summary <- summary[order(names(summary))]

# Convert lists to matrices
summary <- lapply(summary, function(group_data) {
  if (length(group_data) > 0) {
    out <- do.call(cbind, lapply(group_data, function(x) matrix(x[, 1], ncol = 1)))
    colnames(out) <- names(group_data)
    out
  } else {
    NULL
  }
})

save <- TRUE

if (save) {
  for (iii in seq_along(summary)) {
    out <- summary[[iii]]
    if (!is.null(out)) {
      wb <- createWorkbook()
      addWorksheet(wb, 'Sheet1')
      writeData(wb, sheet = 'Sheet1', x = out, startRow = 1, startCol = 1, rowNames = FALSE)

      center_style <- createStyle(halign = 'center', valign = 'center')
      addStyle(wb, sheet = 'Sheet1', style = center_style,
               rows = 1:(nrow(out) + 1), cols = 1:ncol(out),
               gridExpand = TRUE)

      saveWorkbook(wb, file = file.path(xlsx_path, paste0(names(summary)[iii], '.xlsx')), overwrite = TRUE)
    }
  }
}

summary[[1]][1:10,]
#      24o14002   24o14008   24o15004   24729002   24730002   24821001   24822001
# 1   2.4108886  1.5002441  5.2404783  8.3514400 -0.4781087  0.7369995   3.143921
# 2   0.1729329  1.5002441  4.6301267  2.8582762  1.5563964 -1.0940551 -10.283813
# 3  -1.0477701 -0.9411621 -3.3044432  5.9100339 -1.6988118  7.4508663   3.754272
# 4   3.6315916  5.7727048 -3.9147947 10.1824946 -4.5471189  4.7042844   1.312866
# 5  -0.4374186  0.2795410  3.4094237 -0.1934814  0.7425944  7.4508663  -7.842407
# 6   4.6488442 -3.3825682  8.2922359 -8.1280514  2.5736490 -1.7044067 -15.166625
# 7   3.8350421  2.7209471 -4.5251463  2.2479247 -3.3264159 10.8077998  -1.128540
# 8   2.6143391 -0.9411621 -0.8630371  4.0789793  0.1322428  4.3991087  -8.452758
# 9  -0.4374186 15.5383293 -0.2526855 -3.2452391 -3.7333169  0.4318237  -9.673461
# 10  1.5970865  1.5002441  8.9025875  6.5203854  4.4047036  6.5353390  -6.621704

summary[[2]][1:10,]
#      24523004   24523009   24524007   24524011   24530003   24530007   24531005  24613007   24614003   24614006
# 1   1.4967854 -5.0213621  1.8707275 -4.1111245 -6.4617917  1.2576293 -6.2473548 -1.520996 -3.0432127 -0.6980387
# 2  -1.7584228 -3.5972085  1.5655517 -2.2800699  1.6762288  4.9197385 -0.3472900 -6.098633 -1.6190592  4.1847736
# 3  -5.0136309 -2.5799559 -5.1483152 -4.7214760  0.4555257 -3.9303587 -0.7541910 -2.436523  1.0257975  0.5226644
# 4   1.2933349 -1.1558024 -4.2327879  3.0096434 -4.8341876 -0.8786010 -0.9576416 -4.877929 -2.6363117 -1.3083902
# 5   4.7519936 -5.4282631  3.0914305  1.5854898 -1.1720784 -9.1183467 -1.9748941 -4.572754  4.2810057 -1.5118407
# 6   0.2760823 -0.3420003 -1.4862060 -1.8731689  0.6589762 -0.2682495  2.0941161 -8.234863  1.0257975 -1.5118407
# 7   1.9036864 -1.1558024  8.5845943 -8.3835852 -0.7651774 -3.9303587 -2.1783446 -2.436523  1.4326985  2.1502685
# 8   3.7347410 -1.9696044 -4.5379637 -2.8904214 -0.9686279  3.6990355  4.7389728  2.141113  3.2637531  3.9813231
# 9  -0.5377197  0.6752522  0.3448486  0.7716878  3.3038329  3.3938597  0.4665120 -7.014160 -0.1949056 -1.3083902
# 10 -0.7411702 -2.1730549  5.2276609 -0.4490153 -5.4445391  0.3421020 -0.9576416 -7.014160  5.9086097  0.7261149
