# ==============================================
# IMPORTS ABF FILE
# Analyses ABF file and creates all outputs for summaries and subsequent analysis
# ==============================================

# remove all objects from the environment
rm(list = ls(all = TRUE))

save <- FALSE

# load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp', 'reticulate', 'robustbase', 'signal', 'readABF', 'yaml')
load_required_packages(required.packages)

library(reticulate)
reticulate::py_config()

graphics.off()

# load user config which contains local paths
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)
print(config$username)

username <- config$username
path_repository <- config$path_repository
ABF2NWB_repository <- config$ABF2NWB_repository
path_analysis <- config$path_analysis

# construct local paths and identifiers
identifier1 <- 'Figure 11'
file_path1 <- paste0('/Users/', username, path_repository)
file_path2 <- paste0('/Users/', username, path_analysis)
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')
file_path5 <- paste0('/Users/', username, ABF2NWB_repository)

# source(paste0(file_path1, '/nNLS functions.R'))
source(file.path(file_path5, 'ABF2NWB_functions.R'))

# path where all graphs are stored
xlsx_path <- paste0(file_path4, 'xlsx')
if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

svg_path <- paste0(file_path4, 'svg')
if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}

setwd(file_path3)
expt_id <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)

# use most up to date experimental_dict WITH traces2average specified; batch = TRUE and silent = FALSE
experimental_dict <- list(
  'ChI-NGF control' = list( #
    '24d12009' = list(traces2average = list(c(5,6,7,8,9), c(12,13,14,15,16)),                           indices = c(10,17),     levels = c('control','AP5+NBQX')), 
    '24d12013' = list(traces2average = list(c(3,4,5,6,7), c(22,23,24,25,26)),                           indices = c(8,27),      levels = c('control','Mecamylamine')), 
    '25304011' = list(traces2average = list(c(15,17,18,19), c(42,43,44,45)),                            indices = c(20,46),     levels = c('control','Mecamylamine')), 
    '25306008' = list(traces2average = list(c(4,5,6,7,8), c(30,31,32,33,34), c(41,42,43,44,45)),        indices = c(9,35,46),   levels = c('control','Mecamylamine','Scopolamine')), 
    '25306003' = list(traces2average = list(c(1,2,3,4,5), c(26,27,28,29,30)),                           indices = c(6,31),      levels = c('control','Mecamylamine')), 
    '25409011' = list(traces2average = list(c(6,7,8,9,10)),                                             indices = c(11),        levels = c('control')), 
    '25409015' = list(traces2average = list(c(9,10,11,12,13)),                                          indices = c(14),        levels = c('control'))
  ),
  'ChI-NGF 6OHDA' = list(
    '25310000' = list(traces2average = list(c(6,7,8,9,10)),                                             indices = c(11),        levels = c('control')), 
    '25310001' = list(traces2average = list(c(8,9,10,11), c(31,32,33,34)),                              indices = c(13,35),     levels = c('control','Mecamylamine')), 
    '25310005' = list(traces2average = list(c(7,8,9,10), c(21,22,23,24,25)),                            indices = c(11,26),     levels = c('control','Mecamylamine')), 
    '25310006' = list(traces2average = list(c(5,6,7,8,9)),                                              indices = c(10),        levels = c('control')), 
    '25310007' = list(traces2average = list(c(12,13,14,15,16), c(33,34,35,36,37)),                      indices = c(17,38),     levels = c('control','Mecamylamine')), 
    '25310012' = list(traces2average = list(c(5,6,7,8,9), c(24,25,26,27,28)),                           indices = c(10,29),     levels = c('control','Mecamylamine')), 
    '25311000' = list(traces2average = list(c(13,14,15,16,17), c(21,22,23,24,25)),                      indices = c(18,26),     levels = c('control','Mecamylamine')), 
    '25311002' = list(traces2average = list(c(5,6,7,8,9), c(14,15,16,17,18), c(22,23,24,25,26)),        indices = c(10,19,27),  levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25311007' = list(traces2average = list(c(3,4,6,7), c(23,24,25,26,27), c(41,42,43,44,45)),          indices = c(8,28,46),   levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25312000' = list(traces2average = list(c(9,10,11,12,13), c(24,25,26,27), c(38,39,40,41)),          indices = c(14,29,42),  levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25312010' = list(traces2average = list(c(1,2,3,4,5), c(8,9,10,11,12), c(18,19,20,21,22)),          indices = c(6,13,23),   levels = c('control','Mecamylamine','AP5+NBQX')) 
  )
)

force_batch <- TRUE

ABF_summary <- ABF_batch_analysis(experimental_dict, file_path2, file_path3, silent=TRUE, force_batch=force_batch, save=save)

groups <- names(experimental_dict)
n_counts <- sapply(experimental_dict, length)

animal_counts <- sapply(experimental_dict, function(x) {
  ids <- substr(names(x), 1, 5)
  length(unique(ids))
})

n_counts
# ChI-NGF control   ChI-NGF 6OHDA 
#               7              11 animal_counts
animal_counts
# ChI-NGF control   ChI-NGF 6OHDA 
#               4               3 


###################################################################################################
###################################################################################################

# ==============================================
# EXTRACT SUMMARY DATA FROM ABF OUTPUT FILES
# Extracting data for fits to XLSX files
# ==============================================
# remove all objects from the environment
rm(list = ls(all = TRUE))

save <- FALSE

# load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp', 'reticulate', 'robustbase', 'signal', 'readABF', 'yaml')
load_required_packages(required.packages)

graphics.off()

# load user config which contains local paths
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)
print(config$username)

username <- config$username
path_repository <- config$path_repository
ABF2NWB_repository <- config$ABF2NWB_repository
path_analysis <- config$path_analysis

# construct local paths and identifiers
identifier1 <- 'Figure 11'
file_path1 <- paste0('/Users/', username, path_repository)
file_path2 <- paste0('/Users/', username, path_analysis)
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data summaries/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')
file_path5 <- paste0('/Users/', username, ABF2NWB_repository)

# source(paste0(file_path1, '/nNLS functions.R'))
source(file.path(file_path5, 'ABF2NWB_functions.R'))

# path where all graphs are stored
xlsx_path <- paste0(file_path4, 'xlsx')
if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

svg_path <- paste0(file_path4, 'svg')
if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}

# extract baseline adjusted responses for analysis
setwd(file_path3)
expt_id <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)

summary <- lapply(seq_along(expt_id), function(iii) { 
  
  file_path5 <- file.path(file_path3, expt_id[iii], 'xlsx')
  setwd(file_path5)
  folders <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)
  
  if (length(folders) == 0) {
    out <- list()
  } else {
    out <- lapply(seq_along(folders), function(ii) {
      folder <- folders[ii]
      folder_path <- file.path(file_path5, folder)
      setwd(folder_path)
      out1 <- load_data2(wd = folder_path, name = folder, header = TRUE)
      setwd(file_path5)  
      out1$`single examples`[-1]
    })
    names(out) <- folders
  }  
  out

})

names(summary) <- expt_id
  
# save data in separate sheets
if (save){
  invisible(
    sapply(seq_along(names(summary)), function(ii) list2excel(summary[[ii]], paste0(names(summary)[[ii]], '_full_abf.xlsx'), wd = xlsx_path, center_align = TRUE) )
    )
}

# save control traces in single spreadsheet 
summary2 <- lapply(seq_along(summary), function(iii){
  out <- sapply(seq_along(summary[[iii]]), function(ii){
    matrix(summary[[iii]][[ii]][,'control'], ncol=1)
  })
  colnames(out) <- names(summary[[iii]]) 
  out
})
names(summary2) <- expt_id

expt_id
# [1] "ChI-NGF 6OHDA"   "ChI-NGF control"

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
    addStyle(wb, sheet = 'Sheet1', style = center_style,
             rows = 1:(nrow(out) + 1), cols = 1:ncol(out),
             gridExpand = TRUE)

    saveWorkbook(wb, file = file.path(xlsx_path, paste0(names(summary2)[iii], '_abf.xlsx')), overwrite = TRUE)
  }
}

###################################################################################################
###################################################################################################

# ==============================================
# CREATE NWB FILES FOR DANDISET REPOSITORY
# converts ABF to NWB format and adds metadata
# ==============================================

# remove all objects from the environment
rm(list = ls(all = TRUE))

# load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp', 'reticulate', 'robustbase', 'signal', 'readABF', 'yaml')
load_required_packages(required.packages)

env_name <- "NWBenv"
if (!env_name %in% conda_list()$name) {
  conda_create(env_name, python_version = "3.11")  # pin Python version
  conda_install(env_name, packages = c("pynwb", "numpy"), pip = TRUE)
}
use_condaenv(env_name, required = TRUE)

py_config()  # show NWBenv path

graphics.off()

# load user config which contains local paths
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)
print(config$username)

username <- config$username
path_repository <- config$path_repository
ABF2NWB_repository <- config$ABF2NWB_repository
path_analysis <- config$path_analysis

identifier1 <- 'Figure 11'
file_path5 <- paste0('/Users/', username, ABF2NWB_repository)
source(file.path(file_path5, 'ABF2NWB_functions.R'))

# Define experimental data
experimental_dict <- list(
  'ChI-NGF control' = list( #
    '24d12009' = list(traces2average = list(c(5,6,7,8,9), c(12,13,14,15,16)),                           indices = c(10,17),     levels = c('control','AP5+NBQX')), 
    '24d12013' = list(traces2average = list(c(3,4,5,6,7), c(22,23,24,25,26)),                           indices = c(8,27),      levels = c('control','Mecamylamine')), 
    '25304011' = list(traces2average = list(c(15,17,18,19), c(42,43,44,45)),                            indices = c(20,46),     levels = c('control','Mecamylamine')), 
    '25306008' = list(traces2average = list(c(4,5,6,7,8), c(30,31,32,33,34), c(41,42,43,44,45)),        indices = c(9,35,46),   levels = c('control','Mecamylamine','Scopolamine')), 
    '25306003' = list(traces2average = list(c(1,2,3,4,5), c(26,27,28,29,30)),                           indices = c(6,31),      levels = c('control','Mecamylamine')), 
    '25409011' = list(traces2average = list(c(6,7,8,9,10)),                                             indices = c(11),        levels = c('control')), 
    '25409015' = list(traces2average = list(c(9,10,11,12,13)),                                          indices = c(14),        levels = c('control'))
  ),
  'ChI-NGF 6OHDA' = list(
    '25310000' = list(traces2average = list(c(6,7,8,9,10)),                                             indices = c(11),        levels = c('control')), 
    '25310001' = list(traces2average = list(c(8,9,10,11), c(31,32,33,34)),                              indices = c(13,35),     levels = c('control','Mecamylamine')), 
    '25310005' = list(traces2average = list(c(7,8,9,10), c(21,22,23,24,25)),                            indices = c(11,26),     levels = c('control','Mecamylamine')), 
    '25310006' = list(traces2average = list(c(5,6,7,8,9)),                                              indices = c(10),        levels = c('control')), 
    '25310007' = list(traces2average = list(c(12,13,14,15,16), c(33,34,35,36,37)),                      indices = c(17,38),     levels = c('control','Mecamylamine')), 
    '25310012' = list(traces2average = list(c(5,6,7,8,9), c(24,25,26,27,28)),                           indices = c(10,29),     levels = c('control','Mecamylamine')), 
    '25311000' = list(traces2average = list(c(13,14,15,16,17), c(21,22,23,24,25)),                      indices = c(18,26),     levels = c('control','Mecamylamine')), 
    '25311002' = list(traces2average = list(c(5,6,7,8,9), c(14,15,16,17,18), c(22,23,24,25,26)),        indices = c(10,19,27),  levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25311007' = list(traces2average = list(c(3,4,6,7), c(23,24,25,26,27), c(41,42,43,44,45)),          indices = c(8,28,46),   levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25312000' = list(traces2average = list(c(9,10,11,12,13), c(24,25,26,27), c(38,39,40,41)),          indices = c(14,29,42),  levels = c('control','Mecamylamine','AP5+NBQX')), 
    '25312010' = list(traces2average = list(c(1,2,3,4,5), c(8,9,10,11,12), c(18,19,20,21,22)),          indices = c(6,13,23),   levels = c('control','Mecamylamine','AP5+NBQX')) 
  )
)


# Define metadata function for Figure 11
get_figure11_metadata <- function(experiment_name) {
  
  base_metadata <- list(
    species = 'Mus musculus',
    age = 'P80D',
    sex = 'M',
    virus_injection_site = 'dorsolateral striatum',
    experimenter = 'Belal, Marziyeh',
    institution = 'Northwestern University',
    lab = 'Surmeier Lab',
    experiment_description = 'Optogenetic activation of cholinergic interneurons during whole-cell voltage clamp recordings from NDNF+ neurogliaform interneurons to examine ChI-NGF functional connectivity',
    device_description = 'Axon Multiclamp 700B patch clamp amplifier'
  )
  
  if (experiment_name == 'ChI-NGF control') {
    base_metadata$cell_type <- 'NGF'
    base_metadata$cross <- 'ChAT-Cre X NDNF-Flp'
    base_metadata$genotype <- 'ChAT-Cre+/-; NDNF-Flp+/-'
    base_metadata$session_description <- 'Whole-cell voltage clamp recording from NDNF+ neurogliaform interneuron in control mouse during optogenetic stimulation of cholinergic interneurons'
    base_metadata$notes <- 'Whole-cell voltage clamp recording in P80-90 ChAT-Cre X NDNF-Flp mice. NDNF-expressing NGFs were visually identified via Flp-dependent mCherry expression. AAV carrying Cre-dependent Chronos was stereotaxically injected into the DLS along with AAV carrying Flp-dependent mCherry expression construct. Recordings performed 4 weeks post-injection. Recordings were filtered at 2 kHz, and digitized at a sampling rate of 10 kHz using Clampex 10.7 software.'
    base_metadata$virus <- 'AAV-Cre-dependent Chronos (AAV-DIO-Chronos) + AAV-Flp-dependent mCherry (AAV-fDIO-mCherry)'
    base_metadata$keywords <- c('striatum', 'voltage clamp', 'optogenetics', 'Chronos', 'cholinergic interneurons', 'neurogliaform interneurons', 'NGF', 'NDNF', 'nicotinic', 'muscarinic', 'control')
    
  } else if (experiment_name == 'ChI-NGF 6OHDA') {
    base_metadata$cell_type <- 'NGF'
    base_metadata$cross <- 'ChAT-Cre X NDNF-Flp'
    base_metadata$genotype <- 'ChAT-Cre+/-; NDNF-Flp+/-'
    base_metadata$session_description <- 'Whole-cell voltage clamp recording from NDNF+ neurogliaform interneuron in 6-OHDA lesioned mouse during optogenetic stimulation of cholinergic interneurons'
    base_metadata$notes <- 'Whole-cell voltage clamp recording in P80-90 ChAT-Cre X NDNF-Flp mice with 6-OHDA lesion. NDNF-expressing NGFs were visually identified via Flp-dependent mCherry expression. AAV carrying Cre-dependent Chronos was stereotaxically injected into the DLS along with AAV carrying Flp-dependent mCherry expression construct. Unilateral 6-OHDA injection into medial forebrain bundle (MFB) to lesion dopaminergic neurons. Recordings performed 4 weeks post-lesion. Recordings were filtered at 2 kHz, and digitized at a sampling rate of 10 kHz using Clampex 10.7 software.'
    base_metadata$virus <- 'AAV-Cre-dependent Chronos (AAV-DIO-Chronos) + AAV-Flp-dependent mCherry (AAV-fDIO-mCherry)'
    base_metadata$keywords <- c('striatum', 'voltage clamp', 'optogenetics', 'Chronos', 'cholinergic interneurons', 'neurogliaform interneurons', 'NGF', 'NDNF', 'nicotinic', 'muscarinic', '6-OHDA', 'dopamine depletion', 'Parkinson disease model')
  }
  
  return(base_metadata)
}

# Run conversion
convert_abf_to_nwb(
  experimental_dict = experimental_dict,
  username = username,
  path_analysis = path_analysis,
  identifier = identifier1,
  get_metadata_fn = get_figure11_metadata
)


###################################################################################################
###################################################################################################
# ==============================================
# EXTRACT SUMMARY DATA FROM NWB FILES
# Extracting data for fits to XLSX files
# ==============================================

# remove all objects from the environment
rm(list = ls(all = TRUE))

save <- FALSE

# load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp', 'reticulate', 'robustbase', 'signal', 'readABF', 'yaml')
load_required_packages(required.packages)

env_name <- "NWBenv"
if (!env_name %in% conda_list()$name) {
  conda_create(env_name, python_version = "3.11")  # pin Python version
  conda_install(env_name, packages = c("pynwb", "numpy"), pip = TRUE)
}
use_condaenv(env_name, required = TRUE)

py_config()  # show NWBenv path

graphics.off()

# load user config which contains local paths
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)
print(config$username)

username <- config$username
path_repository <- config$path_repository
ABF2NWB_repository <- config$ABF2NWB_repository
path_analysis <- config$path_analysis

# construct local paths and identifiers
identifier1 <- 'Figure 11'
file_path1 <- paste0('/Users/', username, path_repository)
file_path2 <- paste0('/Users/', username, path_analysis)
file_path3 <- paste0('/Users/', username, path_analysis, '/NWB data/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')
file_path5 <- paste0('/Users/', username, ABF2NWB_repository)

# source(paste0(file_path1, '/nNLS functions.R'))
source(file.path(file_path5, 'ABF2NWB_functions.R'))

# path where all graphs are stored
xlsx_path <- paste0(file_path4, 'xlsx')
if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

# Extract baseline adjusted responses from NWB files
setwd(file_path3)

# Define experiment IDs (group names) upfront
expt_id <- c('ChI-NGF control', 'ChI-NGF 6OHDA')

# Get all subject folders (now at root level for DANDI compliance)
all_folders <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)
all_folders <- all_folders[grepl('^sub-', all_folders)]  # Only get subject folders

# Sort like MacOS Finder: letters before numbers
sort_key <- all_folders
sort_key <- gsub('o', '01', sort_key)  # October -> 01
sort_key <- gsub('n', '02', sort_key)  # November -> 02  
sort_key <- gsub('d', '03', sort_key)  # December -> 03
all_folders <- all_folders[order(sort_key)]

# Extract cell IDs
cell_ids <- gsub('^sub-m', '', all_folders)

# Process all NWB files
summary_nwb_all <- lapply(seq_along(all_folders), function(ii) {
  folder <- all_folders[ii]
  
  # Construct NWB filepath
  nwb_filename <- paste0(folder, '_icephys.nwb')
  nwb_filepath <- file.path(file_path3, folder, nwb_filename)
  
  if (!file.exists(nwb_filepath)) {
    cat("WARNING: NWB file not found:", nwb_filepath, "\n")
    return(NULL)
  }
  
  # Load and extract averaged traces
  averaged_traces <- load_nwb_averages(nwb_filepath)
  averaged_traces
})

names(summary_nwb_all) <- cell_ids

# Initialize summary_nwb with named list using expt_id
summary_nwb <- setNames(vector("list", length(expt_id)), expt_id)

for (ii in seq_along(all_folders)) {
  folder <- all_folders[ii]
  cell_id <- cell_ids[ii]
  
  # Read session description from NWB file to determine group
  nwb_filepath <- file.path(file_path3, folder, paste0(folder, '_icephys.nwb'))
  
  py$nwb_filepath <- nwb_filepath
  session_desc <- py_run_string("
from pynwb import NWBHDF5IO
io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()
session_desc = str(nwb.session_description)
io.close()
", convert = TRUE)$session_desc
  
  # Assign to appropriate group based on session description
  if (grepl('6-OHDA', session_desc, ignore.case = TRUE)) {
    summary_nwb[[expt_id[2]]][[cell_id]] <- summary_nwb_all[[ii]]  # 6OHDA
  } else {
    summary_nwb[[expt_id[1]]][[cell_id]] <- summary_nwb_all[[ii]]  # Control
  }
}

# Print summary
cat("\nControl group:", length(summary_nwb[[expt_id[1]]]), "subjects\n")
cat("6-OHDA group:", length(summary_nwb[[expt_id[2]]]), "subjects\n")

# SAVE DATA IN SEPARATE SHEETS
setwd(file_path4)
invisible(
  sapply(seq_along(names(summary_nwb)), function(ii) {
    list2excel(summary_nwb[[ii]], paste0(names(summary_nwb)[[ii]], '_full.xlsx'), 
               wd = xlsx_path, center_align = TRUE)
  })
)

# Save control traces
summary2_nwb <- lapply(seq_along(summary_nwb), function(iii){
  out <- sapply(seq_along(summary_nwb[[iii]]), function(ii){
    matrix(summary_nwb[[iii]][[ii]][,'control'], ncol=1)
  })
  colnames(out) <- names(summary_nwb[[iii]]) 
  out
})
names(summary2_nwb) <- expt_id  

# Sort alphabetically
summary2_nwb <- summary2_nwb[order(names(summary2_nwb))]
   names(summary2_nwb)
[1] "ChI-NGF 6OHDA"   "ChI-NGF control"

summary2_nwb[[1]][1:10,]
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

summary2_nwb[[2]][1:10,]
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
  for (iii in seq_along(summary2_nwb)) {
    out <- summary2_nwb[[iii]]
    wb <- createWorkbook()
    addWorksheet(wb, 'Sheet1')
    writeData(wb, sheet = 'Sheet1', x = out, startRow = 1, startCol = 1, rowNames = FALSE)

    center_style <- createStyle(halign = 'center', valign = 'center')
    addStyle(wb, sheet = 'Sheet1', style = center_style,
             rows = 1:(nrow(out) + 1), cols = 1:ncol(out),
             gridExpand = TRUE)

    saveWorkbook(wb, file = file.path(xlsx_path, paste0(names(summary2_nwb)[iii], '.xlsx')), overwrite = TRUE)
  }
}


###################################################################################################
###################################################################################################

# ==============================================
# DANDI-SPECIFIC VALIDATION
# ==============================================

# use terminal for validation:

# install dandi
# pip install dandi nwbinspector

# validate with DANDI-specific checks with dandi validate or nwbinspector

# cd ~"$(grep path_analysis ~/.abf2nwb_config.yaml | cut -d"'" -f2)/NWB data/Figure 11"
# dandi validate 


# Validate with PyNWB directly
python -c "
from pynwb import validate, NWBHDF5IO
import os

for folder in os.listdir('.'):
    if folder.startswith('sub-'):
        nwb_file = os.path.join(folder, f'{folder}_icephys.nwb')
        if os.path.exists(nwb_file):
            with NWBHDF5IO(nwb_file, 'r') as io:
                errors = validate(io)
                status = '✓' if len(errors) == 0 else f'✗ ({len(errors)} errors)'
                print(f'{folder}: {status}')
"


###################################################################################################
###################################################################################################
# ==============================================
# READ METADATA FROM FIRST NWB FILE
# ==============================================

# remove all objects from the environment
rm(list = ls(all = TRUE))

# load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}

required.packages <- c('dbscan', 'jsonlite', 'minpack.lm', 'openxlsx', 'Rcpp', 'reticulate', 'robustbase', 'signal', 'readABF', 'yaml')
load_required_packages(required.packages)

env_name <- "NWBenv"
if (!env_name %in% conda_list()$name) {
  conda_create(env_name, python_version = "3.11")  # pin Python version
  conda_install(env_name, packages = c("pynwb", "numpy"), pip = TRUE)
}
use_condaenv(env_name, required = TRUE)

py_config()  # show NWBenv path

graphics.off()

# load user config which contains local paths
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)
print(config$username)

username <- config$username
path_repository <- config$path_repository
ABF2NWB_repository <- config$ABF2NWB_repository
path_analysis <- config$path_analysis

identifier1 <- 'Figure 11'

file_path3 <- paste0('/Users/', username, path_analysis, '/NWB data/', identifier1)
setwd(file_path3)
folders <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)

file_path5 <- paste0('/Users/', username, ABF2NWB_repository)

source(file.path(file_path5, 'ABF2NWB_functions.R'))

# Get first NWB file
first_folder <- folders[1]
nwb_filepath <- file.path(file_path3, first_folder, paste0(first_folder, '_icephys.nwb'))

# Pass filepath to Python
py$nwb_filepath <- nwb_filepath

# use print_nwb_metadata() to display metadata
print_nwb_metadata(nwb_filepath)

# === SESSION INFO ===
# Session description: Whole-cell voltage clamp recording from NDNF+ neurogliaform interneuron in control mouse during optogenetic stimulation of cholinergic interneurons. Notes: Whole-cell voltage clamp recording in P80-90 ChAT-Cre X NDNF-Flp mice. NDNF-expressing NGFs were visually identified via Flp-dependent mCherry expression. AAV carrying Cre-dependent Chronos was stereotaxically injected into the DLS along with AAV carrying Flp-dependent mCherry expression construct. Recordings performed 4 weeks post-injection. Recordings were filtered at 2 kHz, and digitized at a sampling rate of 10 kHz using Clampex 10.7 software. Original recording duration: 903.0 seconds.
# Session start time: 2024-12-12 20:33:04-06:00
# Experimenter: Belal, Marziyeh
# Institution: Northwestern University
# Lab: Surmeier Lab

# === SUBJECT INFO ===
# Subject ID: m24d12009
# Species: Mus musculus
# Age: P80D
# Sex: M
# Genotype: ChAT-Cre+/-; NDNF-Flp+/-
# Description: Cross: ChAT-Cre X NDNF-Flp | Virus: AAV-Cre-dependent Chronos (AAV-DIO-Chronos) + AAV-Flp-dependent mCherry (AAV-fDIO-mCherry) | Injection site: dorsolateral striatum

# === EXPERIMENT INFO ===
# Experiment description: Optogenetic activation of cholinergic interneurons during whole-cell voltage clamp recordings from NDNF+ neurogliaform interneurons to examine ChI-NGF functional connectivity
# Keywords: ['striatum' 'voltage clamp' 'optogenetics' 'Chronos'
#  'cholinergic interneurons' 'neurogliaform interneurons' 'NGF' 'NDNF'
#  'nicotinic' 'muscarinic' 'control']

# === DEVICES ===
# Multiclamp 700B Patch clamp amplifier:
#   Description: Axon Multiclamp 700B patch clamp amplifier

# === ELECTRODES ===
# elec0:
#   Cell ID: m24d12009
#   Description: Whole-cell patch-clamp electrode | Original unit: pA | Data stored in amperes | Headstage gain: 0.0005 V/pA | Additional gain: 5.00x | Total gain: 2.50 mV/pA
#   Location: Striatum
#   Resistance: 3-5 MOhm
#   Filtering: 5000.0 Hz low-pass

# === ACQUISITION DATA ===
# Number of sweeps: 10
#   Sweep 4: Voltage clamp recording - sweep 4 (control condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 5: Voltage clamp recording - sweep 5 (control condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 6: Voltage clamp recording - sweep 6 (control condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 7: Voltage clamp recording - sweep 7 (control condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 8: Voltage clamp recording - sweep 8 (control condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 11: Voltage clamp recording - sweep 11 (AP5+NBQX condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 12: Voltage clamp recording - sweep 12 (AP5+NBQX condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 13: Voltage clamp recording - sweep 13 (AP5+NBQX condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 14: Voltage clamp recording - sweep 14 (AP5+NBQX condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes
#   Sweep 15: Voltage clamp recording - sweep 15 (AP5+NBQX condition)
#     Rate: 10000.0 Hz
#     Gain: 2.500000
#     Unit: amperes

# === PROCESSING MODULES ===
# Module: icephys
#   baseline_duration: [100.]
#   level_column_mapping: {'control': [0, 1, 2, 3, 4], 'AP5+NBQX': [5, 6, 7, 8, 9]}
#   samplingIntervalInSec: [1.e-04]
#   stimulation_time: [346.9]
#   traces2save: [ 4  5  6  7  8 11 12 13 14 15]
#   traces_for_averaging: [0 1 2 3 4 5 6 7 8 9]

# === FILE INFO ===
# Identifier: NWB_m24d12009
# File created: 2026-04-01 22:26:45 