# ==============================================
# IMPORTS ABF FILE
# Analyses ABF file and creates all outputs for summaries and subsequent analysis
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
identifier2 <- 'Figure 11 full'
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

# retrieve the names of all folders (directories) in the working directory
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data/', identifier2)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier2, '/')

setwd(file_path3)
expt_id <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)

full_experimental_dict <- list(
  'ChI-NGF control' = list( #
    '24d09004' = list( # OMIT: RUNDOWN no baseline period; HP - 69.1mV
      indices = c(19, 53),     
      levels = c('control', 'Mecamylamine')
      ),
    '24d12009' = list( # KEEP; HP -61.3mV
      indices = c(10,17),     
      levels = c('control', 'AP5+NBQX')
      ),
    '24d12013' = list( # KEEP; HP -70.2mV 
      indices = c(8, 27),     
      levels = c('control', 'Mecamylamine')
      ),
    '25304011' = list( # KEEP; HP -64.3mV 
      indices = c(20, 46),     
      levels = c('control', 'Mecamylamine')
      ),
    '25306008' = list( # KEEP; HP -81.1mV 
      indices = c(9, 35, 46),        
      levels = c('control', 'Mecamylamine', 'Scopolamine')
      ),
    '25306003' = list( # KEEP; HP -72.2mV
      indices = c(6, 31),     
      levels = c('control', 'Mecamylamine')
      ),
    '25409011' = list( # KEEP; HP -80.0mV 
      indices = c(11),     
      levels = c('control')
      ),
    '25409015' = list( # KEEP; HP -84.1mV
      indices = c(14),     
      levels = c('control')
      )                             
  ),
  'ChI-NGF 6OHDA' = list(
    '25310000' = list( # KEEP; HP -81.1mV 
      indices = c(11),
      levels = c('control')
    ),
    '25310001' = list( # KEEP; HP -71.1mV 
      indices = c(13, 35),
      levels = c('control', 'Mecamylamine')
    ),
    '25310005' = list( # KEEP; ? stability; HP -71.1mV 
      indices = c(11, 26),
      levels = c('control', 'Mecamylamine')
    ),
    '25310006' = list( # KEEP; HP -71.1mV  
      indices = c(10),
      levels = c('control')
    ),
    '25310007' = list( # KEEP; HP -71.1mV 
      indices = c(17, 38),
      levels = c('control', 'Mecamylamine')
    ),
    '25310012' = list( # KEEP; HP -81.0mV 
      indices = c(10, 29),
      levels = c('control', 'Mecamylamine')
    ),
    '25311000' = list( # KEEP; HP -82.0mV 
      indices = c(18, 26),
      levels = c('control', 'Mecamylamine')
    ),
    '25311002' = list( # KEEP; HP -82.0mV 
      indices = c(10, 19, 27),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25311007' = list( # KEEP; HP -81.0mV 
      indices = c(8, 28, 46),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25312000' = list( # KEEP; HP - 81.0mV
      indices = c(14, 29, 42),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25312002' = list( # OMIT: leaky; HP -80.0mV
      indices = c(11, 23, 36),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25312010' = list( # KEEP; HP -72.1mV 
      indices = c(6, 13, 23),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    )
  )
)

save <- FALSE
force_batch <-TRUE

width <- 6 * 0.7
height <- 8 * 0.7

graph_settings <- list(width=width, height=height, xlim=NA, ylim1=NA, ylim2=NA, xlab='time (minutes)',
    ylab='|PSC| (pA)', colors='#CD5C5C', xmajor_tick=10, ymajor_tick=100, cex.points=0.8)

ABF_summary <- ABF_batch_analysis(full_experimental_dict, file_path2, file_path3, graph_settings=graph_settings, force_batch=force_batch, save=save)

expt_id <- names(ABF_summary)
HP <- setNames(lapply(1:length(ABF_summary), function(iii){
  sapply(1:length(ABF_summary[[iii]]), function(ii) mean(ABF_summary[[iii]][[ii]]$summary[, 'holding potential (mV)']))
  }), expt_id)

# cycle through all to extract those gains
metadata_summary <- setNames(lapply(expt_id, function(folder) {
  full_folder <- file.path('.', folder)
  # sub_folders <- list.dirs(path = full_folder, full.names = FALSE, recursive = FALSE)
  sub_folders <- names(full_experimental_dict[[folder]])
  sub_result <- setNames(lapply(sub_folders, function(sub_folder) {
    full_sub_folder <- file.path(full_folder, sub_folder)
    file_name <- list.files(full_sub_folder, pattern = '\\.abf$', full.names = FALSE)[1]

    if (!is.na(file_name)) {
      out <- readABF(file.path(full_sub_folder, file_name))

      dig1 <- get_all_digitisation(out$data, col = 1)
      dig2 <- get_all_digitisation(out$data, col = 2)
      out$digitisation <- list(col1 = dig1, col2 = dig2)

      out[names(out) != 'data']
    } else {
      NULL
    }
  }), sub_folders)

  sub_result
}), expt_id)

adc_scaling_summary_all <- lapply(metadata_summary, function(sublist) {
  lapply(sublist, function(meta) {
    get_adc_scaling_summary(meta)
  })
})

adc_scaling_summary_all
# $`ChI-NGF control`
# $`ChI-NGF control`$`24d09004`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`24d12009`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`24d12013`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`25304011`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`25306008`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`25306003`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`25409011`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF control`$`25409015`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE


# $`ChI-NGF 6OHDA`
# $`ChI-NGF 6OHDA`$`25310000`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25310001`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25310005`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25310006`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25310007`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25310012`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25311000`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25311002`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25311007`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25312000`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25312002`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE

# $`ChI-NGF 6OHDA`$`25312010`
#   channel units telegraph enabled scale factor (mV/unit) additional gain total gain (mV/unit) gain units AD min AD max digitisation trace digitisation recorded
# 1       0    pA              TRUE                    0.5               5                  2.5      mV/pA  -4000   4000   0.12207031         0.12207031     TRUE
# 2       1    mV              TRUE                   10.0               2                 20.0      mV/mV   -500    500   0.01525879         0.01525879     TRUE


# setwd(paste0(file_path3, '/NPY Cre X iSPN WT/23404012'))
# setwd(paste0(file_path3, '/NPY Cre X iSPN WT/23d20000'))

# abf_files <- list.files(path='.', pattern = '\\.[Aa][Bb][Ff]$', full.names=FALSE)
# out <- readABFs(abf_files)

HP
# $`ChI-NGF control`
# [1] -69.08234 -61.26321 -70.17874 -64.32392 -81.01107 -72.16824 -80.01010
# [8] -84.06792

# $`ChI-NGF 6OHDA`
#  [1] -81.08182 -71.13773 -71.13724 -71.11176 -71.09311 -80.96605 -82.02505
#  [8] -82.03046 -81.00710 -81.03669 -80.02735 -72.12288

############
# retrieve the names of all folders (directories) in the working directory
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')      
setwd(file_path3)

experimental_dict <- list(
  'ChI-NGF control' = list( #
    '24d12009' = list( # KEEP; HP -61.3mV
      indices = c(10,17),     
      levels = c('control', 'AP5+NBQX')
      ),
    '24d12013' = list( # KEEP; HP -70.2mV 
      indices = c(8, 27),     
      levels = c('control', 'Mecamylamine')
      ),
    '25304011' = list( # KEEP; HP -64.3mV 
      indices = c(20, 46),     
      levels = c('control', 'Mecamylamine')
      ),
    '25306008' = list( # KEEP; HP -81.1mV 
      indices = c(9, 35, 46),        
      levels = c('control', 'Mecamylamine', 'Scopolamine')
      ),
    '25306003' = list( # KEEP; HP -72.2mV
      indices = c(6, 31),     
      levels = c('control', 'Mecamylamine')
      ),
    '25409011' = list( # KEEP; HP -80.0mV 
      indices = c(11),     
      levels = c('control')
      ),
    '25409015' = list( # KEEP; HP -84.1mV
      indices = c(14),     
      levels = c('control')
      )                             
  ),
  'ChI-NGF 6OHDA' = list(
    '25310000' = list( # KEEP; HP -81.1mV 
      indices = c(11),
      levels = c('control')
    ),
    '25310001' = list( # KEEP; HP -71.1mV 
      indices = c(13, 35),
      levels = c('control', 'Mecamylamine')
    ),
    '25310005' = list( # KEEP; ? stability; HP -71.1mV 
      indices = c(11, 26),
      levels = c('control', 'Mecamylamine')
    ),
    '25310006' = list( # KEEP; HP -71.1mV  
      indices = c(10),
      levels = c('control')
    ),
    '25310007' = list( # KEEP; HP -71.1mV 
      indices = c(17, 38),
      levels = c('control', 'Mecamylamine')
    ),
    '25310012' = list( # KEEP; HP -81.0mV 
      indices = c(10, 29),
      levels = c('control', 'Mecamylamine')
    ),
    '25311000' = list( # KEEP; HP -82.0mV 
      indices = c(18, 26),
      levels = c('control', 'Mecamylamine')
    ),
    '25311002' = list( # KEEP; HP -82.0mV 
      indices = c(10, 19, 27),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25311007' = list( # KEEP; HP -81.0mV 
      indices = c(8, 28, 46),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25312000' = list( # KEEP; HP - 81.0mV
      indices = c(14, 29, 42),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    ),
    '25312010' = list( # KEEP; HP -72.1mV 
      indices = c(6, 13, 23),
      levels = c('control', 'Mecamylamine', 'AP5+NBQX')
    )
  )
)

# retrieve the names of all folders (directories) in the working directory
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')      
setwd(file_path3)

width <- 6 * 0.7
height <- 8 * 0.7

graph_settings <- list(width=width, height=height, xlim=NA, ylim1=NA, ylim2=NA, xlab='time (minutes)',
    ylab='|PSC| (pA)', colors='#CD5C5C', xmajor_tick=10, ymajor_tick=100, cex.points=0.8)

force_batch <-FALSE
save <- FALSE

ABF_summary <- ABF_batch_analysis(experimental_dict, file_path2, file_path3, graph_settings=graph_settings, force_batch=force_batch, save=save)

for (expt in names(ABF_summary)) {
  for (folder in names(ABF_summary[[expt]])) {
    result <- ABF_summary[[expt]][[folder]]
    if (!is.null(result$'final traces2average')) {
      experimental_dict[[expt]][[folder]]$traces2average <- result$'final traces2average'
    }
  }
}

format_experimental_dict(experimental_dict)

  #   '24d12009' = list(traces2average = list(c(5,6,7,8,9), c(12,13,14,15,16)),    indices = c(10,17),    levels = c('control','AP5+NBQX')), 
  #   '24d12013' = list(traces2average = list(c(3,4,5,6,7), c(22,23,24,25,26)),    indices = c(8,27),    levels = c('control','Mecamylamine')), 
  #   '25304011' = list(traces2average = list(c(15,17,18,19), c(42,43,44,45)),    indices = c(20,46),    levels = c('control','Mecamylamine')), 
  #   '25306008' = list(traces2average = list(c(4,5,6,7,8), c(30,31,32,33,34), c(41,42,43,44,45)),    indices = c(9,35,46),    levels = c('control','Mecamylamine','Scopolamine')), 
  #   '25306003' = list(traces2average = list(c(1,2,3,4,5), c(26,27,28,29,30)),    indices = c(6,31),    levels = c('control','Mecamylamine')), 
  #   '25409011' = list(traces2average = list(c(6,7,8,9,10)),    indices = c(11),    levels = c('control')), 
  #   '25409015' = list(traces2average = list(c(9,10,11,12,13)),    indices = c(14),    levels = c('control')) 
  # ),
  #   '25310000' = list(traces2average = list(c(6,7,8,9,10)),    indices = c(11),    levels = c('control')), 
  #   '25310001' = list(traces2average = list(c(8,9,10,11), c(31,32,33,34)),    indices = c(13,35),    levels = c('control','Mecamylamine')), 
  #   '25310005' = list(traces2average = list(c(7,8,9,10), c(21,22,23,24,25)),    indices = c(11,26),    levels = c('control','Mecamylamine')), 
  #   '25310006' = list(traces2average = list(c(5,6,7,8,9)),    indices = c(10),    levels = c('control')), 
  #   '25310007' = list(traces2average = list(c(12,13,14,15,16), c(33,34,35,36,37)),    indices = c(17,38),    levels = c('control','Mecamylamine')), 
  #   '25310012' = list(traces2average = list(c(5,6,7,8,9), c(24,25,26,27,28)),    indices = c(10,29),    levels = c('control','Mecamylamine')), 
  #   '25311000' = list(traces2average = list(c(13,14,15,16,17), c(21,22,23,24,25)),    indices = c(18,26),    levels = c('control','Mecamylamine')), 
  #   '25311002' = list(traces2average = list(c(5,6,7,8,9), c(14,15,16,17,18), c(22,23,24,25,26)),    indices = c(10,19,27),    levels = c('control','Mecamylamine','AP5+NBQX')), 
  #   '25311007' = list(traces2average = list(c(3,4,6,7), c(23,24,25,26,27), c(41,42,43,44,45)),    indices = c(8,28,46),    levels = c('control','Mecamylamine','AP5+NBQX')), 
  #   '25312000' = list(traces2average = list(c(9,10,11,12,13), c(24,25,26,27), c(38,39,40,41)),    indices = c(14,29,42),    levels = c('control','Mecamylamine','AP5+NBQX')), 
  #   '25312010' = list(traces2average = list(c(1,2,3,4,5), c(8,9,10,11,12), c(18,19,20,21,22)),    indices = c(6,13,23),    levels = c('control','Mecamylamine','AP5+NBQX')), 
  # ),

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
# retrieve the names of all folders (directories) in the working directory
file_path3 <- paste0('/Users/', username, path_analysis, '/Raw ABF data/', identifier1)
file_path4 <- paste0('/Users/', username, path_analysis, '/',  identifier1, '/')      
setwd(file_path3)

force_batch <- TRUE
save <- TRUE

ABF_summary <- ABF_batch_analysis(experimental_dict, file_path2, file_path3, silent=TRUE, force_batch=force_batch, save=save)

groups <- names(experimental_dict)
n_counts <- sapply(experimental_dict, length)

animal_counts <- sapply(experimental_dict, function(x) {
  ids <- substr(names(x), 1, 5)
  length(unique(ids))
})

n_counts
# ChI-NGF control   ChI-NGF 6OHDA 
#               7              11 
animal_counts
# ChI-NGF control   ChI-NGF 6OHDA 
#               4               3 




