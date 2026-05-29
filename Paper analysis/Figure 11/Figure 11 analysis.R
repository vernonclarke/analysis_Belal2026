# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

# Load required packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
load_required_packages(c('robustbase', 'minpack.lm', 'openxlsx', 'Rcpp', 'signal', 'dbscan'))

# Load user config
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)

# Construct paths
username <- config$username
file_path1 <- paste0('/Users/', username, config$path_repository)
file_path2 <- paste0('/Users/', username, config$path_analysis)

source(paste0(file_path1, '/nNLS functions.R'))

identifier <- 'Figure 11'
analysis_path <- paste0(file_path2, '/', identifier)
xlsx_path <- paste0(analysis_path, '/xlsx')
svg_path <- paste0(analysis_path, '/xlsx')

svg_path <- paste0(analysis_path, '/svg')
if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}

setwd(analysis_path)

dataset_names <- c('ChI-NGF control', 'ChI-NGF 6OHDA')

# data load
ChI_NGF_control_data <- load_data2(wd=xlsx_path, name=dataset_names[1] , header=TRUE)[[1]]
expt_ids <- colnames(ChI_NGF_control_data)

ChI_NGF_6OHDA_data <- load_data2(wd=xlsx_path, name=dataset_names[2] , header=TRUE)[[1]]
expt_ids <- colnames(ChI_NGF_6OHDA_data)


ChI_NGF_control_data[1:10,]
#       24d12009   24d12013   25304011    25306003    25306008    25409011   25409015
# 1  -0.48649900  0.6825439  1.1930847 -1.39470208 -0.01547851 1.034252881 -0.7939941
# 2  -0.36442869 -0.5381592  0.2165222 -0.56462400 -0.06430664 0.008862304 -0.2080566
# 3   0.17268066 -0.2696045 -0.2717590 -0.49138181 -0.84555660 0.008862304  0.4999511
# 4  -0.53532712 -0.8311279  0.2470398 -0.05192871 -0.45493162 0.423901347  0.3290527
# 5  -0.55974118  0.1210205  0.9794616 -0.63786618 -0.77231442 0.423901347 -0.3057129
# 6  -0.51091306 -0.3916748  1.4677429 -0.46696775 -0.57700193 0.448315408 -1.6729003
# 7   0.07502441 -0.2207764  1.3151550  0.21662597 -0.21079101 0.375073224 -1.1113769
# 8   0.49006345  0.1942627  0.9794616  0.33869627 -0.60141599 0.204174795 -1.3066894
# 9   0.09943847 -0.5381592  0.1860046  0.77814938 -0.67465817 0.985424758  0.3046387
# 10 -0.70622555 -1.1240966  0.6437683  1.58381340 -0.50375974 1.351635678  0.4267090

# redo time columns
dt <- 0.1
time <- seq(0, by = dt, length.out = dim(ChI_NGF_control_data)[1])

ChI_NGF_control_data <- cbind(time = time, ChI_NGF_control_data)
ChI_NGF_6OHDA_data <- cbind(time = time, ChI_NGF_6OHDA_data)



# ChI_NGF_control_data is a list of individual experiments n=1.... with at least Control and GABAzine application then Mecamylamine and NBQX-AP5 within subjects design

ChI_NGF_control_data[1:10,]

#    time    24d12009   24d12013   25304011    25306003    25306008    25409011   25409015
# 1   0.0 -0.48649900  0.6825439  1.1930847 -1.39470208 -0.01547851 1.034252881 -0.7939941
# 2   0.1 -0.36442869 -0.5381592  0.2165222 -0.56462400 -0.06430664 0.008862304 -0.2080566
# 3   0.2  0.17268066 -0.2696045 -0.2717590 -0.49138181 -0.84555660 0.008862304  0.4999511
# 4   0.3 -0.53532712 -0.8311279  0.2470398 -0.05192871 -0.45493162 0.423901347  0.3290527
# 5   0.4 -0.55974118  0.1210205  0.9794616 -0.63786618 -0.77231442 0.423901347 -0.3057129
# 6   0.5 -0.51091306 -0.3916748  1.4677429 -0.46696775 -0.57700193 0.448315408 -1.6729003
# 7   0.6  0.07502441 -0.2207764  1.3151550  0.21662597 -0.21079101 0.375073224 -1.1113769
# 8   0.7  0.49006345  0.1942627  0.9794616  0.33869627 -0.60141599 0.204174795 -1.3066894
# 9   0.8  0.09943847 -0.5381592  0.1860046  0.77814938 -0.67465817 0.985424758  0.3046387
# 10  0.9 -0.70622555 -1.1240966  0.6437683  1.58381340 -0.50375974 1.351635678  0.4267090

# metadata
baseline <- 100
stimulation_time <- 100

x <- ChI_NGF_control_data[,1]
y <- ChI_NGF_control_data[,2]
plot(x, y, type='l')
# Add a vertical dotted line at x = 1350
abline(v = stimulation_time, lty = 2)

filter <- FALSE

peak_fun <- function(data_mat, stimulation_time, baseline,
                     detection_window = 50, dt = 0.1, smooth = 50) {
  
  trace_cols <- setdiff(colnames(data_mat), "time")
  
  sapply(trace_cols, function(col) {
    y <- data_mat[[col]]
    t <- data_mat[["time"]]
    
    keep <- !(is.na(y) | is.na(t))
    y <- y[keep]
    t <- t[keep]
    
    idx <- which(t <= (stimulation_time + detection_window))
    
    if (length(idx) == 0) {
      return(NA)
    }
    
    peak.fun(y = y[idx], dt = dt, stimulation_time = stimulation_time, baseline = baseline, smooth = smooth)
  })
}

charge_fun <- function(data_mat, fc = 300, dt = 0.1, tmax = NULL,
                       baseline = 10, filter = TRUE, showplot = TRUE) {
  
  sapply(2:ncol(data_mat), function(ii) {
    x <- data_mat$time
    y <- data_mat[[ii]]
    
    y <- y[!is.na(y)]
    x <- x[!is.na(y)]
    
    charge_transfer_fun(x, y, tmax = tmax, fc = fc, dt = dt, baseline = baseline, filter = filter, showplot = showplot)
  })
}



ChI_NGF_control_peaks <- peak_fun(ChI_NGF_control_data, stimulation_time=stimulation_time, baseline=baseline, detection_window=250, dt=dt, smooth=50)

ChI_NGF_6OHDA_peaks <- peak_fun(ChI_NGF_6OHDA_data, stimulation_time=stimulation_time, baseline=baseline, detection_window=250, dt=dt, smooth=50)

ChI_NGF_control_areas <- -charge_fun(ChI_NGF_control_data,  baseline=baseline, filter=filter, tmax=400)
graphics.off()

ChI_NGF_6OHDA_areas <- -charge_fun(ChI_NGF_6OHDA_data, baseline=baseline, filter=filter, tmax=400)
graphics.off()

names(ChI_NGF_control_peaks) <- names(ChI_NGF_control_areas)
names(ChI_NGF_6OHDA_peaks) <- names(ChI_NGF_6OHDA_areas)

# List of objects to keep
keep_objects <- c('ChI_NGF_control_data', 'ChI_NGF_control_peaks', 'ChI_NGF_control_areas', 'ChI_NGF_6OHDA_data', 'ChI_NGF_6OHDA_peaks', 'ChI_NGF_6OHDA_areas', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(identifier, '.RData')))







