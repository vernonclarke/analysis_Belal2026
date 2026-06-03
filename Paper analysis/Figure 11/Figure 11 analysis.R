# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 11'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

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

ChI_NGF_6OHDA_data[1:10,]
#      25310000    25310001   25310005    25310006    25310007   25310012   25311000    25311002  25311007    25312000     25312010
# 1  -1.9458739 -0.79577633 -0.3421936 -1.59694817 -0.14765624  0.3149414  2.2543212 -1.33498529 -3.186767  0.90334468  1.496166921
# 2  -2.2388427 -0.42956541  0.3291931 -2.01198721  0.46269529  0.6323242  0.2767822 -0.26076659 -2.271240  0.02443848  1.129956001
# 3  -1.3599365 -0.73474118 -0.7694397 -1.40163568  0.36503905  0.2172851 -0.4800537 -0.13869628 -2.393310 -0.24411620  0.226635731
# 4  -1.5796630 -0.49060056  0.2376404 -0.30300292  0.07207031 -0.1245117 -0.2603271  0.03220215 -3.125732  0.19533690 -0.261645495
# 5  -1.2622802  0.18078612  0.5733337 -0.03444824  0.56035154 -0.4663086  0.5941650 -0.48049314 -2.973144  0.85451656 -0.261645495
# 6  -1.0913818 -0.55163572  0.1766052 -0.25417479  0.48710935 -0.4907226  0.4720947  0.12985839 -1.966064  0.75686032  0.006909179
# 7  -1.0425537  0.27233885  1.2447204 -1.01101069 -0.17207030  0.2661133 -1.1392333  1.03317866 -2.362793  0.68361813 -0.066333005
# 8  -0.1392334 -0.39904783  0.6648864 -1.59694817  0.21855468  0.6323242  0.2035400 -0.04104004 -1.752441  0.31740721  0.788159142
# 9  -0.4810303 -0.49060056  0.3597107 -0.62038571  0.46269529  0.5590820  0.1547119 -0.38283690 -1.538818 -0.48825681  2.130932516
# 10 -1.2622802  0.05871582 -0.8915100  0.16086425  0.12089843 -0.1977539 -0.3335693 -0.13869628 -1.660889 -0.87888179  1.984448148

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
keep_objects <- c('ChI_NGF_control_data', 'ChI_NGF_control_peaks', 'ChI_NGF_control_areas', 'ChI_NGF_6OHDA_data', 
  'ChI_NGF_6OHDA_peaks', 'ChI_NGF_6OHDA_areas', 'analysis_path', 'identifier')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(identifier, '.RData')))







