# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 1'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
baseline <- 100
stimulation_time <- 100
dt <- 0.1

setwd(analysis_path)

# Load data
dSPN_data <- load_data2(wd=xlsx_path, name='ChAT-Cre X tdTomato', header=TRUE)
names(dSPN_data) <- seq_along(dSPN_data)

iSPN_data <- load_data2(wd=xlsx_path, name='ChAT-Cre X De eGFP', header=TRUE)
names(iSPN_data) <- seq_along(iSPN_data)

# Redo time columns
fix_time <- function(df) {
  df$time <- seq(0, by = dt, length.out = nrow(df))
  df[, c("time", setdiff(names(df), "time"))]
}
dSPN_data <- lapply(dSPN_data, fix_time)
iSPN_data <- lapply(iSPN_data, fix_time)


# dSPN_data is a list of individual experiments n=1.... with at least Control and GABAzine application then Mecamylamine and NBQX-AP5 within subjects design

dSPN_data[[1]][1:10,]

#    time     control  GABAzine Mecamylamine
# 1   0.0  0.95540360 -4.685669   0.21138508
# 2   0.1  0.34505207 -7.127075   0.41483559
# 3   0.2  0.54850258 -8.551228   1.22863764
# 4   0.3  0.54850258 -8.754679   0.82173661
# 5   0.4 -0.06184896 -7.940877   0.00793457
# 6   0.5 -1.07910151 -8.347778   1.63553866
# 7   0.6 -0.67220049 -8.958129   2.44934070
# 8   0.7 -0.46874998 -8.347778   1.63553866
# 9   0.8 -0.26529947 -7.533976   2.04243968
# 10  0.9 -1.68945304 -7.940877   1.63553866

# Peak detection function
peak_fun <- function(data_list, condition, stimulation_time, baseline, detection_window=50, dt=0.1, smooth=50) {
  sapply(data_list, function(df) {
    if (condition %in% colnames(df)) {
      y <- df[[condition]]
      y <- y[!is.na(y)]
      idx <- (stimulation_time + detection_window) / dt
      peak.fun(y=y[1:idx], dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
    } else {
      NA
    }
  })
}

# Calculate peaks
conditions <- c('control', 'GABAzine', 'Mecamylamine')
peak_args <- list(stimulation_time=stimulation_time, baseline=baseline, detection_window=250, dt=dt, smooth=50)

dSPN_peaks <- do.call(cbind, lapply(conditions, function(cond) 
  peak_fun(data_list=dSPN_data, condition=cond, stimulation_time=peak_args$stimulation_time, 
           baseline=peak_args$baseline, detection_window=peak_args$detection_window, 
           dt=peak_args$dt, smooth=peak_args$smooth)))
colnames(dSPN_peaks) <- conditions

iSPN_peaks <- do.call(cbind, lapply(conditions, function(cond) 
  peak_fun(data_list=iSPN_data, condition=cond, stimulation_time=peak_args$stimulation_time, 
           baseline=peak_args$baseline, detection_window=peak_args$detection_window, 
           dt=peak_args$dt, smooth=peak_args$smooth)))
colnames(iSPN_peaks) <- conditions

# Calculate areas
dSPN_areas <- do.call(cbind, lapply(conditions, function(cond) 
  -charge_fun(data_list=dSPN_data, condition=cond, baseline=baseline, filter=FALSE)))
colnames(dSPN_areas) <- conditions

iSPN_areas <- do.call(cbind, lapply(conditions, function(cond) 
  -charge_fun(data_list=iSPN_data, condition=cond, baseline=baseline, filter=FALSE)))
colnames(iSPN_areas) <- conditions

graphics.off()

# Clean up and save
keep_objects <- c('dSPN_data', 'dSPN_peaks', 'dSPN_areas', 'iSPN_data', 'iSPN_peaks', 'iSPN_areas', 'identifier', 'analysis_path')
rm(list = setdiff(ls(), keep_objects))
save.image(file = file.path(analysis_path, paste0(identifier, '.RData')))

dSPN_peaks
#       control    GABAzine Mecamylamine
# 1  -207.33005  -99.449725    -77.65431
# 2  -167.38608  -91.693132    -43.19042
# 3  -235.58272 -104.946136    -47.99995
# 4  -188.39128 -101.077548           NA
# 5  -246.94094 -165.512406    -58.54199
# 6  -102.27480   -7.535083           NA
# 7  -179.88431  -30.131268           NA
# 8  -180.82227 -157.741041           NA
# 9   -31.05392  -23.394194    -22.26879
# 10 -323.84937  -73.901788    -51.60650

iSPN_peaks
#      control   GABAzine Mecamylamine
# 1 -238.67171  -53.76314           NA
# 2 -213.46519  -33.49574           NA
# 3  -56.75674  -16.42002           NA
# 4 -312.67521 -124.20345    -73.16225
# 5 -217.75717 -134.98472    -61.00537
# 6 -234.50310  -39.86229           NA

dSPN_areas
#      control   GABAzine Mecamylamine
# 1  22.004771  5.2686697     3.735078
# 2  24.344662  4.6795131     1.585975
# 3  52.980218  4.4188235     3.210763
# 4  15.012641  5.4502999           NA
# 5  29.185738  4.8532703     1.828901
# 6   8.496155  0.4595433           NA
# 7  12.042255  1.3983342           NA
# 8  27.356642 19.4382113           NA
# 9   6.068458  4.1780631     4.282606
# 10 65.220723  8.7389735     6.516338

iSPN_areas
#     control GABAzine Mecamylamine
# 1 55.552952 2.758048           NA
# 2 49.880087 2.516400           NA
# 3  6.170359 1.921617           NA
# 4 21.757989 5.008854     2.175596
# 5 48.404705 9.696649    10.940780
# 6 68.098898 5.110898           NA

