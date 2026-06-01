# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

identifier <- 'Figure 3'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'CRISPR control'

CRISPR_control_data <- load_data2(wd=xlsx_path, name=name)[[1]]

CRISPR_control_data[1:10,]
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

# ==============================================
# METADATA FOR ANAYSIS
# ==============================================
n <- 100
stimulation_time <- 100
baseline <- 100
dt <- 0.1
method <- 'BF.LM'
weight_method <- 'none'
downsample <- 1
interval=c(0.2, 0.8) # % rise and decay
# MLE as initial estimate then MCMC and the Metropolis–Hastings algorithm to obtain posterior
MLEsettings <- list(iter=1000, metropolis.scale=1.5, fit.attempts=10, RWm=FALSE)
func <- product2N
# limits 5 xs previous estimates (see English et al., 2012)
fast.decay.limit <- c(30, 500)
fast.decay.limit <- NULL
# fitting to 10% of the peak to avoid any slow components (see English et al., 2012)
rel.decay.fit.limit <- 0.1

# ==============================================
# FITTING DATA
# ==============================================

# to analyse any trace:
# analyse_PSC(response=CRISPR_control_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=CRISPR_control_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=385, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -141.120  1.566   5.357  2.722  1.037   7.939 10.140      7.549  1256.568
# slow -211.089 13.003 100.599 30.553 11.071 140.705  8.138    107.909 28771.237

out2  <- analyse_PSC(response=CRISPR_control_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=365, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -78.429  1.852   5.747  3.094  1.184   8.641 10.657      8.414   772.207
# slow -99.447 15.454 100.255 34.163 12.540 141.031  5.007    113.173 14018.046

out3  <- analyse_PSC(response=CRISPR_control_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=365, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500)) 

#            A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast -143.396  2.038  4.137  2.844  1.104   6.821 8.973      7.245  1179.683
# slow -298.921 12.815 83.320 28.351 10.405 117.197 7.356     93.993 35000.992

out4  <- analyse_PSC(response=CRISPR_control_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=285, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -208.987  1.297  4.164  2.198  0.840   6.228 11.636      6.016  1475.159
# slow -241.229 14.317 70.111 28.581 10.685 100.151  4.613     86.655 25424.875

out5  <- analyse_PSC(response=CRISPR_control_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1 τrise τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast  -74.823 1.993  1.994  1.993  0.781   4.326 8.224      4.877   405.435
# slow -143.314 9.403 85.980 23.366  8.354 119.827 2.946     88.223 16169.945

out6  <- analyse_PSC(response=CRISPR_control_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=475, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -296.552 22.099 22.138 22.119  8.670  47.995  7.403     54.111 17830.02
# slow -257.915 78.690 78.788 78.739 30.865 170.858 23.757    192.626 55202.95

out7  <- analyse_PSC(response=CRISPR_control_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=495, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast  -99.433  2.688   2.690  2.689  1.054   5.834 11.06      6.578   726.744
# slow -176.352 11.752 105.063 28.985 10.383 146.481  7.18    108.425 24414.292

# organise outputs
CRISPR_control_summary <- list(out1, out2, out3, out4, out5, out6, out7)

names(CRISPR_control_summary) <- 1:length(CRISPR_control_summary)

CRISPR_control_fits <- t(sapply(1:length(CRISPR_control_summary), function(ii){
  X <- CRISPR_control_summary[[ii]]$output
  X <- if (dim(X)[1] == 1) c(rep(NA, dim(X)[2]), as.vector(t(X))) else as.vector(t(X))
  as.vector(t(X))
  })
)
# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(CRISPR_control_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(CRISPR_control_fits) <- new_colnames
rownames(CRISPR_control_fits) <- 1:length(CRISPR_control_summary)


CRISPR_control_peaks <- sapply(1:length(CRISPR_control_data), function(ii){
  y <- CRISPR_control_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y[x<900], dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


CRISPR_control_areas <- sapply(1:length(CRISPR_control_data), function(ii){
  y <- CRISPR_control_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, baseline=baseline, tmax=900)
  })


CRISPR_control_fits[, "A1"][is.na(CRISPR_control_fits[, "A1"])] <- 0
CRISPR_control_fits[, "area1"][is.na(CRISPR_control_fits[, "area1"])] <- 0
CRISPR_control_fits
# using n <- 100 and method <- 'BF.LM':
#         A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -141.120  1.566  5.357  2.722  1.037  7.939 10.140      7.549  1256.568 -211.089 13.003 100.599 30.553 11.071 140.705  8.138    107.909 28771.24
# 2  -78.429  1.852  5.747  3.094  1.184  8.641 10.657      8.414   772.207  -99.447 15.454 100.255 34.163 12.540 141.031  5.007    113.173 14018.05
# 3 -143.396  2.038  4.137  2.844  1.104  6.821  8.973      7.245  1179.683 -298.921 12.815  83.320 28.351 10.405 117.197  7.356     93.993 35000.99
# 4 -208.987  1.297  4.164  2.198  0.840  6.228 11.636      6.016  1475.159 -241.229 14.317  70.111 28.581 10.685 100.151  4.613     86.655 25424.88
# 5  -74.823  1.993  1.994  1.993  0.781  4.325  8.224      4.877   405.432 -143.314  9.403  85.980 23.366  8.354 119.827  2.946     88.223 16169.95
# 6 -296.552 22.098 22.139 22.119  8.670 47.995  7.403     54.111 17830.022 -257.915 78.696  78.783 78.739 30.865 170.858 23.757    192.626 55202.95
# 7  -99.433  2.688  2.690  2.689  1.054  5.834 11.060      6.578   726.744 -176.352 11.752 105.063 28.985 10.383 146.481  7.180    108.425 24414.29

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
data <- cbind(-CRISPR_control_fits[,1], -CRISPR_control_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps = 150
setwd(svg_path)
DBSCAN_analyse(data, eps=150, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('CRISPR_control_peaks', 'CRISPR_control_areas', 'CRISPR_control_fits', 
  'CRISPR_control_summary', 'CRISPR_control_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, ".RData")))

# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

# Paths
repo_root <- normalizePath(
  '~/Documents/Repositories/analysis_Belal2026',
  mustWork = TRUE
)

identifier <- 'Figure 3'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

name <- 'CRISPR control'

load(paste0(analysis_path, '/', name, '.RData'))

CRISPR_control_fits

#         A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -141.120  1.566  5.357  2.722  1.037  7.939 10.140      7.549  1256.568 -211.089 13.003 100.599 30.553 11.071 140.705  8.138    107.909 28771.24
# 2  -78.429  1.852  5.747  3.094  1.184  8.641 10.657      8.414   772.207  -99.447 15.454 100.255 34.163 12.540 141.031  5.007    113.173 14018.05
# 3 -143.396  2.038  4.137  2.844  1.104  6.821  8.973      7.245  1179.683 -298.921 12.815  83.320 28.351 10.405 117.197  7.356     93.993 35000.99
# 4 -208.987  1.297  4.164  2.198  0.840  6.228 11.636      6.016  1475.159 -241.229 14.317  70.111 28.581 10.685 100.151  4.613     86.655 25424.88
# 5  -74.823  1.993  1.994  1.993  0.781  4.326  8.224      4.877   405.435 -143.314  9.403  85.980 23.366  8.354 119.827  2.946     88.223 16169.94
# 6 -296.552 22.099 22.138 22.119  8.670 47.995  7.403     54.111 17830.019 -257.915 78.690  78.788 78.739 30.865 170.858 23.757    192.626 55202.95
# 7  -99.433  2.688  2.690  2.689  1.054  5.834 11.060      6.578   726.744 -176.352 11.752 105.063 28.985 10.383 146.481  7.180    108.425 24414.29

# save
setwd(analysis_path)

# path where all graphs are stored
svg_path <- paste0(analysis_path, '/svg')
xlsx_path <- paste0(analysis_path, '/xlsx')

if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}
if (!dir.exists(xlsx_path)) {
  dir.create(xlsx_path, recursive = TRUE)
}

save <- TRUE
setwd(svg_path)
  for (ii in 1:length(CRISPR_control_summary)){
  traces <- CRISPR_control_summary[[ii]]$traces
  func <- if (dim(traces)[2]==4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


