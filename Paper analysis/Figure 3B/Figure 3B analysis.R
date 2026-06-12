# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 3B'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'NDNF GABA PSCs'

NDNF_GABA_PSCs_data <- load_data2(wd=xlsx_path, name=name)[[1]]
NDNF_GABA_PSCs_data[1:10,]

#       24d02000 24d02006  24d03000   24d03002   24d05002    24d05007    26211001    26211007
# 1   2.55715930 1.481018 -2.659241 -1.9159667  1.5502013  0.12875366  0.10556640  0.99190262
# 2   1.51956170 1.786194 -1.804748 -0.9882324 -0.7080993 -0.23745726 -0.35830076  1.39880364
# 3  -0.18942260 1.877746 -2.079407 -2.1601073 -1.3794860  2.44808948  0.61826169  1.11397293
# 4  -0.03683472 1.328430 -2.354065 -1.0614746 -0.6775818 -1.64126579  1.05771479  0.21879068
# 5   0.20730590 2.152405 -3.269592  0.1592285 -0.8301696 -2.89248643 -0.30947264 -0.39156085
# 6  -0.64718625 1.206360 -3.391662 -0.9394043  1.8248595 -2.46524036 -0.40712889 -0.26949055
# 7  -1.28805536 2.854309 -3.117004 -1.8671386  1.3365783 -3.07559189 -0.08974609  0.05603027
# 8  -0.86080929 3.006897 -2.476135 -1.3056152  1.8553771 -0.75625607 -0.28505858  0.95121252
# 9  -0.28097533 2.549133 -1.163879 -0.9394043  3.1371153  0.52548215  0.44736326  0.46293129
# 10  0.20730590 1.236877 -1.621643 -1.4765136  2.9234923  0.09823608  0.88681636 -0.71708167

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
# analyse_PSC(response=NDNF_GABA_PSCs_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=NDNF_GABA_PSCs_data[1:15000,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=740, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, 
  fast.decay.limit=c(30, 500), latency.limit=25, downsample=downsample) 

out1a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[1:15000,1], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=740, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, 
  latency.limit=25, downsample=downsample) 

out2  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=915, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25, 
  downsample=downsample)  

out2a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=915, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, 
  latency.limit=25, downsample=downsample)  

out3  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,3], dt=dt, n=n, func=func, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1060, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25, 
  fast.decay.limit=c(30, 500), half_width_fit_limit=1000, downsample=10) 

out3a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,3], dt=dt, n=n, func=product1N, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1060, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25, 
  half_width_fit_limit=1000, downsample=10) 

out4  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=950, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
  fast.decay.limit=c(30, 500), downsample=downsample) 

out4a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=950, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
  downsample=downsample) 

out5  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,5], dt=dt, n=10, func=func, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1120, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, 
  latency.limit=25, fast.decay.limit=c(30, 500), half_width_fit_limit = 1000, downsample=10)

out5a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,5], dt=dt, n=10, func=product1N, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1120, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, downsample=10)

out6  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=615, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency=25,  
  downsample=downsample)

out6a  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,6], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=615, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency=25,  
  downsample=downsample)

out7  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,7], dt=dt, n=1000, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=915, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
  fast.decay.limit=c(30, 500), downsample=10) 

out7a <- analyse_PSC(response=NDNF_GABA_PSCs_data[,7], dt=dt, n=1000, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=915, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
   downsample=10) 

out8  <- analyse_PSC(response=NDNF_GABA_PSCs_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=705, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
  downsample=downsample) 

out8a <- analyse_PSC(response=NDNF_GABA_PSCs_data[,8], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=705, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25,
  downsample=downsample) 

out1$BIC < out1a$BIC 
out2$BIC < out2a$BIC 
out3$BIC < out3a$BIC 
out4$BIC < out4a$BIC 
out5$BIC < out5a$BIC 
out6$BIC < out6a$BIC 
out7$BIC < out7a$BIC 
out8$BIC < out8a$BIC 

# organise outputs
NDNF_GABA_PSCs_summary <- list(out1a, out2a, out3, out4, out5, out6, out7, out8)

names(NDNF_GABA_PSCs_summary) <- 1:length(NDNF_GABA_PSCs_summary)

NDNF_GABA_PSCs_fits <- t(sapply(1:length(NDNF_GABA_PSCs_summary), function(ii){
  X <- NDNF_GABA_PSCs_summary[[ii]]$output
  X <- if (dim(X)[1] == 1) c(rep(NA, dim(X)[2]), as.vector(t(X))) else as.vector(t(X))
  as.vector(t(X))
  })
)
# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(NDNF_GABA_PSCs_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(NDNF_GABA_PSCs_fits) <- new_colnames
rownames(NDNF_GABA_PSCs_fits) <- 1:length(NDNF_GABA_PSCs_summary)


NDNF_GABA_PSCs_peaks <- sapply(1:length(NDNF_GABA_PSCs_data), function(ii){
  y <- NDNF_GABA_PSCs_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y[x<2000], dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


NDNF_GABA_PSCs_areas <- sapply(1:length(NDNF_GABA_PSCs_data), function(ii){
  y <- NDNF_GABA_PSCs_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y[x<2000], x=x[x<2000], dt=dt, baseline=baseline)
  })


NDNF_GABA_PSCs_fits[, "A1"][is.na(NDNF_GABA_PSCs_fits[, "A1"])] <- 0
NDNF_GABA_PSCs_fits[, "area1"][is.na(NDNF_GABA_PSCs_fits[, "area1"])] <- 0
NDNF_GABA_PSCs_fits
# using n <- 100 and method <- 'BF.LM':
#         A1 τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area1
# 1    0.000    NA     NA     NA     NA     NA     NA         NA     0.000  -83.338 27.970 220.487  66.139 23.931 308.240 22.490    235.288  24802.818
# 2    0.000    NA     NA     NA     NA     NA     NA         NA     0.000 -132.380 16.370 294.826  50.107 16.809 408.856 20.454    262.131  46259.139
# 3 -199.125 7.096 30.000 13.400  5.052 43.382 16.204     39.047  9337.491 -350.525 41.871 384.282 104.169 37.235 535.528 12.749    394.165 176642.619
# 4 -171.350 2.722 26.226  6.880  2.449 36.520 15.170     26.559  5841.865 -241.536 39.880 302.065  93.029 33.763 422.730 17.000    325.936  99273.913
# 5 -536.365 1.412 30.000  4.529  1.492 41.595 14.628     25.960 18712.888 -696.715 31.644 425.388  88.833 30.688 590.484 12.905    399.624 365202.220
# 6  -23.617 8.849  8.868  8.858  3.472 19.221 20.328     21.670   568.676  -47.292 20.882 135.323  46.145 16.940 190.369 16.019    152.810   9000.177
# 7  -81.882 7.095  8.363  7.694  3.014 16.786 19.314     18.863  1718.308 -162.242 30.230 313.420  78.245 27.680 436.008 18.000    311.818  65269.794
# 8 -153.808 5.357  5.363  5.360  2.101 11.631 19.129     13.113  2241.004 -146.833 26.691 259.408  67.659 24.069 361.181 16.042    262.146  49440.441

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
data <- cbind(-NDNF_GABA_PSCs_fits[,1], -NDNF_GABA_PSCs_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps = 100
setwd(svg_path)
DBSCAN_analyse(data, eps=400, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('NDNF_GABA_PSCs_peaks', 'NDNF_GABA_PSCs_areas', 'NDNF_GABA_PSCs_fits', 
  'NDNF_GABA_PSCs_summary', 'NDNF_GABA_PSCs_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, '.RData')))

# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 3B'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

name <- 'NDNF GABA PSCs'

rdata_path <- file.path(analysis_path, paste0(name, '.RData'))
if (!file.exists(rdata_path)) {
  stop('RData file not found at: ', rdata_path)
}

load(rdata_path)

NDNF_GABA_PSCs_fits

#         A1 τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area1
# 1    0.000    NA     NA     NA     NA     NA     NA         NA     0.000  -83.338 27.970 220.487  66.139 23.931 308.240 22.490    235.288  24802.818
# 2    0.000    NA     NA     NA     NA     NA     NA         NA     0.000 -132.380 16.370 294.826  50.107 16.809 408.856 20.454    262.131  46259.139
# 3 -199.125 7.096 30.000 13.400  5.052 43.382 16.204     39.047  9337.491 -350.525 41.871 384.282 104.169 37.235 535.528 12.749    394.165 176642.619
# 4 -171.350 2.722 26.226  6.880  2.449 36.520 15.170     26.559  5841.865 -241.536 39.880 302.065  93.029 33.763 422.730 17.000    325.936  99273.913
# 5 -536.365 1.412 30.000  4.529  1.492 41.595 14.628     25.960 18712.888 -696.715 31.644 425.388  88.833 30.688 590.484 12.905    399.624 365202.220
# 6  -23.617 8.850  8.866  8.858  3.472 19.221 20.328     21.670   568.675  -47.292 20.882 135.323  46.145 16.939 190.369 16.019    152.810   9000.183
# 7  -81.882 7.095  8.363  7.694  3.014 16.786 19.314     18.863  1718.308 -162.242 30.230 313.420  78.245 27.680 436.008 18.000    311.818  65269.794
# 8 -153.808 5.357  5.363  5.360  2.101 11.631 19.129     13.113  2241.005 -146.833 26.691 259.408  67.659 24.069 361.181 16.042    262.146  49440.443

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
  for (ii in 1:length(NDNF_GABA_PSCs_summary)){
  traces <- NDNF_GABA_PSCs_summary[[ii]]$traces
  func <- if (dim(traces)[2]==4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


