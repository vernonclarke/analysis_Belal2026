# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

# Settings
identifier <- 'Figure 10'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)
name <- 'NPY Cre X iSPN 6OHDA' # iSPN data
NPY_Cre_X_iSPN_6OHDA_data <- load_data2(wd=xlsx_path, name=name)[[1]]

NPY_Cre_X_iSPN_6OHDA_data[1:10,]
#       23412000   23412004    23412005   23412007   23412010    23413000    23413001  23413005   23n02013
# 1   0.40589598  0.4155477  0.40637205  1.0885772  1.8497314  1.38387038 -0.42947386 1.2524617 -1.6004638
# 2  -0.03355713  0.4969279  0.04016113  0.5392608  1.2393798  0.79386389  0.30294798 1.3541869 -3.5047606
# 3  -0.36314696  0.9038289 -0.32604979  0.4782257  0.7816162  0.06144205  0.57760617 1.6593627 -3.7489012
# 4   0.18616942  1.0869344  0.38602700  0.5697784  0.7205810 -0.65063473  0.54708860 1.5372924 -1.4295654
# 5   0.57679441  1.5141804  0.12154134  0.5697784  1.1173095 -0.91512040  0.18087768 1.2931518 -0.9168701
# 6   0.69886471  1.5141804 -0.40742999  1.2259063  0.7510986 -0.67097979 -0.42947386 1.7814330 -2.7235106
# 7   0.57679441  1.3921101 -0.38708494  0.9359893 -0.2864990  0.91593420 -0.33792113 0.7031453 -3.4071043
# 8   0.40589598  0.8224487 -0.69226071  0.3561554  0.2017822  2.58422839 -0.58206174 0.2352091 -2.1375731
# 9   0.49134519 -0.2354940 -0.30570474 -0.1779022  0.3848877  1.91284171  0.05880737 0.6421102 -2.1375731
# 10  0.58900144 -1.3544718 -0.36673989 -0.2847137 -0.2254639  1.89249665 -0.36843870 0.6014201 -2.3084716

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
# analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=860, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, half_width_fit_limit=1000, fast.decay.limit=c(100,500), latency.limit=20)

#           A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# fast -45.758  3.080  91.172  10.799  3.427 126.393 4.577     75.291  4696.459
# slow -28.192 63.403 326.170 128.905 48.046 464.351 1.962    405.532 13652.173

out2  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#            A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width     area1
# fast -183.839 1.711  10.275 3.681  1.358  14.502 5.179     11.879  2702.604
# slow  -96.653 1.130 118.579 5.309  1.429 164.385 3.787     87.898 11985.771

out3  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=475, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast  -9.302  6.473   6.485  6.479  2.540  14.059 4.383     15.850  163.830
# slow -37.810 13.205 139.381 34.376 12.142 193.856 6.032    138.085 6744.163

out4  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -31.284 5.492   5.502  5.497  2.155  11.928 4.981     13.447  467.451
# slow -49.662 3.613 151.573 13.829  4.205 210.125 3.365    120.321 8246.571

out5  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=325, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -159.597  5.225   5.242  5.234  2.052  11.357  2.746     12.804 2270.651
# slow  -56.981 10.220 138.260 28.746  9.924 191.913 13.640    129.718 9698.934

out6  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=280, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20) 

#            A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width     area1
# fast -234.010 1.833   8.699 3.617  1.355  12.456 4.161     10.869  3085.194
# slow -113.134 1.137 135.556 5.483  1.450 187.920 5.550     99.837 15968.887

out7  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=475, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -42.266 1.971   9.667  3.936  1.471  13.807 4.337     11.941  613.927
# slow -54.131 6.116 158.054 20.689  6.670 219.119 3.464    132.890 9752.172

out8  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(100,500), latency.limit=20) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -440.901  2.841  76.153  9.705  3.115 105.574 4.635     63.710 38139.38
# slow -150.636 31.571 189.214 67.854 25.043 267.087 6.800    218.872 40796.85

out9  <- analyse_PSC(response=NPY_Cre_X_iSPN_6OHDA_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=750, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, half_width_fit_limit=1000, latency.limit=20)

#           A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width    area1
# fast -47.117  7.552 138.952  23.257  7.784 192.688 7.207    123.091 7739.849
# slow -13.979 42.129 522.430 115.374 40.162 725.550 0.298    502.572 9107.897

# organise outputs
NPY_Cre_X_iSPN_6OHDA_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8, out9)

names(NPY_Cre_X_iSPN_6OHDA_summary) <- 1:length(NPY_Cre_X_iSPN_6OHDA_summary)

NPY_Cre_X_iSPN_6OHDA_fits <- t(sapply(1:length(NPY_Cre_X_iSPN_6OHDA_summary), function(ii){
  X <- NPY_Cre_X_iSPN_6OHDA_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(NPY_Cre_X_iSPN_6OHDA_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(NPY_Cre_X_iSPN_6OHDA_fits) <- new_colnames
rownames(NPY_Cre_X_iSPN_6OHDA_fits) <- 1:length(NPY_Cre_X_iSPN_6OHDA_summary)

NPY_Cre_X_iSPN_6OHDA_peaks <- sapply(1:length(NPY_Cre_X_iSPN_6OHDA_data), function(ii){
  y <- NPY_Cre_X_iSPN_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

NPY_Cre_X_iSPN_6OHDA_areas <- sapply(1:length(NPY_Cre_X_iSPN_6OHDA_data), function(ii){
  y <- NPY_Cre_X_iSPN_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]

  tmax <- 1000
  charge_transfer_fun(y=y[x<tmax], x=x[x<tmax], dt=dt, baseline=baseline)
  })


NPY_Cre_X_iSPN_6OHDA_fits
# using n <- 100 and method <- 'BF.LM':

#         A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -45.758 3.080  91.172 10.799  3.427 126.393 4.577     75.291  4696.459  -28.192 63.403 326.170 128.905 48.046 464.351  1.962    405.532 13652.173
# 2 -183.839 1.711  10.275  3.681  1.358  14.502 5.179     11.879  2702.604  -96.653  1.130 118.579   5.309  1.429 164.385  3.787     87.898 11985.771
# 3   -6.016 6.164   6.171  6.168  2.418  13.383 6.361     15.088   100.854  -38.008 13.531 139.373  34.949 12.371 193.902  4.385    138.876  6807.087
# 4  -31.284 5.492   5.502  5.497  2.155  11.928 4.981     13.447   467.453  -49.662  3.613 151.573  13.829  4.205 210.125  3.365    120.321  8246.571
# 5 -159.597 5.225   5.242  5.234  2.052  11.357 2.746     12.804  2270.651  -56.981 10.220 138.260  28.746  9.924 191.913 13.640    129.718  9698.934
# 6 -234.010 1.833   8.699  3.617  1.355  12.456 4.161     10.869  3085.194 -113.134  1.137 135.556   5.483  1.450 187.920  5.550     99.837 15968.887
# 7  -42.266 1.971   9.667  3.936  1.471  13.807 4.337     11.941   613.927  -54.131  6.116 158.054  20.689  6.670 219.119  3.464    132.890  9752.172
# 8 -440.901 2.841  76.153  9.705  3.115 105.574 4.635     63.710 38139.375 -150.636 31.571 189.214  67.854 25.043 267.087  6.800    218.872 40796.846
# 9  -47.117 7.552 138.952 23.257  7.784 192.688 7.207    123.091  7739.849  -13.979 42.129 522.430 115.374 40.162 725.550  0.298    502.572  9107.897

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-NPY_Cre_X_iSPN_6OHDA_fits[,1], -NPY_Cre_X_iSPN_6OHDA_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=100
setwd(svg_path)
DBSCAN_analyse(data=data, eps=80, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('NPY_Cre_X_iSPN_6OHDA_peaks', 'NPY_Cre_X_iSPN_6OHDA_areas', 'NPY_Cre_X_iSPN_6OHDA_fits', 
  'NPY_Cre_X_iSPN_6OHDA_summary', 'NPY_Cre_X_iSPN_6OHDA_data', 'name', 'analysis_path')

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

identifier <- 'Figure 10'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

name <- 'NPY Cre X iSPN 6OHDA' # iSPN data

load(paste0(analysis_path, '/', name, '.RData'))

NPY_Cre_X_iSPN_6OHDA_fits        A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -45.758 3.080  91.172 10.799  3.427 126.393 4.577     75.291  4696.459  -28.192 63.403 326.170 128.905 48.046 464.351  1.962    405.532 13652.173
# 2 -183.839 1.711  10.275  3.681  1.358  14.502 5.179     11.879  2702.604  -96.653  1.130 118.579   5.309  1.429 164.385  3.787     87.898 11985.771
# 3   -6.016 6.164   6.171  6.168  2.418  13.383 6.361     15.088   100.854  -38.008 13.531 139.373  34.949 12.371 193.902  4.385    138.876  6807.087
# 4  -31.284 5.492   5.502  5.497  2.155  11.928 4.981     13.447   467.453  -49.662  3.613 151.573  13.829  4.205 210.125  3.365    120.321  8246.571
# 5 -159.597 5.225   5.242  5.234  2.052  11.357 2.746     12.804  2270.651  -56.981 10.220 138.260  28.746  9.924 191.913 13.640    129.718  9698.934
# 6 -234.010 1.833   8.699  3.617  1.355  12.456 4.161     10.869  3085.194 -113.134  1.137 135.556   5.483  1.450 187.920  5.550     99.837 15968.887
# 7  -42.266 1.971   9.667  3.936  1.471  13.807 4.337     11.941   613.927  -54.131  6.116 158.054  20.689  6.670 219.119  3.464    132.890  9752.172
# 8 -440.901 2.841  76.153  9.705  3.115 105.574 4.635     63.710 38139.375 -150.636 31.571 189.214  67.854 25.043 267.087  6.800    218.872 40796.846
# 9  -47.117 7.552 138.952 23.257  7.784 192.688 7.207    123.091  7739.849  -13.979 42.129 522.430 115.374 40.162 725.550  0.298    502.572  9107.897

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
for (ii in 1:length(NPY_Cre_X_iSPN_6OHDA_summary)){
  fit_plot(traces=NPY_Cre_X_iSPN_6OHDA_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


