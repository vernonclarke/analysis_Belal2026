# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

identifier <- 'Figure 2'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)
name <- 'ChAT-Cre X tdTomato'

ChAT_Cre_X_tdTomato_data <- load_data2(wd=xlsx_path, name=name)

# pad matrix as traces are different lengths
max_len <- max(sapply(ChAT_Cre_X_tdTomato_data, nrow))
ChAT_Cre_X_tdTomato_data <- do.call(cbind, setNames(lapply(ChAT_Cre_X_tdTomato_data, function(x) {
  c(x[,1], rep(NA, max_len - nrow(x)))
}), names(ChAT_Cre_X_tdTomato_data)))

ChAT_Cre_X_tdTomato_data[1:10,]

#       22405005    22405015   22406005   22406013  22411000   22411004  22412005    24o21006    24o21010   24o22003
#  [1,] 5.641072  0.94885249  1.5003356 -1.8062743 -4.265137  1.9136352  8.364705 -0.40327960 -1.55700676 -0.5395914
#  [2,] 7.472127 -2.26566558  1.9886168 -1.9283446 -3.135986  1.6084594  7.632283 -0.52534991 -0.78389482 -2.1265054
#  [3,] 9.099731 -1.24841303  1.7749938 -1.7248941 -4.661865  0.6318969  8.608846 -0.52534991  0.23335774 -2.2078856
#  [4,] 9.303182  0.54195147  0.7679138 -1.1145426 -3.074951 -0.8939819  8.039184 -1.66467277  1.16923009 -3.3065184
#  [5,] 7.879028  0.94885249  2.1106872 -1.2366129 -3.013916 -1.4432983  7.795044 -0.89156083 -0.49906410 -3.4285887
#  [6,] 7.268676 -0.71944170  1.9886168 -0.9517822 -3.532715  0.2656860  9.951619 -0.68811032 -1.55700676 -2.1671955
#  [7,] 8.285929  0.05367024  0.7679138 -0.7076416 -2.525635 -0.7108764 10.114379 -0.03706868 -1.15010574 -2.3706460
#  [8,] 7.879028 -1.61462395 -1.2157287 -0.4228109 -3.105469  0.1436157 10.114379  0.41052244 -0.53975421 -0.3361409
#  [9,] 7.268676 -1.24841303 -0.6053772  0.1061605 -4.295654  0.5098266  9.666788  0.16638183 -0.09216308 -1.2313232
# [10,] 6.251424 -2.67256660  0.4017029  0.8385823 -4.143066 -0.2836304  8.324015 -0.32189940 -0.37699380 -2.6147867

# ==============================================
# METADATA FOR ANAYSIS
# ==============================================

n <- 100
stimulation_time <- 100
baseline <- 100
dt <- 0.1
method <- 'BF.LM'
weight_method <- 'none'
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
# analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)  

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast  -88.305  3.935  3.942  3.938  1.544   8.546 11.360      9.635   945.38
# slow -117.851 11.549 97.056 27.905 10.046 135.483 13.294    101.797 15248.34

out2 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,2], dt=dt, n=10, func=func, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=665, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit = c(30, 500))  

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -88.606 6.407   7.100  6.742  2.642  14.660 25.720     16.507  1625.942
# slow -70.643 7.041 259.004 26.094  8.065 359.056 31.497    208.470 20236.213

out3 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=570, stimulation_time=stimulation_time, baseline=baseline) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -79.712  1.442  31.399  4.655   1.53  43.534 14.273     27.066  2902.916
# slow -210.653 21.857 157.532 50.125  18.26 220.782 27.329    172.395 45616.726

out4 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=405, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit = c(30, 500)) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -69.261 2.956  13.383  5.730  2.152  19.235 10.288     16.999 1422.304
# slow -41.612 3.047 186.091 12.738  3.689 257.977 11.467    142.867 8292.307

out5 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=395, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -132.419  2.720  2.723  2.721  1.067   5.905 14.385      6.658   979.584
# slow -177.149 10.615 90.084 25.732  9.257 125.724 16.978     94.239 21234.360

out6 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=325, stimulation_time=stimulation_time, baseline=baseline) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -51.147 6.279   6.280  6.279  2.461  13.626 22.778     15.362  873.039
# slow -45.042 3.075 117.152 11.494  3.537 162.408 18.579     93.934 5820.659

out7 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline) 

#            A1 τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -101.555 1.565  7.671  3.125  1.168  10.958 10.127      9.478 1170.867
# slow  -86.024 8.587 75.015 21.018  7.543 104.633  7.961     77.880 8539.756

out8 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=960, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25)

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -23.844 7.794   8.711  8.236  3.228  17.915 20.381     20.168  534.621
# slow -10.815 6.774 595.574 30.671  8.458 825.641 18.800    445.916 6781.617

out9 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=690, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit = c(30, 500))

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast  -3.217  0.778   2.969  1.412  0.535   4.340 24.680      4.013   15.368
# slow -10.282 18.944 230.012 51.541 17.977 319.489 10.651    220.773 2958.884

out10 <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,10], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=485, stimulation_time=stimulation_time, baseline=baseline)

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -127.557  3.980   3.991  3.986  1.562   8.648 15.773      9.750  1381.94
# slow -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.69


# organise outputs
ChAT_Cre_X_tdTomato_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8, out9, out10)

names(ChAT_Cre_X_tdTomato_summary) <- 1:length(ChAT_Cre_X_tdTomato_summary)

ChAT_Cre_X_tdTomato_fits <- t(sapply(1:length(ChAT_Cre_X_tdTomato_summary), function(ii){
  X <- ChAT_Cre_X_tdTomato_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(ChAT_Cre_X_tdTomato_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(ChAT_Cre_X_tdTomato_fits) <- new_colnames
rownames(ChAT_Cre_X_tdTomato_fits) <- 1:length(ChAT_Cre_X_tdTomato_summary)

ChAT_Cre_X_tdTomato_peaks <- sapply(1:dim(ChAT_Cre_X_tdTomato_data)[2], function(ii){
  y <- ChAT_Cre_X_tdTomato_data[,ii]
  x = seq(y)*dt - dt 
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

ChAT_Cre_X_tdTomato_areas <- sapply(1:dim(ChAT_Cre_X_tdTomato_data)[2], function(ii){
  y <- ChAT_Cre_X_tdTomato_data[,ii]
  x = seq(y)*dt - dt 
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, baseline=baseline)
  })


ChAT_Cre_X_tdTomato_fits
# using n <- 100 and method <- 'BF.LM':

#          A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -88.305 3.935  3.942 3.938  1.544  8.546 11.360      9.635  945.380 -117.851 11.549  97.056 27.905 10.046 135.483 13.294    101.797 15248.338
# 2   -88.606 6.407  7.100 6.742  2.642 14.660 25.720     16.507 1625.942  -70.643  7.041 259.004 26.094  8.065 359.056 31.497    208.470 20236.213
# 3   -79.712 1.442 31.399 4.655  1.530 43.534 14.273     27.066 2902.916 -210.653 21.857 157.532 50.125 18.260 220.782 27.329    172.395 45616.726
# 4   -69.261 2.956 13.383 5.730  2.152 19.235 10.288     16.999 1422.304  -41.612  3.047 186.091 12.738  3.689 257.977 11.467    142.867  8292.307
# 5  -132.419 2.720  2.723 2.721  1.067  5.905 14.385      6.658  979.584 -177.149 10.615  90.084 25.732  9.257 125.724 16.978     94.239 21234.360
# 6   -51.147 6.279  6.280 6.279  2.461 13.626 22.778     15.362  873.039  -45.042  3.075 117.152 11.494  3.537 162.408 18.579     93.934  5820.659
# 7  -101.555 1.565  7.671 3.125  1.168 10.958 10.127      9.478 1170.867  -86.024  8.587  75.015 21.018  7.543 104.633  7.961     77.880  8539.756
# 8   -23.844 7.794  8.711 8.236  3.228 17.915 20.381     20.168  534.621  -10.815  6.774 595.574 30.671  8.458 825.641 18.800    445.916  6781.617
# 9    -3.217 0.778  2.969 1.412  0.535  4.340 24.680      4.013   15.368  -10.282 18.944 230.012 51.541 17.977 319.489 10.651    220.773  2958.884
# 10 -127.557 3.980  3.991 3.986  1.562  8.648 15.773      9.750 1381.940 -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.689# 10 -127.557 3.980  3.991 3.986  1.562  8.648 15.773      9.750 1381.940 -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.689

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-ChAT_Cre_X_tdTomato_fits[,1], -ChAT_Cre_X_tdTomato_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=100
setwd(svg_path)
DBSCAN_analyse(data=data, eps=100, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('ChAT_Cre_X_tdTomato_peaks', 'ChAT_Cre_X_tdTomato_areas', 'ChAT_Cre_X_tdTomato_fits', 
  'ChAT_Cre_X_tdTomato_summary', 'ChAT_Cre_X_tdTomato_data', 'name', 'analysis_path')

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

identifier <- 'Figure 2'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

# Settings
name <- 'ChAT-Cre X tdTomato'

rdata_path <- file.path(analysis_path, paste0(name, '.RData'))
if (!file.exists(rdata_path)) {
  stop('RData file not found at: ', rdata_path)
}

load(rdata_path)

ChAT_Cre_X_tdTomato_fits
#          A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -88.305 3.935  3.942 3.938  1.544  8.546 11.360      9.635  945.380 -117.851 11.549  97.056 27.905 10.046 135.483 13.294    101.797 15248.338
# 2   -88.606 6.407  7.100 6.742  2.642 14.660 25.720     16.507 1625.942  -70.643  7.041 259.004 26.094  8.065 359.056 31.497    208.470 20236.213
# 3   -79.712 1.442 31.399 4.655  1.530 43.534 14.273     27.066 2902.916 -210.653 21.857 157.532 50.125 18.260 220.782 27.329    172.395 45616.726
# 4   -69.261 2.956 13.383 5.730  2.152 19.235 10.288     16.999 1422.304  -41.612  3.047 186.091 12.738  3.689 257.977 11.467    142.867  8292.307
# 5  -132.419 2.720  2.723 2.721  1.067  5.905 14.385      6.658  979.584 -177.149 10.615  90.084 25.732  9.257 125.724 16.978     94.239 21234.360
# 6   -51.147 6.279  6.280 6.279  2.461 13.626 22.778     15.362  873.039  -45.042  3.075 117.152 11.494  3.537 162.408 18.579     93.934  5820.659
# 7  -101.555 1.565  7.671 3.125  1.168 10.958 10.127      9.478 1170.867  -86.024  8.587  75.015 21.018  7.543 104.633  7.961     77.880  8539.756
# 8   -23.844 7.794  8.711 8.236  3.228 17.915 20.381     20.168  534.621  -10.815  6.774 595.574 30.671  8.458 825.641 18.800    445.916  6781.617
# 9    -3.217 0.778  2.969 1.412  0.535  4.340 24.680      4.013   15.368  -10.282 18.944 230.012 51.541 17.977 319.489 10.651    220.773  2958.884
# 10 -127.557 3.980  3.991 3.986  1.562  8.648 15.773      9.750 1381.940 -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.689

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
for (ii in 1:length(ChAT_Cre_X_tdTomato_summary)){
  fit_plot(traces=ChAT_Cre_X_tdTomato_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=TRUE)
}


