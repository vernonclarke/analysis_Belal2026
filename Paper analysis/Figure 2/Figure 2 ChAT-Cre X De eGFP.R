# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 2'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'ChAT-Cre X De eGFP'

# load data
ChAT_Cre_X_De_eGFP_data <- load_data2(wd=xlsx_path, name=name)

# pad matrix as traces are different lengths
max_len <- max(sapply(ChAT_Cre_X_De_eGFP_data, nrow))
ChAT_Cre_X_De_eGFP_data <- do.call(cbind, setNames(lapply(ChAT_Cre_X_De_eGFP_data, function(x) {
  c(x[,1], rep(NA, max_len - nrow(x)))
}), names(ChAT_Cre_X_De_eGFP_data)))

ChAT_Cre_X_De_eGFP_data[1:10,]

#        22316001   22322006   22324018  22329005   24o21000   24o22005
#  [1,] -2.667175 -3.4880980 -0.4484456 -2.482568 -1.1648559  0.6208191
#  [2,] -7.122741 -3.4880980  1.5046793 -3.337060 -4.2166136  1.5973815
#  [3,] -3.948913 -2.5725707  2.5626219 -3.190576 -6.9631955 -0.9660949
#  [4,]  3.009094 -1.7180785  2.2777912 -3.752099 -4.5217893 -2.5835265
#  [5,] -1.385437 -2.5115355  1.9929605 -2.629053 -1.1648559 -0.2947082
#  [6,] -8.526550 -0.1311645  2.3998615 -1.896631 -0.5545044  1.5058288
#  [7,]  2.581848 -0.7415161  3.7019448 -1.676904 -2.0803832  0.2240906
#  [8,] -3.521667 -0.9856567  4.0274656 -2.751123 -1.1648559 -0.2641907
#  [9,] -5.902038 -1.4129028  2.9695230 -3.068506 -1.1648559  0.8954773
# [10,]  2.459778 -1.5960082  2.4405516 -3.630029 -4.2166136 -0.6914367

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
if (method=='MLE'){
  MLEsettings <- list(iter=1000, metropolis.scale=1.5, fit.attempts=10, RWm=TRUE)
  n <- 10
}
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
# analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=605, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit = c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -65.476  9.646  10.032  9.836  3.856  21.350  8.937     24.066  1751.023
# slow -221.043 23.193 168.395 53.323 19.414 235.953 12.500    183.894 51088.809

out2 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=625, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit = c(30, 500)) 
 
#             A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -73.551  4.966   4.972  4.969  1.948  10.782 17.850     12.156   993.475
# slow -205.394 19.079 172.792 47.259 16.911 240.854 20.759    177.732 46654.256

out3 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -26.094 6.899   6.910  6.905  2.707  14.983  5.336     16.892  489.755
# slow -19.980 2.996 140.667 11.783  3.532 195.006 11.173    110.451 3056.172

out4 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#            A1 τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1
# fast  -81.963 1.274  6.060  2.516  0.942  8.675 12.862      7.565   752.274
# slow -162.985 4.713 71.321 13.711  4.682 98.950  6.609     65.450 14088.337

out5 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=560, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500))

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -80.510  4.679   4.685  4.682  1.835  10.159 17.191     11.454  1024.649
# slow -140.766 18.813 141.792 43.811 15.906 198.461 21.985    153.210 27185.635

out6 <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=770, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -80.510  4.679   4.685  4.682  1.835  10.159 17.191     11.454  1024.649
# slow -140.766 18.813 141.792 43.811 15.906 198.461 21.985    153.210 27185.635


# organise outputs
ChAT_Cre_X_De_eGFP_summary <- list(out1, out2, out3, out4, out5, out6)

names(ChAT_Cre_X_De_eGFP_summary) <- 1:length(ChAT_Cre_X_De_eGFP_summary)

ChAT_Cre_X_De_eGFP_fits <- t(sapply(1:length(ChAT_Cre_X_De_eGFP_summary), function(ii){
  X <- ChAT_Cre_X_De_eGFP_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(ChAT_Cre_X_De_eGFP_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(ChAT_Cre_X_De_eGFP_fits) <- new_colnames
rownames(ChAT_Cre_X_De_eGFP_fits) <- 1:length(ChAT_Cre_X_De_eGFP_summary)


ChAT_Cre_X_De_eGFP_peaks <- sapply(1:dim(ChAT_Cre_X_De_eGFP_data)[2], function(ii){
  y <- ChAT_Cre_X_De_eGFP_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


ChAT_Cre_X_De_eGFP_areas <- sapply(1:dim(ChAT_Cre_X_De_eGFP_data)[2], function(ii){
  y <- ChAT_Cre_X_De_eGFP_data[,ii]
  x = seq(y)*dt - dt  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, baseline=baseline)
  })

ChAT_Cre_X_De_eGFP_fits
# using n <- 100 and method <- 'BF.LM':
#        A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -65.476 9.646 10.032 9.836  3.856 21.350  8.937     24.066 1751.023 -221.043 23.193 168.395 53.323 19.414 235.953 12.500    183.894 51088.809
# 2 -73.551 4.967  4.971 4.969  1.948 10.782 17.850     12.156  993.478 -205.394 19.079 172.792 47.259 16.911 240.854 20.759    177.732 46654.255
# 3 -26.094 6.899  6.910 6.905  2.707 14.983  5.336     16.892  489.755  -19.980  2.996 140.667 11.783  3.532 195.006 11.173    110.451  3056.172
# 4 -81.963 1.274  6.060 2.516  0.942  8.675 12.862      7.565  752.274 -162.985  4.713  71.321 13.711  4.682  98.950  6.609     65.450 14088.337
# 5 -80.942 4.700  4.712 4.706  1.845 10.212 17.195     11.513 1035.417 -140.761 18.779 141.793 43.759 15.885 198.452 22.078    153.133 27174.772
# 6 -88.302 5.124  5.128 5.126  2.009 11.123 18.266     12.540 1230.388 -177.971 27.222 232.987 66.176 23.789 325.109  5.281    243.187 55085.135

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
data <- cbind(-ChAT_Cre_X_De_eGFP_fits[,1], -ChAT_Cre_X_De_eGFP_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps = 100
setwd(svg_path)
DBSCAN_analyse(data, eps=125, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('ChAT_Cre_X_De_eGFP_peaks', 'ChAT_Cre_X_De_eGFP_areas', 'ChAT_Cre_X_De_eGFP_fits', 
  'ChAT_Cre_X_De_eGFP_summary', 'ChAT_Cre_X_De_eGFP_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, ".RData")))

# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 2'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

# Settings
name <- 'ChAT-Cre X De eGFP'

rdata_path <- file.path(analysis_path, paste0(name, '.RData'))
if (!file.exists(rdata_path)) {
  stop('RData file not found at: ', rdata_path)
}

load(rdata_path)

ChAT_Cre_X_De_eGFP_fits

#        A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -65.476 9.646 10.032 9.836  3.856 21.350  8.937     24.066 1751.023 -221.043 23.193 168.395 53.323 19.414 235.953 12.500    183.894 51088.809
# 2 -73.551 4.967  4.971 4.969  1.948 10.782 17.850     12.156  993.478 -205.394 19.079 172.792 47.259 16.911 240.854 20.759    177.732 46654.255
# 3 -26.094 6.899  6.910 6.905  2.707 14.983  5.336     16.892  489.755  -19.980  2.996 140.667 11.783  3.532 195.006 11.173    110.451  3056.172
# 4 -81.963 1.274  6.060 2.516  0.942  8.675 12.862      7.565  752.274 -162.985  4.713  71.321 13.711  4.682  98.950  6.609     65.450 14088.337
# 5 -80.942 4.700  4.712 4.706  1.845 10.212 17.195     11.513 1035.417 -140.761 18.779 141.793 43.759 15.885 198.452 22.078    153.133 27174.772
# 6 -88.302 5.124  5.128 5.126  2.009 11.123 18.266     12.540 1230.388 -177.971 27.222 232.987 66.176 23.789 325.109  5.281    243.187 55085.135
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
  for (ii in 1:length(ChAT_Cre_X_De_eGFP_summary)){
  fit_plot(traces=ChAT_Cre_X_De_eGFP_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


