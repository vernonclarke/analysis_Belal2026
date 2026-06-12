# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 8'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'Control for MCI-Park' # ctrl data
ctrl_data <- load_data2(wd=xlsx_path, name=name)[[1]]
ctrl_data[1:10,]

#      24115002    24115010   24306011   24312006   24312007  24313006    24313010   24313012
# 1   1.1040527 -0.02106933  1.7751769  0.0723470  1.5287109 -4.120044  0.01062012 -2.0571532
# 2  -0.2387207 -0.48493650  2.0193175  1.6592610  0.9671875 -2.801684 -1.30163568 -0.6411377
# 3  -2.6557128  1.19963373  0.5239563  0.1130371  1.2113281 -1.825122 -1.51525872  1.0922607
# 4  -1.8744628  1.34611810 -0.3610534  0.1130371  1.5775390 -2.288989  0.04113769 -0.4214111
# 5  -2.1430175 -0.99763179 -0.6357116  0.2757975  1.2357421 -4.242114  0.86511226 -0.8364502
# 6  -2.4604003  0.46721189 -1.0934753  0.8454589  2.9203124 -6.048755  2.14685048 -1.0805908
# 7  -2.5580565  1.39494622 -0.6357116  2.4730630  1.9437499 -5.365161  2.78771959  0.3354248
# 8  -0.4096191  0.44279783  0.3103332  2.3509927 -1.1080078 -4.193286  1.65856926  1.8979247
# 9   1.4946777  0.56486814  0.1577454  0.5606282  0.4300781 -3.143481  0.37683104  0.7992920
# 10  2.2515136  1.07756343 -0.4220886 -1.1483561  1.6263671 -3.021411 -0.29455565 -0.6167236

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
# analyse_PSC(response=ctrl_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=ctrl_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=625, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.constraint=TRUE)  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -100.262  2.407   7.228  3.968  1.521  10.927 17.439     10.724  1254.885
# slow -148.460 23.715 181.256 55.492 20.126 253.599 14.358    195.086 36547.865

out1a  <- analyse_PSC(response=ctrl_data[,1], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=625, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  

#        A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width   area1
# 1 -158.58 13.938 199.918 39.903 13.698 277.421 11.475    185.378 38706.5

out2  <- analyse_PSC(response=ctrl_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=835, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  fast.decay.limit=c(30,500), fast.constraint=TRUE)  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -89.222  2.283  16.275  5.217  1.902  22.818 17.147     17.868  2000.751
# slow -224.536 38.027 259.234 85.538 31.287 363.991 14.334    288.332 80961.666

out2a  <- analyse_PSC(response=ctrl_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=835, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  

#         A1  τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# 1 -231.402 32.885 270.221 78.86 28.441 377.395 8.804    285.126 83719.97

out3  <- analyse_PSC(response=ctrl_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=470, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  fast.decay.limit=c(30,500), latency.limit=25, fast.constraint=TRUE)  

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -44.309  6.049   6.053  6.051  2.372  13.130 20.968     14.803   728.791
# slow -98.728 20.298 109.626 42.012 15.612 155.599 24.370    131.118 15877.831

out3a  <- analyse_PSC(response=ctrl_data[,3], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=470, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  latency.limit=25)  

#         A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -102.178 19.301 114.131 41.283  15.25 161.207 18.12    132.581 16743.74

out4  <- analyse_PSC(response=ctrl_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=840, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  fast.decay.limit=c(30,500), latency.limit=25, fast.constraint=TRUE)  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -71.638  4.324   9.477  6.241  2.418  15.310 15.765     16.048  1311.611
# slow -144.206 23.122 237.756 59.689 21.131 330.783 12.739    237.006 44070.227

out4a <- analyse_PSC(response=ctrl_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=840, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  latency.limit=25)  

#         A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -156.938 10.352 257.804 34.673 11.225 357.411 12.312    217.882 46283.67

out5  <- analyse_PSC(response=ctrl_data[,5], dt=dt, n=n, func=func, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=940, stimulation_time=stimulation_time, baseline=baseline, downsample=10, 
  fast.decay.limit=c(30,500), latency.limit=25, fast.constraint=TRUE)  

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -30.310  6.983   7.024  7.004  2.745  15.197 22.551     17.134   577.041
# slow -60.592 24.346 241.547 62.129 22.063 336.214 18.367    242.868 18928.709

out5a  <- analyse_PSC(response=ctrl_data[,5], dt=dt, n=n, func=product1N, method='MLE', weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=940, stimulation_time=stimulation_time, baseline=baseline, downsample=10, 
  latency.limit=25)  

#        A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -66.722 9.353 264.118 32.394 10.336 366.153    19    219.442 19922.04

out6  <- analyse_PSC(response=ctrl_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=665, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -112.088  5.167   5.174  5.171  2.027  11.220 16.033     12.649  1575.402
# slow -255.573 25.873 187.838 59.483 21.657 263.198 18.679    205.132 65891.187

out6a  <- analyse_PSC(response=ctrl_data[,6], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=665, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  

#         A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -255.964 28.236 189.081 63.119 23.117 265.661 11.184    211.422 67577.58

out7 <- analyse_PSC(response=ctrl_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=650, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30,500))  

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -54.341  2.729   2.730  2.729  1.070   5.922 12.671      6.677   403.156
# slow -88.457 12.766 164.094 35.350 12.264 227.840 15.420    155.596 18004.560

out7a <- analyse_PSC(response=ctrl_data[,7], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=650, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  

#        A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -88.091 14.762 165.07 39.141 13.754 229.445 9.411    161.367 18432.27

out8 <- analyse_PSC(response=ctrl_data[,8], dt=dt, n=1000, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=695, stimulation_time=stimulation_time, baseline=baseline, downsample=10, 
  latency.limit=25, fast.decay.limit=c(30,500), fast.constraint=TRUE)  

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -42.51  9.455   9.536  9.495  3.722  20.604 13.302     23.229  1097.229
# slow -108.88 18.795 198.672 48.951 17.288 276.316 23.355    196.756 27675.170

out8a <- analyse_PSC(response=ctrl_data[,8], dt=dt, n=1000, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=695, stimulation_time=stimulation_time, baseline=baseline, downsample=10, 
  latency.limit=25)  

#         A1  τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# 1 -110.664 21.805 197.614 54.023  19.33 275.45 12.744    203.228 28744.27

out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC
out6$BIC < out6a$BIC
out7$BIC < out7a$BIC
out8$BIC < out8a$BIC

# organise outputs
ctrl_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8)

names(ctrl_summary) <- 1:length(ctrl_summary)

ctrl_fits <- t(sapply(1:length(ctrl_summary), function(ii){
  X <- ctrl_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(ctrl_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(ctrl_fits) <- new_colnames
rownames(ctrl_fits) <- 1:length(ctrl_summary)

ctrl_peaks <- sapply(1:length(ctrl_data), function(ii){
  y <- ctrl_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

ctrl_areas <- sapply(1:length(ctrl_data), function(ii){
  y <- ctrl_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, tmax=1500, dt=dt, baseline=baseline)
  })


ctrl_fits
# using n <- 100 and method <- 'BF.LM':

#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -100.262 2.407  7.228 3.968  1.521 10.927 17.439     10.724 1254.885 -148.460 23.715 181.256 55.492 20.126 253.599 14.358    195.086 36547.86
# 2  -89.222 2.283 16.275 5.217  1.902 22.818 17.147     17.868 2000.751 -224.536 38.027 259.234 85.538 31.287 363.991 14.334    288.332 80961.67
# 3  -44.309 6.049  6.053 6.051  2.372 13.130 20.968     14.803  728.790  -98.728 20.298 109.626 42.012 15.612 155.599 24.370    131.118 15877.83
# 4  -71.638 4.324  9.477 6.241  2.418 15.310 15.765     16.048 1311.611 -144.206 23.122 237.756 59.689 21.131 330.783 12.739    237.006 44070.23
# 5  -30.310 6.983  7.024 7.004  2.745 15.197 22.551     17.134  577.041  -60.592 24.346 241.547 62.129 22.063 336.214 18.367    242.868 18928.71
# 6 -112.088 5.167  5.174 5.171  2.027 11.220 16.033     12.649 1575.402 -255.573 25.873 187.838 59.483 21.657 263.198 18.679    205.132 65891.19
# 7  -54.341 2.729  2.730 2.729  1.070  5.922 12.671      6.677  403.156  -88.457 12.766 164.094 35.350 12.264 227.840 15.420    155.596 18004.56
# 8  -42.510 9.455  9.536 9.495  3.722 20.604 13.302     23.229 1097.229 -108.880 18.795 198.672 48.951 17.288 276.316 23.355    196.756 27675.17

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-ctrl_fits[,1], -ctrl_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=80
setwd(svg_path)
DBSCAN_analyse(data=data, eps=80, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('ctrl_peaks', 'ctrl_areas', 'ctrl_fits', 
  'ctrl_summary', 'ctrl_data', 'name', 'analysis_path')

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

identifier <- 'Figure 8'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

name <- 'Control for MCI-Park'

load(paste0(analysis_path, '/', name, '.RData'))

ctrl_fits
#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -100.262 2.407  7.228 3.968  1.521 10.927 17.439     10.724 1254.885 -148.460 23.715 181.256 55.492 20.126 253.599 14.358    195.086 36547.86
# 2  -89.222 2.283 16.275 5.217  1.902 22.818 17.147     17.868 2000.751 -224.536 38.027 259.234 85.538 31.287 363.991 14.334    288.332 80961.67
# 3  -44.309 6.049  6.053 6.051  2.372 13.130 20.968     14.803  728.790  -98.728 20.298 109.626 42.012 15.612 155.599 24.370    131.118 15877.83
# 4  -71.638 4.324  9.477 6.241  2.418 15.310 15.765     16.048 1311.611 -144.206 23.122 237.756 59.689 21.131 330.783 12.739    237.006 44070.23
# 5  -30.310 6.983  7.024 7.004  2.745 15.197 22.551     17.134  577.041  -60.592 24.346 241.547 62.129 22.063 336.214 18.367    242.868 18928.71
# 6 -112.088 5.167  5.174 5.171  2.027 11.220 16.033     12.649 1575.402 -255.573 25.873 187.838 59.483 21.657 263.198 18.679    205.132 65891.19
# 7  -54.341 2.729  2.730 2.729  1.070  5.922 12.671      6.677  403.156  -88.457 12.766 164.094 35.350 12.264 227.840 15.420    155.596 18004.56
# 8  -42.510 9.455  9.536 9.495  3.722 20.604 13.302     23.229 1097.229 -108.880 18.795 198.672 48.951 17.288 276.316 23.355    196.756 27675.17

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
for (ii in 1:length(ctrl_summary)){
  fit_plot(traces=ctrl_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


