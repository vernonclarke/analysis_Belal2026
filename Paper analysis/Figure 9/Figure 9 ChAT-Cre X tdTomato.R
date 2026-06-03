# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

# Settings
identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)

name <- 'ChAT-Cre X tdTomato' # dSPN data
ChAT_Cre_X_tdTomato_data <- load_data2(wd=xlsx_path, name=name)[[1]]

ChAT_Cre_X_tdTomato_data[1:10,]
#     24610004    24611000    24611006    24613001    24613002    24613003    24710004    24710006    24711000    24711005
# 1  1.3216857 -0.48500974  0.26745604  0.08459472  1.41210931 -0.08306884  0.86849971 -0.65002438 -0.25374348  0.42807615
# 2  1.8404845 -0.33852537 -0.14758300  0.35314940  0.28906249 -0.38824461  0.13607787 -0.24312336 -0.41650389 -0.13344726
# 3  0.2840881 -0.26528319  0.78015133  0.10900878 -0.29687499 -1.30377191 -0.56582639 -0.08036295  0.64143877  0.89194332
# 4  1.3522033  0.36948240  0.48718259  0.52404783 -0.49218748 -0.75445553 -1.38980096 -0.56864418  1.86214184  1.59995110
# 5  0.3756409  0.17416991 -0.88000484 -0.55017087 -0.56542966  0.03900146 -1.51187127 -1.05692541  1.69938143  0.89194332
# 6  0.7418518  0.66245114  0.73132321 -0.98962398  0.50878904  0.25262450 -0.50479124 -0.32450357 -0.05029297  0.03745117
# 7  0.1315002  0.66245114 -0.04992676 -0.74548336  0.92382808  0.52728269 -0.04702759  0.44860838  0.35660806  1.23374018
# 8  0.5587463  0.63803708 -0.02551269  0.52404783 -0.15039062  0.46624754 -0.65737912  0.48929848  0.60074867  1.99057608
# 9  1.4742736 -0.04555664  0.02331543  0.71936032 -0.02832031  0.92401119 -1.05410762  1.22172032  0.31591795  0.30600584
# 10 0.9554748  0.12534179 -0.26965331  0.40197752  1.11914057  0.64935300 -2.03067007  1.54724114  0.84488928 -0.62172849

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
# analyse_PSC(response=ChAT_Cre_X_tdTomato_data[[1]][,'AP5+NBQX+CGP55845A'], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=770, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30, 500))  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -56.038  3.829   3.832  3.830  1.501   8.311 18.404      9.370   583.438
# slow -210.479 24.793 226.340 61.573 22.019 315.451 11.586    232.336 62533.635

out1a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,1], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=770, stimulation_time=stimulation_time, baseline=baseline) 

#         A1  τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# 1 -212.465 23.904 229.231 60.33 21.486 319.227 9.057    232.415 63366.62

out2  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=510, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30, 500))  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -255.239  3.819  16.579  7.285  2.743  23.916 15.642     21.374  6566.786
# slow -433.480 26.791 136.782 54.317 20.255 194.827 17.900    166.685 88198.384

out2a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=510, stimulation_time=stimulation_time, baseline=baseline)  

#         A1  τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# 1 -488.032 16.044 152.631 40.388 14.393 212.58 12.085     155.06 97053.41

out3  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,3], dt=dt, n=300, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=20, fit.limits=615, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -35.507  2.041   2.046  2.043  0.801   4.434 22.969      4.998   197.202
# slow -145.341 24.325 152.419 53.117 19.545 214.718 11.392    173.890 31388.739

out3a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,3], dt=dt, n=300, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=20, fit.limits=615, stimulation_time=stimulation_time, baseline=baseline) 

#        A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -146.91 22.985 154.576 51.458  18.84 217.148 10.813    172.622 31678.96

out4  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,4], dt=dt, n=n, func=func, method=method, weight_method='~y', interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=755, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30, 500)) 

#            A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# fast -291.374 19.742  23.672  21.588  8.457  47.155     0     52.952  17169.17
# slow -532.137 49.593 256.724 101.061 37.653 365.337     0    311.411 202512.35

out4a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method='~y', interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=755, stimulation_time=stimulation_time, baseline=baseline)

#         A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -640.242 14.669 324.671 47.581  15.61 450.139     0    279.182 240676.6

out5  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,5], dt=dt, n=n,func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=495, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -336.909  1.399  12.237  3.426  1.229  17.068 13.833     12.700  5454.644
# slow -508.551 20.904 141.027 46.850 17.149 198.091 10.827    157.342 99980.051

out5a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,5], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=495, stimulation_time=stimulation_time, baseline=baseline) 

#         A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -596.126 5.123 164.191 18.336  5.763 227.619 11.779     134.27 109442.7

out6  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=475, stimulation_time=stimulation_time, baseline=baseline) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -180.028  2.521   6.184  3.819  1.475   9.724 12.428      9.976  2064.51
# slow -238.322 14.853 129.534 36.334 13.042 180.683  9.087    134.542 40866.58

out6a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,6], dt=dt, n=n,func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=475, stimulation_time=stimulation_time, baseline=baseline) 

#         A1 τrise  τdecay tpeak r20_80  d80_20  delay half_width    area1
# 1 -291.883 1.502 143.399  6.92  1.887 198.793 11.944     106.85 43925.24

out7  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,7], dt=dt, n=n,func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=435, stimulation_time=stimulation_time, baseline=baseline, latency.limit=20) 

#           A1  τrise  τdecay tpeak r20_80  d80_20  delay half_width     area1
# fast -80.479  2.686   8.855   4.6  1.756  13.192 17.954     12.659  1198.019
# slow -83.755 17.660 138.627  41.7 15.093 193.821 14.254    148.104 15685.412

out7a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,7], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=435, stimulation_time=stimulation_time, baseline=baseline, latency.limit=20) 

#         A1 τrise  τdecay tpeak r20_80  d80_20  delay half_width    area1
# 1 -110.851 1.399 147.629 6.581   1.77 204.657 17.546    109.401 17110.78

out8  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=1015, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30,500), latency.limit=25, fast.constraint=TRUE)

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -11.571  8.829   9.152  8.989  3.523  19.510  7.458     21.992   282.773
# slow -56.619 12.564 339.811 43.021 13.796 471.092 17.429    283.945 21836.571

out8a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,8], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=1015, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)

#        A1 τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# 1 -57.134 13.04 338.803 44.176 14.233  469.7 13.578    284.652 22053.06

out9  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=685, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=c(30,500), latency.limit=25)

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -8.834 4.909   5.082  4.995  1.958  10.840  7.265     12.220   119.949
# slow -44.599 4.020 208.204 16.180  4.789 288.632 17.469    162.036 10035.976

out9a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,9], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=685, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)

#        A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -44.799 4.358 207.116 17.188  5.144 287.124 16.181     162.44 10081.35

out10  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,10], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=250, stimulation_time=stimulation_time, baseline=baseline)

#           A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -27.649  4.184 10.879  6.497  2.503  16.896 6.483     17.124  546.552
# slow  -8.092 15.383 95.482 33.478 12.326 134.563 0.499    109.250 1097.130

out10a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_data[,10], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, downsample=downsample, fit.limits=250, stimulation_time=stimulation_time, baseline=baseline)

#        A1 τrise τdecay tpeak r20_80 d80_20 delay half_width   area1
# 1 -30.364 1.688 36.705  5.45  1.791 50.891 6.476     31.649 1292.92

out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC
out6$BIC < out6a$BIC
out7$BIC < out7a$BIC
out8$BIC < out8a$BIC
out9$BIC < out9a$BIC
out10$BIC < out10a$BIC

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

ChAT_Cre_X_tdTomato_peaks <- sapply(1:length(ChAT_Cre_X_tdTomato_data), function(ii){
  y <- ChAT_Cre_X_tdTomato_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

ChAT_Cre_X_tdTomato_areas <- sapply(1:length(ChAT_Cre_X_tdTomato_data), function(ii){
  y <- ChAT_Cre_X_tdTomato_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, tmax=2000, baseline=baseline)
  })


ChAT_Cre_X_tdTomato_fits
# using n <- 100 and method <- 'BF.LM':

#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -56.038  3.829  3.832  3.830  1.501  8.311 18.404      9.370   583.438 -210.479 24.793 226.340  61.573 22.019 315.451 11.586    232.336  62533.64
# 2  -255.239  3.819 16.579  7.285  2.743 23.916 15.642     21.374  6566.786 -433.480 26.791 136.782  54.317 20.255 194.827 17.900    166.685  88198.38
# 3   -35.507  2.041  2.046  2.043  0.801  4.434 22.969      4.998   197.202 -145.341 24.325 152.419  53.117 19.545 214.718 11.392    173.890  31388.74
# 4  -291.374 19.742 23.672 21.588  8.457 47.155  0.000     52.952 17169.174 -532.137 49.593 256.724 101.061 37.653 365.337  0.000    311.411 202512.35
# 5  -336.909  1.399 12.237  3.426  1.229 17.068 13.833     12.700  5454.644 -508.551 20.904 141.027  46.850 17.149 198.091 10.827    157.342  99980.05
# 6  -180.028  2.521  6.184  3.819  1.475  9.724 12.428      9.976  2064.510 -238.322 14.853 129.534  36.334 13.042 180.683  9.087    134.542  40866.58
# 7   -80.479  2.686  8.855  4.600  1.756 13.192 17.954     12.659  1198.019  -83.755 17.660 138.627  41.700 15.093 193.821 14.254    148.104  15685.41
# 8   -11.571  8.829  9.152  8.989  3.523 19.510  7.458     21.992   282.773  -56.619 12.564 339.811  43.021 13.796 471.092 17.429    283.945  21836.57
# 9    -8.834  4.909  5.082  4.995  1.958 10.840  7.265     12.220   119.949  -44.599  4.020 208.204  16.180  4.789 288.632 17.469    162.036  10035.98
# 10  -27.649  4.184 10.879  6.497  2.503 16.896  6.483     17.124   546.552   -8.092 15.383  95.482  33.478 12.326 134.563  0.499    109.250   1097.13

(ChAT_Cre_X_tdTomato_fits[,9] + ChAT_Cre_X_tdTomato_fits[,18])/1e3
 #         1          2          3          4          5          6          7          8          9         10 
 # 63.115143  94.765170  31.585931 219.717200 105.434695  42.931089  16.883442  22.119345  10.155925   1.643682 

ChAT_Cre_X_tdTomato_areas
# [1]  -68.398005 -103.043193  -40.404778 -214.003300 -111.293480  -47.052977  -18.031562  -23.121942  -10.447296   -2.513391


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

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- 'ChAT-Cre X tdTomato'

load(paste0(analysis_path, '/', name, '.RData'))

ChAT_Cre_X_tdTomato_fits

#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -56.038  3.829  3.832  3.830  1.501  8.311 18.404      9.370   583.438 -210.479 24.793 226.340  61.573 22.019 315.451 11.586    232.336  62533.64
# 2  -255.239  3.819 16.579  7.285  2.743 23.916 15.642     21.374  6566.786 -433.480 26.791 136.782  54.317 20.255 194.827 17.900    166.685  88198.38
# 3   -35.507  2.041  2.046  2.043  0.801  4.434 22.969      4.998   197.202 -145.341 24.325 152.419  53.117 19.545 214.718 11.392    173.890  31388.74
# 4  -291.374 19.742 23.672 21.588  8.457 47.155  0.000     52.952 17169.174 -532.137 49.593 256.724 101.061 37.653 365.337  0.000    311.411 202512.35
# 5  -336.909  1.399 12.237  3.426  1.229 17.068 13.833     12.700  5454.644 -508.551 20.904 141.027  46.850 17.149 198.091 10.827    157.342  99980.05
# 6  -180.028  2.521  6.184  3.819  1.475  9.724 12.428      9.976  2064.510 -238.322 14.853 129.534  36.334 13.042 180.683  9.087    134.542  40866.58
# 7   -80.479  2.686  8.855  4.600  1.756 13.192 17.954     12.659  1198.019  -83.755 17.660 138.627  41.700 15.093 193.821 14.254    148.104  15685.41
# 8   -11.571  8.829  9.152  8.989  3.523 19.510  7.458     21.992   282.773  -56.619 12.564 339.811  43.021 13.796 471.092 17.429    283.945  21836.57
# 9    -8.834  4.909  5.082  4.995  1.958 10.840  7.265     12.220   119.949  -44.599  4.020 208.204  16.180  4.789 288.632 17.469    162.036  10035.98
# 10  -27.649  4.184 10.879  6.497  2.503 16.896  6.483     17.124   546.552   -8.092 15.383  95.482  33.478 12.326 134.563  0.499    109.250   1097.13

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

# setwd(xlsx_path)
# if (save){
#   write.csv(ChAT_Cre_X_tdTomato_fits, file = 'ChAT_Cre_X_tdTomato_fits.csv', row.names = FALSE)
#   write.csv(ChAT_Cre_X_tdTomato_peaks, file = 'ChAT_Cre_X_tdTomato_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
for (ii in 1:length(ChAT_Cre_X_tdTomato_summary)){
  fit_plot(traces=ChAT_Cre_X_tdTomato_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=TRUE)
}


