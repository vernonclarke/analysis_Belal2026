# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

# Settings
identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)

name <- 'ChAT-Cre X tdTomato 6OHDA'
ChAT_Cre_X_tdTomato_6OHDA_data <- load_data2(wd=xlsx_path, name=name)[[1]]
ChAT_Cre_X_tdTomato_6OHDA_data[1:10,]

#      23n08005   23n08007   23n08010    23n08020   23n09001   23n09013   23n16010    23n22029     23n22031
# 1  -0.4714762 -0.8011474  0.1453369  0.06722005  1.5424194  0.1277262 -0.2062988 -0.24798583  0.440592427
# 2   0.3830159  0.3381754 -0.6359131  0.02652995  1.2372436  0.2497965 -1.3049316 -0.55316159  0.542317683
# 3   1.1968180 -0.3128662 -0.3917724 -0.09554036  0.2606811  0.2497965  1.1364746  0.17926025  0.359212223
# 4   0.8306071 -1.3301188  1.0486572 -0.42106118 -0.4107055  0.6160075  1.5026855  0.11822509  0.684733040
# 5   0.5864664 -0.4756266  0.4871338 -1.31624343  2.0307006  2.0808511  1.9909667  0.05718994  0.094726558
# 6   0.5864664 -1.2080484 -0.3673584 -1.47900384  3.2514037  0.1277262 -0.2062988 -1.89593497 -0.291829413
# 7   0.7085368 -1.2080484 -1.3683349  0.14860025  0.5048218 -0.4419352 -1.4270019 -1.28558344  0.176106762
# 8   0.9119873 -1.4928792 -1.0265380  0.92171220  1.2982787 -0.3605550 -0.6335449  1.09478755  0.562662734
# 9  -0.7156168 -0.7197672 -0.4161865  0.39274087  1.0541381 -0.8081461  1.2585449  1.52203362 -0.006998698
# 10 -0.6749267 -0.2721761 -0.2452881 -0.33968097 -0.4107055 -0.6046956  2.5402831 -0.67523190 -0.047688800

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
# analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, downsample=downsample, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

# nb avoid spontaneous event
out1  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,1], dt=dt, n=1000, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=750, half_width_fit_limit=1000,  downsample=10, stimulation_time=stimulation_time+20, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500), latency.limit=25) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -12.480  6.041   6.078  6.059  2.375  13.149 10.384     14.824   205.556
# slow -53.077 26.912 276.899 69.487 24.598 385.238 19.517    275.984 18888.979

out1a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,1], dt=dt, n=1000, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=750, half_width_fit_limit=1000,  downsample=10, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#        A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# 1 -49.902 127.711 127.733 127.722 50.066 277.147 14.206    319.105 17325.28

out2  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=775, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500), latency.limit=25)  
 
#             A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width    area1
# fast   -5.910   0.00   1.462   0.000  0.000   2.027 5.500      1.013     8.64
# slow -102.055 136.82 161.160 148.327 58.114 323.594 0.474    382.595 41285.98

out2a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=775, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)  

#         A1   τrise  τdecay   tpeak r20_80  d80_20 delay half_width    area1
# 1 -104.272 145.881 145.949 145.915 57.198 316.623 1.332    374.162 41358.18

out3  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
   MLEsettings=MLEsettings, fit.limits=450, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25)

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -13.210  6.387   6.394  6.391  2.505  13.867  5.154     15.634  229.482
# slow -27.725 18.756 133.055 42.778 15.602 186.571 12.826    146.276 5087.855

out3a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,3], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
   MLEsettings=MLEsettings, fit.limits=450, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)

#        A1 τrise  τdecay  tpeak r20_80 d80_20 delay half_width    area1
# 1 -27.785 23.39 131.535 49.131 18.209 186.27 2.536    155.268 5309.651

out4  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
   MLEsettings=MLEsettings, fit.limits=300, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)

#          A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -5.708 1.087  10.934  2.786  0.988  15.216  3.815     10.956   80.526
# slow -2.411 4.501 423.278 20.672  5.648 586.788 13.059    315.665 1071.680

out4a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
   MLEsettings=MLEsettings, fit.limits=300, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline)

#       A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width   area1
# 1 -3.685  0.31 188.582 1.993   0.42 261.431 3.831    132.807 702.259

out5  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500), latency.limit=25)

#           A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# fast -13.978   3.41   3.414   3.412  1.338   7.404 5.992      8.348   129.651
# slow -29.445 144.23 144.512 144.371 56.593 313.274 0.000    369.327 11555.395

out5a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,5], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25)

#        A1   τrise  τdecay   tpeak r20_80  d80_20 delay half_width    area1
# 1 -29.592 135.487 152.597 143.703 56.316 312.712     0     367.66 11579.75

out6  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -22.285  6.343   6.344  6.344  2.487  13.765  4.363     15.519  384.272
# slow -28.608 19.438 266.530 54.897 18.928 369.934 24.772    249.379 9368.700

out6a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,6], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#       A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -26.95 37.115 260.044 84.286 30.768 364.785     0    286.917 9691.075

out7  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30, 500), latency.limit=25) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -11.396  9.582   9.815  9.698  3.801  21.046 5.664     23.725  300.440
# slow  -5.356 62.251 113.101 82.675 32.193 192.663 3.600    208.138 1258.278

out7a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,7], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#        A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -10.046 2.574 137.394 10.433  3.076 190.469 4.977    106.649 1489.097

out8  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=200, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1
# fast -6.091  3.959  4.612  4.269  1.673  9.307  3.736     10.463  70.883
# slow -3.474 28.442 36.264 32.037 12.545 70.347 10.500     78.745 304.736

out8a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,8], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=200, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#       A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width   area1
# 1 -4.744 0.711 100.329 3.544  0.915 139.085 3.805      73.33 493.064

out9  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#          A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width   area1
# fast -0.972  20.145  20.207  20.176  7.909  43.780  4.530     49.358  53.287
# slow -1.444 114.863 115.927 115.393 45.233 250.399 18.086    284.434 453.001

out9a  <- analyse_PSC(response=ChAT_Cre_X_tdTomato_6OHDA_data[,9], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#      A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width   area1
# 1 -1.56 26.181 391.284 75.881 25.944 542.886   1.1     359.93 741.008

out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC
out6$BIC < out6a$BIC
out7$BIC < out7a$BIC
out8$BIC < out8a$BIC
out9$BIC < out9a$BIC

# organise outputs
ChAT_Cre_X_tdTomato_6OHDA_summary <- list(out1, out2a, out3, out4, out5, out6, out7, out8, out9)
names(ChAT_Cre_X_tdTomato_6OHDA_summary) <- 1:length(ChAT_Cre_X_tdTomato_6OHDA_summary)

ChAT_Cre_X_tdTomato_6OHDA_fits <- t(sapply(1:length(ChAT_Cre_X_tdTomato_6OHDA_summary), function(ii){
  out <- ChAT_Cre_X_tdTomato_6OHDA_summary[[ii]]
  
  if (identical(out, NA)) {
    return(rep(NA, 18))
  }
  
  X <- out$output
  if (dim(X)[1] == 1){
    X <- if (X[3] > 30) c(rep(NA, dim(X)[2]), as.vector(t(X))) else c(as.vector(t(X)), rep(NA, dim(X)[2]))
  }
  as.vector(t(X))
}))


# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(ChAT_Cre_X_tdTomato_6OHDA_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(ChAT_Cre_X_tdTomato_6OHDA_fits) <- new_colnames
rownames(ChAT_Cre_X_tdTomato_6OHDA_fits) <- 1:length(ChAT_Cre_X_tdTomato_6OHDA_summary)


ChAT_Cre_X_tdTomato_6OHDA_peaks <- sapply(1:length(ChAT_Cre_X_tdTomato_6OHDA_data), function(ii){
  y <- ChAT_Cre_X_tdTomato_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y[x<600], dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


tlimits <- rep(1000, length(ChAT_Cre_X_tdTomato_6OHDA_summary))
tlimits[c(2, 4, 5, 6, 7,8)] <- c(1500, 400, 900, 1500, 400, 300) # limit fit for 1st trace to avoid incorrect sign

ChAT_Cre_X_tdTomato_6OHDA_areas <- sapply(1:length(ChAT_Cre_X_tdTomato_6OHDA_data), function(ii){
  y <- ChAT_Cre_X_tdTomato_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y[x<tlimits[ii]], x=x[x<tlimits[ii]], dt=dt, baseline=baseline)

  })

cols_to_zero <- colnames(ChAT_Cre_X_tdTomato_6OHDA_fits) %in% c("A1", "area1")
ChAT_Cre_X_tdTomato_6OHDA_fits[, cols_to_zero][is.na(ChAT_Cre_X_tdTomato_6OHDA_fits[, cols_to_zero])] <- 0

ChAT_Cre_X_tdTomato_6OHDA_fits

#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1 -12.480  6.041  6.078  6.059  2.375 13.149 10.384     14.824 205.556  -53.077  26.912 276.899  69.487 24.598 385.238 19.517    275.984 18888.979
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000 -104.272 145.881 145.949 145.915 57.198 316.623  1.332    374.162 41358.180
# 3 -13.210  6.387  6.394  6.391  2.505 13.867  5.154     15.634 229.482  -27.725  18.756 133.055  42.778 15.602 186.571 12.826    146.276  5087.855
# 4  -5.708  1.087 10.934  2.786  0.988 15.216  3.815     10.956  80.526   -2.411   4.501 423.278  20.672  5.648 586.788 13.059    315.665  1071.680
# 5 -13.978  3.410  3.414  3.412  1.338  7.404  5.992      8.348 129.651  -29.445 144.230 144.512 144.371 56.593 313.274  0.000    369.327 11555.395
# 6 -22.285  6.343  6.344  6.344  2.487 13.765  4.363     15.519 384.272  -28.608  19.438 266.530  54.897 18.928 369.934 24.772    249.379  9368.700
# 7 -11.396  9.582  9.815  9.698  3.801 21.046  5.664     23.725 300.440   -5.356  62.251 113.101  82.675 32.193 192.663  3.600    208.138  1258.278
# 8  -6.091  3.959  4.612  4.269  1.673  9.307  3.736     10.463  70.883   -3.474  28.442  36.264  32.037 12.545  70.347 10.500     78.745   304.736
# 9  -0.972 20.145 20.207 20.176  7.909 43.780  4.530     49.358  53.287   -1.444 114.863 115.927 115.393 45.233 250.399 18.086    284.434   453.001


# # DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
# data <- cbind(-ChAT_Cre_X_tdTomato_6OHDA_fits[,1], -ChAT_Cre_X_tdTomato_6OHDA_fits[,10])
# colnames(data) <- c('Afast', 'Aslow')

# # DBSCAN_analyse(data) # eps = 100
# setwd(svg_path)
# DBSCAN_analyse(data, eps=125, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
# setwd(analysis_path)

# List of objects to keep
keep_objects <- c('ChAT_Cre_X_tdTomato_6OHDA_peaks', 'ChAT_Cre_X_tdTomato_6OHDA_areas', 'ChAT_Cre_X_tdTomato_6OHDA_fits', 
  'ChAT_Cre_X_tdTomato_6OHDA_summary', 'ChAT_Cre_X_tdTomato_6OHDA_data', 'name', 'analysis_path')

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

identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- 'ChAT-Cre X tdTomato 6OHDA'

load(paste0(analysis_path, '/', name, '.RData'))

ChAT_Cre_X_tdTomato_6OHDA_fits

#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1 -12.480  6.041  6.078  6.059  2.375 13.149 10.384     14.824 205.556  -53.077  26.912 276.899  69.487 24.598 385.238 19.517    275.984 18888.979
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000 -104.272 145.881 145.949 145.915 57.198 316.623  1.332    374.162 41358.180
# 3 -13.210  6.387  6.394  6.391  2.505 13.867  5.154     15.634 229.482  -27.725  18.756 133.055  42.778 15.602 186.571 12.826    146.276  5087.855
# 4  -5.708  1.087 10.934  2.786  0.988 15.216  3.815     10.956  80.526   -2.411   4.501 423.278  20.672  5.648 586.788 13.059    315.665  1071.680
# 5 -13.978  3.410  3.414  3.412  1.338  7.404  5.992      8.348 129.651  -29.445 144.230 144.512 144.371 56.593 313.274  0.000    369.327 11555.395
# 6 -22.285  6.343  6.344  6.344  2.487 13.765  4.363     15.519 384.272  -28.608  19.438 266.530  54.897 18.928 369.934 24.772    249.379  9368.700
# 7 -11.396  9.582  9.815  9.698  3.801 21.046  5.664     23.725 300.440   -5.356  62.251 113.101  82.675 32.193 192.663  3.600    208.138  1258.278
# 8  -6.091  3.959  4.612  4.269  1.673  9.307  3.736     10.463  70.883   -3.474  28.442  36.264  32.037 12.545  70.347 10.500     78.745   304.736
# 9  -0.972 20.145 20.207 20.176  7.909 43.780  4.530     49.358  53.287   -1.444 114.863 115.927 115.393 45.233 250.399 18.086    284.434   453.001

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
#   write.csv(ChAT_Cre_X_tdTomato_6OHDA_fits, file = 'ChAT_Cre_X_tdTomato_6OHDA_fits.csv', row.names = FALSE)
#   write.csv(ChAT_Cre_X_tdTomato_6OHDA_peaks, file = 'ChAT_Cre_X_tdTomato_6OHDA_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
for (ii in seq_along(ChAT_Cre_X_tdTomato_6OHDA_summary)) {
  out <- ChAT_Cre_X_tdTomato_6OHDA_summary[[ii]]
  # skip if entry is not a list or doesn't contain 'traces'
  if (!is.list(out) || is.null(out$traces)) next
  traces <- out$traces
  func <- if (ncol(traces) == 4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


