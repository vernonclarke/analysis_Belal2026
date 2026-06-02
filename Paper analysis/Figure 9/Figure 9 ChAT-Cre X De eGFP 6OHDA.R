# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')

source(file.path(root_dir, 'R functions', 'setup.R'))

# Settings
identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)

name <- 'ChAT-Cre X De eGFP 6OHDA'
ChAT_Cre_X_De_eGFP_6OHDA_data <- load_data2(wd=xlsx_path, name=name)[[1]]

ChAT_Cre_X_De_eGFP_6OHDA_data[1:10,]

#      23n13001   23n13007   23n14002    23n16015   23d04007
# 1  -1.5887450 -0.5537719 -4.1129148 -0.11442057  0.3717041
# 2  -2.5042723 -0.4927368 -2.3428954  1.59456373 -0.2996826
# 3  -2.0159911 -0.3706665 -0.8780517  0.86214189  0.7989502
# 4  -1.9549560 -0.1875610 -1.9156493 -1.82340486  1.2261962
# 5  -1.8939208 -0.5537719 -3.0142821 -1.70133455  2.2027587
# 6  -2.4432372  0.6669311 -2.9532469 -1.25374343  1.8365478
# 7  -2.3822020  0.7279663 -3.8687742 -1.00960282  1.4093017
# 8  -1.8939208  1.2772827 -2.4039305 -0.48063149  0.6768799
# 9  -0.5511474  0.9721069 -2.6480712  0.04833984 -0.1165771
# 10 -1.1004638 -0.2485962 -3.6856688 -1.00960282 -0.2386474

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
MLEsettings <- list(iter=1000, metropolis.scale=1.5, fit.attempts=10, RWm=TRUE)
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
# analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, downsample=downsample, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 


out1  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#          A1 τrise τdecay tpeak r20_80 d80_20  delay half_width   area1
# fast -1.667 0.805  0.842 0.823  0.323  1.788 14.338      2.015   3.732
# slow -6.113 1.226 21.006 3.699  1.248 29.133  4.475     18.841 153.144

out1a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,1], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#       A1 τrise τdecay tpeak r20_80 d80_20 delay half_width  area1
# 1 -6.366 1.327 20.132 3.863  1.319  27.93 4.541     18.465 155.26

out2  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25, fast.constraint=TRUE) 

#          A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -2.542 0.000  12.605  0.000  0.000  17.474  9.031      8.737   32.039
# slow -8.104 2.252 241.799 10.629  2.851 335.205 24.727    179.020 2047.551

out2a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, half_width_fit_limit=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#       A1 τrise  τdecay tpeak r20_80  d80_20  delay half_width    area1
# 1 -8.239 1.974 239.117 9.546  2.519 331.487 24.628    175.974 2050.384

out3  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -12.837  2.219   2.223  2.221  0.871   4.820 18.845      5.434   77.504
# slow -18.058 24.719 108.719 47.388 17.826 156.642 15.309    139.508 3035.804

out3a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,3], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#        A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -18.325 23.561 113.214 46.703 17.479 161.946 11.457     140.84 3134.006

out4  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500))

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -19.727  4.849   4.849  4.849  1.901  10.522  4.213     11.863  260.018
# slow -21.597 20.611 219.211 53.786 18.985 304.860 12.857    216.785 6050.889

out4a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline)

#        A1  τrise  τdecay  tpeak r20_80 d80_20 delay half_width    area1
# 1 -21.581 22.446 228.904 57.792 20.474  318.5     0    228.634 6358.716

out5  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=25) 

#           A1   τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast  -7.150  22.028  22.073  22.05  8.644  47.848 21.671     53.944  428.592
# slow -14.631 125.706 125.874 125.79 49.309 272.955  4.712    313.522 5002.983

out5a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_6OHDA_data[,5], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, latency.limit=25) 

#        A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -17.583 21.267 297.04 60.398 20.788 412.243 13.838    276.894 6400.404

out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC

# nb zero rise time choose out2a

# organise outputs
ChAT_Cre_X_De_eGFP_6OHDA_summary <- list(out1a, out2a, out3, out4, out5)
names(ChAT_Cre_X_De_eGFP_6OHDA_summary) <- 1:length(ChAT_Cre_X_De_eGFP_6OHDA_summary)

ChAT_Cre_X_De_eGFP_6OHDA_fits <- t(sapply(1:length(ChAT_Cre_X_De_eGFP_6OHDA_summary), function(ii){
  X <- ChAT_Cre_X_De_eGFP_6OHDA_summary[[ii]]$output
  if (dim(X)[1] == 1){
    X <- if (X[3] > 30) c(rep(NA, dim(X)[2]), as.vector(t(X))) else c(as.vector(t(X)), rep(NA, dim(X)[2]))
  }
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(ChAT_Cre_X_De_eGFP_6OHDA_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(ChAT_Cre_X_De_eGFP_6OHDA_fits) <- new_colnames
rownames(ChAT_Cre_X_De_eGFP_6OHDA_fits) <- 1:length(ChAT_Cre_X_De_eGFP_6OHDA_summary)


ChAT_Cre_X_De_eGFP_6OHDA_peaks <- sapply(1:length(ChAT_Cre_X_De_eGFP_6OHDA_data), function(ii){
  y <- ChAT_Cre_X_De_eGFP_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y[x<600], dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


tlimits <- rep(800, length(ChAT_Cre_X_De_eGFP_6OHDA_summary))
tlimits[1] <- 250 # limit fit for 1st trace to avoid incorrect sign

ChAT_Cre_X_De_eGFP_6OHDA_areas <- sapply(1:length(ChAT_Cre_X_De_eGFP_6OHDA_data), function(ii){
  y <- ChAT_Cre_X_De_eGFP_6OHDA_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y[x<tlimits[ii]], x=x[x<tlimits[ii]], dt=dt, baseline=baseline)

  })

cols_to_zero <- colnames(ChAT_Cre_X_De_eGFP_6OHDA_fits) %in% c("A1", "area1")
ChAT_Cre_X_De_eGFP_6OHDA_fits[, cols_to_zero][is.na(ChAT_Cre_X_De_eGFP_6OHDA_fits[, cols_to_zero])] <- 0
ChAT_Cre_X_De_eGFP_6OHDA_fits

#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# 1  -6.366  1.327 20.132  3.863  1.319 27.930  4.541     18.465 155.260   0.000      NA      NA      NA     NA      NA     NA         NA    0.000
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000  -8.239   1.974 239.117   9.546  2.519 331.487 24.628    175.974 2050.384
# 3 -12.837  2.219  2.223  2.221  0.871  4.820 18.845      5.434  77.504 -18.058  24.719 108.719  47.388 17.826 156.642 15.309    139.508 3035.804
# 4 -19.727  4.849  4.849  4.849  1.901 10.522  4.213     11.863 260.018 -21.597  20.611 219.211  53.786 18.985 304.860 12.857    216.785 6050.889
# 5  -7.150 22.028 22.073 22.050  8.644 47.848 21.671     53.944 428.592 -14.631 125.706 125.874 125.790 49.309 272.955  4.712    313.522 5002.983


# # DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
# data <- cbind(-ChAT_Cre_X_De_eGFP_6OHDA_fits[,1], -ChAT_Cre_X_De_eGFP_6OHDA_fits[,10])
# colnames(data) <- c('Afast', 'Aslow')

# # DBSCAN_analyse(data) # eps = 100
# setwd(svg_path)
# DBSCAN_analyse(data, eps=125, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
# setwd(analysis_path)

# List of objects to keep
keep_objects <- c('ChAT_Cre_X_De_eGFP_6OHDA_peaks', 'ChAT_Cre_X_De_eGFP_6OHDA_areas', 'ChAT_Cre_X_De_eGFP_6OHDA_fits', 
  'ChAT_Cre_X_De_eGFP_6OHDA_summary', 'ChAT_Cre_X_De_eGFP_6OHDA_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, ".RData")))

# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- 'ChAT-Cre X De eGFP 6OHDA'

load(paste0(analysis_path, '/', name, '.RData'))
ChAT_Cre_X_De_eGFP_6OHDA_fits

#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# 1  -6.366  1.327 20.132  3.863  1.319 27.930  4.541     18.465 155.260   0.000      NA      NA      NA     NA      NA     NA         NA    0.000
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000  -8.239   1.974 239.117   9.546  2.519 331.487 24.628    175.974 2050.384
# 3 -12.837  2.219  2.223  2.221  0.871  4.820 18.845      5.434  77.504 -18.058  24.719 108.719  47.388 17.826 156.642 15.309    139.508 3035.804
# 4 -19.727  4.849  4.849  4.849  1.901 10.522  4.213     11.863 260.018 -21.597  20.611 219.211  53.786 18.985 304.860 12.857    216.785 6050.889
# 5  -7.150 22.028 22.073 22.050  8.644 47.848 21.671     53.944 428.592 -14.631 125.706 125.874 125.790 49.309 272.955  4.712    313.522 5002.983

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
#   write.csv(ChAT_Cre_X_De_eGFP_6OHDA_fits, file = 'ChAT_Cre_X_De_eGFP_6OHDA_fits.csv', row.names = FALSE)
#   write.csv(ChAT_Cre_X_De_eGFP_6OHDA_peaks, file = 'ChAT_Cre_X_De_eGFP_6OHDA_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
  for (ii in 1:length(ChAT_Cre_X_De_eGFP_6OHDA_summary)){
  traces <- ChAT_Cre_X_De_eGFP_6OHDA_summary[[ii]]$traces
  func <- if (dim(traces)[2]==4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


