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

name <- 'ChAT-Cre X De eGFP'
ChAT_Cre_X_De_eGFP_data <- load_data2(wd=xlsx_path, name=name)[[1]]
ChAT_Cre_X_De_eGFP_data[1:10,]

#      24610000    24610005    24610006    24610007     24611003   24618001
# 1   1.3520914 -0.26652831  1.80830883 -1.33229974  0.922949175 -1.5257975
# 2   2.1658935  0.07526855  0.42484536 -1.60085442  1.460058524 -0.9561360
# 3   0.8638102 -0.19328612 -1.48758945 -1.13698725  0.166113273 -0.5899251
# 4   1.5962320  0.29499510  0.14001464 -0.16042480 -0.004785156 -0.5899251
# 5   1.9624429 -0.09562988  1.23864740  0.03488769  0.629980439 -0.3457845
# 6   0.3755290 -0.68156735 -0.06343587 -0.77077633  0.434667948 -0.3457845
# 7  -0.4789632  0.12409667 -0.67378740 -1.96706534  0.996191359 -0.4678548
# 8   0.5789795  0.34382323 -0.02274577 -2.28444813  1.069433543 -0.1423340
# 9   3.0610757  0.78327633  0.62829587 -0.96608882  0.703222623 -0.4271647
# 10  3.9969481  0.61237790 -0.02274577 -0.57546384  0.581152316 -1.5257975

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
# analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 
# fits merely as guides to measurements

out1  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,1], dt=dt, n=1000, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=730, downsample=10, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -163.623  4.117  25.722  8.981  3.305  36.240 14.935     29.372   5967.36
# slow -502.946 30.551 223.181 70.389 25.616 312.661 17.000    243.298 153869.15

out1a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,1], dt=dt, n=1000, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=730, downsample=10, stimulation_time=stimulation_time, baseline=baseline) 

#         A1  τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# 1 -534.725 24.685 230.348 61.748 22.041 320.92 11.638    235.144 161040.3

out2  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500))  

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area1
# fast  -84.247  6.614   7.668  7.115  2.788  15.507 10.754     17.436   1633.832
# slow -365.648 19.393 255.620 54.116 18.730 354.866 13.300    241.114 115504.782

out2a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline)  

#        A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -369.91 18.927 258.69 53.402 18.419 359.059 9.241    242.202 117632.9

out3  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast -11.687  1.341   5.154  2.440  0.925   7.527 15.678      6.947    96.709
# slow -69.100 34.474 245.407 78.721 28.703 344.074  4.687    269.522 23371.084

out3a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,3], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline) 

#        A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -69.334 33.402 249.572 77.556 28.177 349.402 3.669     270.33 23610.31

out4  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500)) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -72.115  4.205   4.208  4.206  1.649   9.127 20.066      10.29   824.548
# slow -160.548 19.992 306.514 58.385 19.911 425.235 15.842     280.58 59536.032

out4a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1000, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline) 

#         A1  τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# 1 -164.792 15.948 314.707 50.102 16.644 436.37 13.155    275.553 60811.21

out5  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=460, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500))

#            A1  τrise  τdecay tpeak r20_80  d80_20  delay half_width    area1
# fast  -73.373  3.637   3.642  3.64  1.427   7.898 12.663      8.904   725.95
# slow -241.665 22.674 104.786 44.29 16.614 150.345 15.409    132.131 38643.95

out5a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,5], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=460, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline)

#         A1  τrise  τdecay tpeak r20_80  d80_20  delay half_width    area1
# 1 -238.247 27.507 101.319 49.23 18.692 148.713 10.082    138.776 39240.65

out6  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=725, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit = c(30, 500))

#            A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -55.194  3.369   3.372  3.370  1.321   7.313 14.062      8.245   505.645
# slow -143.395 10.360 165.569 30.628 10.402 229.669 16.154    150.386 28566.182

out6a  <- analyse_PSC(response=ChAT_Cre_X_De_eGFP_data[,6], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=725, downsample=downsample, stimulation_time=stimulation_time, baseline=baseline)

#         A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -144.508 10.749 167.052 31.518 10.734 231.745 12.291    152.527 29152.99

out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC
out6$BIC < out6a$BIC

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


ChAT_Cre_X_De_eGFP_peaks <- sapply(1:length(ChAT_Cre_X_De_eGFP_data), function(ii){
  y <- ChAT_Cre_X_De_eGFP_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


ChAT_Cre_X_De_eGFP_areas <- sapply(1:length(ChAT_Cre_X_De_eGFP_data), function(ii){
  y <- ChAT_Cre_X_De_eGFP_data[,ii]
  x = seq(y)*dt - dt
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, tmax=2000,baseline=baseline)
  })



ChAT_Cre_X_De_eGFP_fits

#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -163.623 4.117 25.722 8.981  3.305 36.240 14.935     29.372 5967.360 -502.946 30.551 223.181 70.389 25.616 312.661 17.000    243.298 153869.15
# 2  -84.247 6.614  7.668 7.115  2.788 15.507 10.754     17.436 1633.832 -365.648 19.393 255.620 54.116 18.730 354.866 13.300    241.114 115504.78
# 3  -11.687 1.341  5.154 2.440  0.925  7.527 15.678      6.947   96.709  -69.100 34.474 245.407 78.721 28.703 344.074  4.687    269.522  23371.08
# 4  -72.115 4.205  4.208 4.206  1.649  9.127 20.066     10.290  824.548 -160.548 19.992 306.514 58.385 19.911 425.235 15.842    280.580  59536.03
# 5  -73.373 3.637  3.642 3.640  1.427  7.898 12.663      8.904  725.950 -241.665 22.674 104.786 44.290 16.614 150.345 15.409    132.131  38643.95
# 6  -55.194 3.369  3.372 3.370  1.321  7.313 14.062      8.245  505.645 -143.395 10.360 165.569 30.628 10.402 229.669 16.154    150.386  28566.18


(ChAT_Cre_X_De_eGFP_fits[,9] + ChAT_Cre_X_De_eGFP_fits[,18])/1e3
#         1         2         3         4         5         6 
# 159.94819 117.13861  23.46779  60.36056  39.36990  29.07183

ChAT_Cre_X_De_eGFP_areas
# [1] -168.50505 -148.41056  -30.94637  -67.50860  -46.30922  -38.03937

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
data <- cbind(-ChAT_Cre_X_De_eGFP_fits[,1], -ChAT_Cre_X_De_eGFP_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps = 100
setwd(svg_path)
DBSCAN_analyse(data, eps=200, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
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

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- 'ChAT-Cre X De eGFP'

load(paste0(analysis_path, '/', name, '.RData'))

ChAT_Cre_X_De_eGFP_fits

#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -163.623 4.117 25.722 8.981  3.305 36.240 14.935     29.372 5967.360 -502.946 30.551 223.181 70.389 25.616 312.661 17.000    243.298 153869.15
# 2  -84.247 6.614  7.668 7.115  2.788 15.507 10.754     17.436 1633.832 -365.648 19.393 255.620 54.116 18.730 354.866 13.300    241.114 115504.78
# 3  -11.687 1.341  5.154 2.440  0.925  7.527 15.678      6.947   96.709  -69.100 34.474 245.407 78.721 28.703 344.074  4.687    269.522  23371.08
# 4  -72.115 4.205  4.208 4.206  1.649  9.127 20.066     10.290  824.548 -160.548 19.992 306.514 58.385 19.911 425.235 15.842    280.580  59536.03
# 5  -73.373 3.637  3.642 3.640  1.427  7.898 12.663      8.904  725.950 -241.665 22.674 104.786 44.290 16.614 150.345 15.409    132.131  38643.95
# 6  -55.194 3.369  3.372 3.370  1.321  7.313 14.062      8.245  505.645 -143.395 10.360 165.569 30.628 10.402 229.669 16.154    150.386  28566.18

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
#   write.csv(ChAT_Cre_X_De_eGFP_fits, file = 'ChAT_Cre_X_De_eGFP_fits.csv', row.names = FALSE)
#   write.csv(ChAT_Cre_X_De_eGFP_peaks, file = 'ChAT_Cre_X_De_eGFP_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
  for (ii in 1:length(ChAT_Cre_X_De_eGFP_summary)){
  fit_plot(traces=ChAT_Cre_X_De_eGFP_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


