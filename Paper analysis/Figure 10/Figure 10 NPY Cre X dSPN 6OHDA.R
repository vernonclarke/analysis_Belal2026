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

name <- 'NPY Cre X dSPN 6OHDA' # dSPN data
NPY_Cre_X_dSPN_6OHDA_data <- load_data2(wd=xlsx_path, name=name)[[1]]

NPY_Cre_X_dSPN_6OHDA_data[1:10,]

#       23406002    23406004   23406005   23406007   23406008   23406010   23o30006 23o31021
# 1   0.29165648 -1.32666009 -0.1302490 -0.1090698 -0.6156158 -0.1953430 -0.5953613 2.859090
# 2   0.74942013 -1.40804030 -0.7100830 -0.3532104 -1.2564849  0.1098328 -0.5221191 2.859090
# 3   0.41372679 -1.02148433 -0.2828369  1.0505981 -1.0886383  0.7201843 -0.8639160 3.265991
# 4  -0.07455444 -0.22802733 -0.6795654  0.8064575 -0.9970855  0.4760437 -0.7662597 2.859090
# 5   0.16958617 -0.06526692 -0.6795654 -0.6787313 -1.0123443 -0.6225891  0.4056152 3.469442
# 6   0.29165648 -0.87906897 -1.1373290 -0.6380411 -0.5850982 -1.8127746  0.6497558 2.899780
# 7  -0.04403686 -1.08251948 -1.3509521  0.8268025 -0.3409576 -2.0569152  0.9427246 2.818400
# 8  -0.22714232 -0.65527341 -1.0152587  0.9895629 -0.9360504 -1.1719055  1.0159667 2.859090
# 9  -0.28817748 -0.12630208 -0.5880127 -0.2311401 -1.2717437 -0.2563782 -0.4244629 2.981160
# 10 -0.47128294 -0.77734371 -0.1912842  0.2367960 -1.3022613 -1.2634582 -0.5465332 3.143921

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
# analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=540, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -36.672  6.523  62.230 16.436  5.856  86.669 6.720     63.177 2971.913
# slow -15.339 20.347 276.321 57.297 19.774 383.542 3.825    259.046 5215.046

out2  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=710, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20) 

#           A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# fast -37.024 27.708  59.162  39.533 15.331  96.217 6.835    101.324  4272.994
# slow -27.455 59.497 241.521 110.603 41.787 350.658 4.932    320.232 10482.424

out3  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast -48.507  9.638   9.660  9.649  3.782  20.937 5.275     23.605  1272.256
# slow -92.606 16.090 111.171 36.363 13.288 156.023 8.568    123.163 14278.678

out4  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast  -9.857  3.183   3.189  3.186  1.249   6.914 13.253      7.795   85.374
# slow -48.987 14.055 131.626 35.199 12.560 183.370  4.138    134.245 8424.807

out5  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=250, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#          A1  τrise τdecay  tpeak r20_80 d80_20 delay half_width   area1
# fast -7.445  9.234 10.721  9.940  3.895 21.667 6.066     24.361 201.719
# slow -4.769 34.549 39.853 37.075 14.528 80.780 8.501     90.847 481.862

out6  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -16.864  3.178  62.249  9.963  3.312  86.315 5.942     54.570 1231.948
# slow  -5.829 22.011 497.027 71.789 23.499 689.094 0.000    426.097 3347.283

out7  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=420, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE) 

#            A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast  -58.909  7.084  10.927  8.730  3.410  19.669 8.634     21.679  1430.995
# slow -251.516 16.411 119.770 37.798 13.757 167.795 5.015    130.601 41302.221

out8  <- analyse_PSC(response=NPY_Cre_X_dSPN_6OHDA_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(100,500))

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -43.639 19.155  49.441 29.651 11.429  76.900 8.476     78.063  3930.24
# slow -69.250 56.866 183.765 96.592 36.899 274.567 4.844    264.776 21525.76

# organise outputs
NPY_Cre_X_dSPN_6OHDA_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8)

names(NPY_Cre_X_dSPN_6OHDA_summary) <- 1:length(NPY_Cre_X_dSPN_6OHDA_summary)

NPY_Cre_X_dSPN_6OHDA_fits <- t(sapply(1:length(NPY_Cre_X_dSPN_6OHDA_summary), function(ii){
  X <- NPY_Cre_X_dSPN_6OHDA_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(NPY_Cre_X_dSPN_6OHDA_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(NPY_Cre_X_dSPN_6OHDA_fits) <- new_colnames
rownames(NPY_Cre_X_dSPN_6OHDA_fits) <- 1:length(NPY_Cre_X_dSPN_6OHDA_summary)

NPY_Cre_X_dSPN_6OHDA_peaks <- sapply(1:length(NPY_Cre_X_dSPN_6OHDA_data), function(ii){
  y <- NPY_Cre_X_dSPN_6OHDA_data[,ii]
  x = seq(y)*dt - dt  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

NPY_Cre_X_dSPN_6OHDA_areas <- sapply(1:length(NPY_Cre_X_dSPN_6OHDA_data), function(ii){
  y <- NPY_Cre_X_dSPN_6OHDA_data[,ii]
  x = seq(y)*dt - dt  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]

  tmax <- if (ii==5) 400 else 1000
  charge_transfer_fun(y=y[x<tmax], x=x[x<tmax], dt=dt, baseline=baseline)
  })


NPY_Cre_X_dSPN_6OHDA_fits
# using n <- 100 and method <- 'BF.LM':

#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# 1 -36.672  6.523 62.230 16.436  5.856 86.669  6.720     63.177 2971.913  -15.339 20.347 276.321  57.297 19.774 383.542 3.825    259.046  5215.046
# 2 -37.024 27.708 59.162 39.533 15.331 96.217  6.835    101.324 4272.994  -27.455 59.497 241.521 110.603 41.787 350.658 4.932    320.232 10482.424
# 3 -48.507  9.638  9.660  9.649  3.782 20.937  5.275     23.605 1272.256  -92.606 16.090 111.171  36.363 13.288 156.023 8.568    123.163 14278.678
# 4  -9.857  3.183  3.189  3.186  1.249  6.914 13.253      7.795   85.374  -48.987 14.055 131.626  35.199 12.560 183.370 4.138    134.245  8424.807
# 5  -7.445  9.234 10.721  9.940  3.895 21.667  6.066     24.361  201.719   -4.769 34.549  39.853  37.075 14.528  80.780 8.501     90.847   481.862
# 6 -16.864  3.178 62.249  9.963  3.312 86.315  5.942     54.570 1231.948   -5.829 22.011 497.027  71.789 23.499 689.094 0.000    426.097  3347.283
# 7 -58.909  7.084 10.927  8.730  3.410 19.669  8.634     21.679 1430.995 -251.516 16.411 119.770  37.798 13.757 167.795 5.015    130.601 41302.221
# 8 -43.639 19.155 49.441 29.651 11.429 76.900  8.476     78.063 3930.240  -69.250 56.866 183.765  96.592 36.899 274.567 4.844    264.776 21525.759


# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-NPY_Cre_X_dSPN_6OHDA_fits[,1], -NPY_Cre_X_dSPN_6OHDA_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=100
setwd(svg_path)
DBSCAN_analyse(data=data, eps=75, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('NPY_Cre_X_dSPN_6OHDA_peaks', 'NPY_Cre_X_dSPN_6OHDA_areas', 'NPY_Cre_X_dSPN_6OHDA_fits', 
  'NPY_Cre_X_dSPN_6OHDA_summary', 'NPY_Cre_X_dSPN_6OHDA_data', 'name', 'analysis_path')

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

name <- 'NPY Cre X dSPN 6OHDA' # dSPN data

load(paste0(analysis_path, '/', name, '.RData'))

NPY_Cre_X_dSPN_6OHDA_fits
#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# 1 -36.672  6.523 62.230 16.436  5.856 86.669  6.720     63.177 2971.913  -15.339 20.347 276.321  57.297 19.774 383.542 3.825    259.046  5215.046
# 2 -37.024 27.708 59.162 39.533 15.331 96.217  6.835    101.324 4272.994  -27.455 59.497 241.521 110.603 41.787 350.658 4.932    320.232 10482.424
# 3 -48.507  9.638  9.660  9.649  3.782 20.937  5.275     23.605 1272.256  -92.606 16.090 111.171  36.363 13.288 156.023 8.568    123.163 14278.678
# 4  -9.857  3.183  3.189  3.186  1.249  6.914 13.253      7.795   85.374  -48.987 14.055 131.626  35.199 12.560 183.370 4.138    134.245  8424.807
# 5  -7.445  9.234 10.721  9.940  3.895 21.667  6.066     24.361  201.719   -4.769 34.549  39.853  37.075 14.528  80.780 8.501     90.847   481.862
# 6 -16.864  3.178 62.249  9.963  3.312 86.315  5.942     54.570 1231.948   -5.829 22.011 497.027  71.789 23.499 689.094 0.000    426.097  3347.283
# 7 -58.909  7.084 10.927  8.730  3.410 19.669  8.634     21.679 1430.995 -251.516 16.411 119.770  37.798 13.757 167.795 5.015    130.601 41302.221
# 8 -43.639 19.155 49.441 29.651 11.429 76.900  8.476     78.063 3930.240  -69.250 56.866 183.765  96.592 36.899 274.567 4.844    264.776 21525.759

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
#   write.csv(NPY_Cre_X_dSPN_6OHDA_fits, file = 'NPY_Cre_X_dSPN_6OHDA_fits.csv', row.names = FALSE)
#   write.csv(NPY_Cre_X_dSPN_6OHDA_peaks, file = 'NPY_Cre_X_dSPN_6OHDA_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
for (ii in 1:length(NPY_Cre_X_dSPN_6OHDA_summary)){
  fit_plot(traces=NPY_Cre_X_dSPN_6OHDA_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


