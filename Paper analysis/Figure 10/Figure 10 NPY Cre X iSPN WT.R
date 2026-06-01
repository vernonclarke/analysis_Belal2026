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
name <- 'NPY Cre X iSPN WT' # iSPN data
NPY_Cre_X_iSPN_WT_data <- load_data2(wd=xlsx_path, name=name)[[1]]

NPY_Cre_X_iSPN_WT_data[1:10,]
#       23331001    23331004   23331006   23d20000   23d20003   23d20010
# 1   0.61978146 -0.05711364  1.1205749 -0.2655640  0.6435791 -0.7521565
# 2   0.54348752 -0.24021910  0.8459167  0.1616821  1.4248290 -0.3045654
# 3   0.23831176  0.11073303 -0.4358215 -0.4486694  0.7168213  1.6078694
# 4  -0.08212280 -0.37754820  0.1440124 -0.2350464 -0.2841553  0.3057861
# 5  -0.34152220 -0.33177183  0.4491882  0.4363403 -0.1132568  0.4278564
# 6  -0.28048705 -0.10289001 -0.5884094 -0.1434936  0.1552978  1.2823486
# 7  -0.05160522 -0.27073668 -1.2597961 -0.2655640  0.4238525  1.0382080
# 8   0.17727660 -0.31651305 -0.7715149  0.1311645  0.4726806 -0.5080159
# 9   0.10098266 -0.53013608 -0.3442688  0.9246215 -0.3329834 -0.5487060
# 10 -0.17367553 -0.22496032  0.1440124  0.1006470 -0.2109131 -0.4266357

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
# analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[[1]][,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=390, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -24.765 3.283   6.864  4.641  1.801  11.226 5.335     11.864  334.248
# slow -18.260 3.599 129.618 13.266  4.111 179.689 9.946    104.572 2621.938

out2 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=320, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -38.914 2.987   7.263 4.508  1.741  11.444 4.673     11.760  525.756
# slow -29.206 1.450 133.800 6.632  1.817 185.487 2.680     99.892 4106.368

out3 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=220, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1
# fast -44.072  3.819  3.824  3.822  1.498  8.293  4.043      9.349 457.844
# slow  -2.458 26.678 26.959 26.818 10.512 58.193 17.682     65.607 179.176

out4 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=440, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)  

#           A1  τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -32.985  2.757  21.545  6.50  2.354  30.127 4.951     23.047  960.950
# slow  -8.676 20.492 293.323 58.63 20.131 407.041 0.900    272.099 3107.883

out5 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=700, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(100, 500), latency.limit=20)  

#            A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width    area1
# fast -191.684  2.162  99.992   8.474  2.545 138.619 6.051     78.625 20861.98
# slow  -62.773 64.643 343.276 132.971 49.464 487.726 0.000    424.679 31742.64

out6 <- analyse_PSC(response=NPY_Cre_X_iSPN_WT_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=375, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)

#           A1 τrise  τdecay tpeak r20_80 d80_20 delay half_width    area1
# fast -41.619 3.400  15.686 6.637  2.490  22.51 5.677     19.791  996.712
# slow -22.638 2.196 181.160 9.809  2.728 251.14 4.084    136.170 4329.315

# organise outputs
NPY_Cre_X_iSPN_WT_summary <- list(out1, out2, out3, out4, out5, out6)

names(NPY_Cre_X_iSPN_WT_summary) <- 1:length(NPY_Cre_X_iSPN_WT_summary)

NPY_Cre_X_iSPN_WT_fits <- t(sapply(1:length(NPY_Cre_X_iSPN_WT_summary), function(ii){
  X <- NPY_Cre_X_iSPN_WT_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(NPY_Cre_X_iSPN_WT_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(NPY_Cre_X_iSPN_WT_fits) <- new_colnames
rownames(NPY_Cre_X_iSPN_WT_fits) <- 1:length(NPY_Cre_X_iSPN_WT_summary)

NPY_Cre_X_iSPN_WT_peaks <- sapply(1:length(NPY_Cre_X_iSPN_WT_data), function(ii){
  y <- NPY_Cre_X_iSPN_WT_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

NPY_Cre_X_iSPN_WT_areas <- sapply(1:length(NPY_Cre_X_iSPN_WT_data), function(ii){
  y <- NPY_Cre_X_iSPN_WT_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y[x<1000], x=x[x<1000], dt=dt, baseline=baseline)
  })


NPY_Cre_X_iSPN_WT_fits
#         A1 τrise τdecay tpeak r20_80  d80_20 delay half_width     area1      A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -24.765 3.283  6.864 4.641  1.801  11.226 5.335     11.864   334.248 -18.260  3.599 129.618  13.266  4.111 179.689  9.946    104.572  2621.938
# 2  -38.914 2.987  7.263 4.508  1.741  11.444 4.673     11.760   525.756 -29.206  1.450 133.800   6.632  1.817 185.487  2.680     99.892  4106.368
# 3  -44.072 3.819  3.824 3.822  1.498   8.293 4.043      9.349   457.844  -2.458 26.678  26.959  26.818 10.512  58.193 17.682     65.607   179.176
# 4  -32.985 2.757 21.545 6.500  2.354  30.127 4.951     23.047   960.950  -8.676 20.492 293.323  58.630 20.131 407.041  0.900    272.099  3107.883
# 5 -191.684 2.162 99.992 8.474  2.545 138.619 6.051     78.625 20861.978 -62.773 64.643 343.276 132.971 49.464 487.726  0.000    424.679 31742.637
# 6  -41.619 3.400 15.686 6.637  2.490  22.510 5.677     19.791   996.712 -22.638  2.196 181.160   9.809  2.728 251.140  4.084    136.170  4329.315

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-NPY_Cre_X_iSPN_WT_fits[,1], -NPY_Cre_X_iSPN_WT_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=100
setwd(svg_path)
DBSCAN_analyse(data=data, eps=50, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('NPY_Cre_X_iSPN_WT_peaks', 'NPY_Cre_X_iSPN_WT_areas', 'NPY_Cre_X_iSPN_WT_fits', 
  'NPY_Cre_X_iSPN_WT_summary', 'NPY_Cre_X_iSPN_WT_data', 'name', 'analysis_path')

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

name <- 'NPY Cre X iSPN WT' # iSPN data

load(paste0(analysis_path, '/', name, '.RData'))

NPY_Cre_X_iSPN_WT_fits
#         A1 τrise τdecay tpeak r20_80  d80_20 delay half_width     area1      A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -24.765 3.283  6.864 4.641  1.801  11.226 5.335     11.864   334.248 -18.260  3.599 129.618  13.266  4.111 179.689  9.946    104.572  2621.938
# 2  -38.914 2.987  7.263 4.508  1.741  11.444 4.673     11.760   525.756 -29.206  1.450 133.800   6.632  1.817 185.487  2.680     99.892  4106.368
# 3  -44.072 3.819  3.824 3.822  1.498   8.293 4.043      9.349   457.844  -2.458 26.678  26.959  26.818 10.512  58.193 17.682     65.607   179.176
# 4  -32.985 2.757 21.545 6.500  2.354  30.127 4.951     23.047   960.950  -8.676 20.492 293.323  58.630 20.131 407.041  0.900    272.099  3107.883
# 5 -191.684 2.162 99.992 8.474  2.545 138.619 6.051     78.625 20861.978 -62.773 64.643 343.276 132.971 49.464 487.726  0.000    424.679 31742.637
# 6  -41.619 3.400 15.686 6.637  2.490  22.510 5.677     19.791   996.712 -22.638  2.196 181.160   9.809  2.728 251.140  4.084    136.170  4329.315

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
for (ii in 1:length(NPY_Cre_X_iSPN_WT_summary)){
  fit_plot(traces=NPY_Cre_X_iSPN_WT_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


