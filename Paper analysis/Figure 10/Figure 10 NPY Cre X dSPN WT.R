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

name <- 'NPY Cre X dSPN WT' # dSPN data
NPY_Cre_X_dSPN_WT_data <- load_data2(wd=xlsx_path, name=name)[[1]]

NPY_Cre_X_dSPN_WT_data[1:10,]
#       23315002    23315010    23315012   23315013   23316006    23316014   23330003   23330008    23d18010      23d18011   23d18013    23d20005   23d20008
# 1  -0.79193111  1.08734126 -1.36856073  2.7578734 -0.8133951  0.24200438  0.7878418 -1.5223632  0.32849120 -0.8242797460 -0.6147868 -0.44318846 -0.4968261
# 2  -1.50400790 -0.05706787 -1.13967890  1.3845825 -1.2202962  0.42510984  0.3605957 -1.3514648  1.24401850 -0.4885864026  0.1583252 -0.17463378 -1.6442870
# 3  -1.09710688 -0.05706787 -0.68191525 -1.0568237 -1.2202962 -0.00213623 -0.1276855 -1.2782226  0.02331543 -0.3054809425  0.6872965 -0.78498531 -2.0837401
# 4  -0.89365637  0.47698972 -0.14785766 -0.9042358 -1.7289224 -0.24627684 -0.4023437 -0.7655273 -0.40393064  0.0607299776  0.9721272  0.04509277  0.3820801
# 5   0.22532144  1.08734126  1.07284541 -1.5145873 -2.0340982 -0.79559322 -0.1887207 -0.3871094  0.87780758 -0.0003051758  0.2803955 -0.07697754  1.4562988
# 6  -0.07985433  0.62957761  0.08102417 -2.1249389 -1.8306477 -0.97869868 -0.1124268 -0.7899414  2.03747549 -0.8242797460 -0.3706461 -0.32111815  1.1633300
# 7   0.02187093  1.62139885 -0.14785766 -2.1249389 -0.9151204 -0.73455807  0.0859375 -1.4369140 -0.46496580 -1.1294555128  0.4024658  1.02165522  0.6018066
# 8  -0.07985433  2.76580797  0.38619993 -0.2938843 -1.8306477  0.54718015 -0.1887207 -1.8763671 -0.40393064 -1.5872191629 -0.2485758  1.04606929  0.3332519
# 9  -0.07985433  2.07916250  0.30990599  0.6216430 -1.6271972  1.03546138 -0.2955322 -1.4857421 -0.03771972 -0.3359985192 -0.8182373 -0.17463378 -0.5456543
# 10 -0.07985433  1.16363520  0.76766964  1.3845825 -1.2202962  0.66925046  0.1164551 -1.2293945 -0.28186034  0.5795287811 -0.1265055  0.04509277 -0.4235840

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
# analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, return.output=FALSE) 

out1  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=330, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -20.740 7.432   7.451  7.441  2.917  16.147 5.666     18.204  419.512
# slow -20.298 3.010 110.691 11.153  3.447 153.451 8.069     89.097 2485.035

out2  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -40.701 3.945  8.127 5.540  2.151  13.344 4.696     14.138  654.040
# slow -11.121 2.316 85.005 8.577  2.652 117.842 0.000     68.435 1045.684

out3  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=200, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#            A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -123.390 1.685  11.723 3.817  1.394  16.448 6.436      12.96 2003.107
# slow  -32.411 1.911 116.823 7.990  2.314 161.951 4.362      89.68 4054.373

out4  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=290, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -179.594  3.031 16.531  6.296  2.338  23.449  6.353     19.708 4344.955
# slow  -37.949 54.203 54.352 54.277 21.276 117.777 19.643    132.783 5599.033

out5  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=325, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#            A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -138.271 2.182   6.505  3.586  1.375   9.845 5.142      9.678 1560.914
# slow  -66.797 2.606 134.418 10.478  3.103 186.343 9.084    104.648 9706.684

out6  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,6], dt=dt, n=10, func=func, method='MLE', weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=250, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)

#           A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width   area1
# fast -15.969  7.562  7.726  7.643  2.996  16.587  4.15     18.699 331.793
# slow  -2.981 20.166 75.856 36.392 13.803 111.060 20.00    103.080 365.297

out7  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=310, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -31.382 3.527  11.098  5.926  2.266  16.649 5.082     16.159  594.049
# slow -20.130 3.323 120.873 12.281  3.801 167.566 6.558     97.411 2693.428

out8  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -16.752 11.340  11.773 11.554  4.529  25.078  4.400     28.268  526.204
# slow -10.555 16.807 117.288 38.115 13.918 164.551 19.846    129.560 1713.322

out9  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
   downsample=downsample, MLEsettings=MLEsettings, fit.limits=550, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, half_width_fit_limit = 1000)

#           A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# fast -20.204  17.458  60.535  30.505 11.616  89.542  7.147     84.853 2024.401
# slow  -6.449 173.471 174.316 173.893 68.165 377.335 37.014    466.359 3048.475

out10 <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,10], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=990, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(100,500))  

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast -19.817  5.441   5.448  5.445  2.134  11.815 4.869     13.320   293.298
# slow -86.816 17.612 317.013 53.898 18.082 439.625 6.030    281.885 32622.547

out11 <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,11], dt=dt, n=10, func=func, method='MLE', weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=650, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)

#           A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast  -6.155  6.563   8.606  7.492  2.933  16.499 6.040     18.436  126.511
# slow -22.931 16.730 228.626 47.201 16.280 317.329 3.767    214.061 6444.854

out12  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,12], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=230, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE)  

#           A1 τrise τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -73.848 2.042 10.417 4.138  1.543  14.838 3.622     12.697 1144.483
# slow -44.866 1.292 81.211 5.438  1.569 112.583 5.364     62.211 3895.942

out13  <- analyse_PSC(response=NPY_Cre_X_dSPN_WT_data[,13], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  downsample=downsample, MLEsettings=MLEsettings, fit.limits=230, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, latency.limit=20)

#           A1 τrise τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -76.051 2.489 11.443 4.853  1.821  16.426 3.981     14.456 1329.930
# slow -32.462 1.959 93.919 7.743  2.314 130.200 5.719     73.601 3310.801

# organise outputs
NPY_Cre_X_dSPN_WT_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12, out13)

names(NPY_Cre_X_dSPN_WT_summary) <- 1:length(NPY_Cre_X_dSPN_WT_summary)

NPY_Cre_X_dSPN_WT_fits <- t(sapply(1:length(NPY_Cre_X_dSPN_WT_summary), function(ii){
  X <- NPY_Cre_X_dSPN_WT_summary[[ii]]$output
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(NPY_Cre_X_dSPN_WT_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(NPY_Cre_X_dSPN_WT_fits) <- new_colnames
rownames(NPY_Cre_X_dSPN_WT_fits) <- 1:length(NPY_Cre_X_dSPN_WT_summary)

NPY_Cre_X_dSPN_WT_peaks <- sapply(1:length(NPY_Cre_X_dSPN_WT_data), function(ii){
  y <- NPY_Cre_X_dSPN_WT_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

NPY_Cre_X_dSPN_WT_areas <- sapply(1:length(NPY_Cre_X_dSPN_WT_data), function(ii){
  y <- NPY_Cre_X_dSPN_WT_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y[x<1000], x=x[x<1000], dt=dt, baseline=baseline)
  })


NPY_Cre_X_dSPN_WT_fits
# using n <- 100 and method <- 'BF.LM':

#          A1  τrise τdecay  tpeak r20_80 d80_20 delay half_width    area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -20.740  7.432  7.451  7.441  2.917 16.147 5.666     18.204  419.512 -20.298   3.010 110.691  11.153  3.447 153.451  8.069     89.097  2485.035
# 2   -40.701  3.945  8.127  5.540  2.151 13.344 4.696     14.138  654.040 -11.121   2.316  85.005   8.577  2.652 117.842  0.000     68.435  1045.684
# 3  -123.390  1.685 11.723  3.817  1.394 16.448 6.436     12.960 2003.107 -32.411   1.911 116.823   7.990  2.314 161.951  4.362     89.680  4054.373
# 4  -179.594  3.031 16.531  6.296  2.338 23.449 6.353     19.708 4344.955 -37.949  54.203  54.352  54.277 21.276 117.777 19.643    132.783  5599.033
# 5  -138.271  2.182  6.505  3.586  1.375  9.845 5.142      9.678 1560.914 -66.797   2.606 134.418  10.478  3.103 186.343  9.084    104.648  9706.684
# 6   -15.969  7.543  7.745  7.643  2.996 16.587 4.150     18.699  331.793  -2.981  20.168  75.825  36.388 13.802 111.021 20.000    103.057   365.236
# 7   -31.382  3.527 11.098  5.926  2.266 16.649 5.082     16.159  594.049 -20.130   3.323 120.873  12.281  3.801 167.566  6.558     97.411  2693.428
# 8   -16.752 11.340 11.773 11.554  4.529 25.078 4.400     28.268  526.204 -10.555  16.807 117.288  38.115 13.918 164.551 19.846    129.560  1713.322
# 9   -20.204 17.458 60.535 30.505 11.616 89.542 7.147     84.853 2024.401  -6.449 173.471 174.316 173.893 68.165 377.335 37.014    466.359  3048.475
# 10  -19.817  5.441  5.448  5.445  2.134 11.815 4.869     13.320  293.298 -86.816  17.612 317.013  53.898 18.082 439.625  6.030    281.885 32622.547
# 11   -6.104  7.287  7.671  7.476  2.930 16.230 6.043     18.292  124.072 -22.956  16.621 228.719  46.993 16.197 317.447  3.774    213.847  6448.170
# 12  -73.848  2.042 10.417  4.138  1.543 14.838 3.622     12.697 1144.483 -44.866   1.292  81.211   5.438  1.569 112.583  5.364     62.211  3895.942
# 13  -76.051  2.489 11.443  4.853  1.821 16.426 3.981     14.456 1329.930 -32.462   1.959  93.919   7.743  2.314 130.200  5.719     73.601  3310.801

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-NPY_Cre_X_dSPN_WT_fits[,1], -NPY_Cre_X_dSPN_WT_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=100
setwd(svg_path)
DBSCAN_analyse(data=data, eps=20, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('NPY_Cre_X_dSPN_WT_peaks', 'NPY_Cre_X_dSPN_WT_areas', 'NPY_Cre_X_dSPN_WT_fits', 
  'NPY_Cre_X_dSPN_WT_summary', 'NPY_Cre_X_dSPN_WT_data', 'name', 'analysis_path')

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

name <- 'NPY Cre X dSPN WT' # dSPN data

load(paste0(analysis_path, '/', name, '.RData'))

NPY_Cre_X_dSPN_WT_fits
#          A1  τrise τdecay  tpeak r20_80 d80_20 delay half_width    area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -20.740  7.432  7.451  7.441  2.917 16.147 5.666     18.204  419.512 -20.298   3.010 110.691  11.153  3.447 153.451  8.069     89.097  2485.035
# 2   -40.701  3.945  8.127  5.540  2.151 13.344 4.696     14.138  654.040 -11.121   2.316  85.005   8.577  2.652 117.842  0.000     68.435  1045.684
# 3  -123.390  1.685 11.723  3.817  1.394 16.448 6.436     12.960 2003.107 -32.411   1.911 116.823   7.990  2.314 161.951  4.362     89.680  4054.373
# 4  -179.594  3.031 16.531  6.296  2.338 23.449 6.353     19.708 4344.955 -37.949  54.203  54.352  54.277 21.276 117.777 19.643    132.783  5599.033
# 5  -138.271  2.182  6.505  3.586  1.375  9.845 5.142      9.678 1560.914 -66.797   2.606 134.418  10.478  3.103 186.343  9.084    104.648  9706.684
# 6   -15.969  7.543  7.745  7.643  2.996 16.587 4.150     18.699  331.793  -2.981  20.168  75.825  36.388 13.802 111.021 20.000    103.057   365.236
# 7   -31.382  3.527 11.098  5.926  2.266 16.649 5.082     16.159  594.049 -20.130   3.323 120.873  12.281  3.801 167.566  6.558     97.411  2693.428
# 8   -16.752 11.340 11.773 11.554  4.529 25.078 4.400     28.268  526.204 -10.555  16.807 117.288  38.115 13.918 164.551 19.846    129.560  1713.322
# 9   -20.204 17.458 60.535 30.505 11.616 89.542 7.147     84.853 2024.401  -6.449 173.471 174.316 173.893 68.165 377.335 37.014    466.359  3048.475
# 10  -19.817  5.441  5.448  5.445  2.134 11.815 4.869     13.320  293.298 -86.816  17.612 317.013  53.898 18.082 439.625  6.030    281.885 32622.547
# 11   -6.104  7.287  7.671  7.476  2.930 16.230 6.043     18.292  124.072 -22.956  16.621 228.719  46.993 16.197 317.447  3.774    213.847  6448.170
# 12  -73.848  2.042 10.417  4.138  1.543 14.838 3.622     12.697 1144.483 -44.866   1.292  81.211   5.438  1.569 112.583  5.364     62.211  3895.942
# 13  -76.051  2.489 11.443  4.853  1.821 16.426 3.981     14.456 1329.930 -32.462   1.959  93.919   7.743  2.314 130.200  5.719     73.601  3310.801
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
#   write.csv(NPY_Cre_X_dSPN_WT_fits, file = 'NPY_Cre_X_dSPN_WT_fits.csv', row.names = FALSE)
#   write.csv(NPY_Cre_X_dSPN_WT_peaks, file = 'NPY_Cre_X_dSPN_WT_peaks.csv', row.names = FALSE)
#   }

setwd(svg_path)
for (ii in 1:length(NPY_Cre_X_dSPN_WT_summary)){
  fit_plot(traces=NPY_Cre_X_dSPN_WT_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


