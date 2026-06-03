# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 8'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK' # ctrl data
test_data <- load_data2(wd=xlsx_path, name=name)[[1]]
test_data [1:10,]
#       24112003   24112011   24125001 24125007    24125011    24209016
# 1   0.49273679  0.1373779 -0.7027588 1.715963  0.95855708  0.56657712
# 2   0.01462809 -0.7903564  0.7620849 1.858378  0.13458251  1.17692866
# 3  -0.39227293  0.3815185 -0.5196533 2.204244  0.40924070 -0.33674315
# 4  -1.35866286  1.5289794 -1.4962157 2.102519 -1.26922601 -1.16682123
# 5  -1.50107822  1.6266357 -2.8389891 2.407694 -0.81146236 -1.21564936
# 6  -0.85003658 -0.1311768 -3.3883055 1.838033 -0.04852295 -0.11701660
# 7   0.56394447 -0.9368408 -0.2755127 1.207336  0.95855708 -0.33674315
# 8   0.60463457  0.5035888 -0.5806884 1.044576 -0.75042721  0.05388183
# 9  -0.20916747  1.2115966 -1.7403564 1.736308 -0.17059325  0.61540524
# 10 -1.12469477  0.5035888  0.1517334 1.003886  0.13458251  0.73747555

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
# analyse_PSC(response=test_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=test_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500))  # 360 to avoid noise in decay

#          A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -2.656 15.077  16.123 15.588  6.110  33.856 19.802     38.148  112.611
# slow -7.196 18.706 118.182 40.967 15.065 166.429 42.435    134.490 1202.775

out1a  <- analyse_PSC(response=test_data[,1], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)  # 350 to avoid noise in decay

#       A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -7.271 54.588 67.335 60.516 23.703 132.482 23.284    148.567 1202.699

out2  <- analyse_PSC(response=test_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#          A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -4.075 11.410  11.427 11.419  4.476  24.778 17.958     27.935  126.479
# slow -4.534 37.215 351.483 93.460 33.327 489.590  4.865    357.719 2079.268

out2a  <- analyse_PSC(response=test_data[,2], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=700, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25) 

#       A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -5.791 4.669 361.256 20.569  5.772 500.806 14.157    272.668 2214.704

out3  <- analyse_PSC(response=test_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=25)  

#           A1  τrise  τdecay  tpeak r20_80 d80_20  delay half_width    area1
# fast  -3.612  7.977   9.643  8.757  3.430  19.14 13.700     21.485   86.373
# slow -14.271 10.174 248.706 33.908 11.001 344.80 24.525    210.750 4067.683
 
out3a  <- analyse_PSC(response=test_data[,3], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=600, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#        A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -14.357 12.875 243.371 39.957 13.336 337.473 18.261    214.611 4117.506

out4  <- analyse_PSC(response=test_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#           A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast  -9.975  9.982   9.993  9.988  3.915  21.672 24.042     24.433  270.813
# slow -12.413 24.149 161.634 53.974 19.768 227.102 11.682    180.759 2801.849

out4a  <- analyse_PSC(response=test_data[,4], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#        A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -17.402 5.163 153.533 18.125  5.748 212.845 19.343    126.718 3006.494

out5  <- analyse_PSC(response=test_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, 
  half_width_fit_limit=1000, fast.constraint=TRUE, latency.limit=25)  

#          A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -3.465  1.972   8.401  3.735  1.407  12.140  7.893     10.904   45.399
# slow -6.858 17.669 168.749 44.536 15.866 235.014 23.021    171.268 1506.739

out5a  <- analyse_PSC(response=test_data[,5], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#       A1  τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# 1 -6.385 39.963 142.644 70.639 26.861 210.193 4.873    197.733 1494.424

out6  <- analyse_PSC(response=test_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=275, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#          A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width   area1
# fast -1.963  1.935  1.939  1.937  0.759   4.203 10.173      4.739  10.333
# slow -3.243 26.890 80.293 44.227 16.952 121.496 22.870    119.396 451.760

out6a  <- analyse_PSC(response=test_data[,6], dt=dt, n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=275, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)  

#       A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width   area1
# 1 -3.242 27.054 80.032 44.327 16.997 121.295  22.8    119.459 451.526  


out1$BIC < out1a$BIC
out2$BIC < out2a$BIC
out3$BIC < out3a$BIC
out4$BIC < out4a$BIC
out5$BIC < out5a$BIC
out6$BIC < out6a$BIC


# organise outputs
test_summary <- list(out1, out2, out3, out4, out5, out6)

names(test_summary) <- 1:length(test_summary)

test_fits <- t(sapply(1:length(test_summary), function(ii){
  X <- test_summary[[ii]]$output
  X <- if (dim(X)[1] == 1) c(rep(NA, dim(X)[2]), as.vector(t(X))) else as.vector(t(X))
  as.vector(t(X))
  })
)

# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(test_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(test_fits) <- new_colnames
rownames(test_fits) <- 1:length(test_summary)

test_peaks <- sapply(1:length(test_data), function(ii){
  y <- test_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })

test_areas <- sapply(1:length(test_data), function(ii){
  y <- test_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, tmax=1000, dt=dt, baseline=baseline)
  })

test_fits[, "A1"][is.na(test_fits[, "A1"])] <- 0
test_fits[, "area1"][is.na(test_fits[, "area1"])] <- 0
tests_fits

test_fits
# using n <- 100 and method <- 'BF.LM':

#       A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -2.656 15.077 16.123 15.588  6.110 33.856 19.802     38.148 112.611  -7.196 18.706 118.182 40.967 15.065 166.429 42.435    134.490 1202.775
# 2 -4.075 11.410 11.427 11.419  4.476 24.778 17.958     27.935 126.479  -4.534 37.215 351.483 93.460 33.327 489.590  4.865    357.719 2079.268
# 3 -3.612  7.977  9.643  8.757  3.430 19.140 13.700     21.485  86.373 -14.271 10.174 248.706 33.908 11.001 344.800 24.525    210.750 4067.683
# 4 -9.975  9.982  9.993  9.988  3.915 21.672 24.042     24.433 270.813 -12.413 24.149 161.634 53.974 19.768 227.102 11.682    180.759 2801.849
# 5 -3.465  1.972  8.401  3.735  1.407 12.140  7.893     10.904  45.399  -6.858 17.669 168.749 44.536 15.866 235.014 23.021    171.268 1506.739
# 6 -1.963  1.935  1.939  1.937  0.759  4.203 10.173      4.739  10.333  -3.243 26.890  80.293 44.227 16.952 121.496 22.870    119.396  451.760


# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-test_fits[,1], -test_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=6
setwd(svg_path)
DBSCAN_analyse(data=data, eps=6, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('test_peaks', 'test_areas', 'test_fits', 
  'test_summary', 'test_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, '.RData')))


# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

# Paths
repo_root <- normalizePath(
  '~/Documents/Repositories/analysis_Belal2026',
  mustWork = TRUE
)

identifier <- 'Figure 8'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

if (!dir.exists(analysis_path)) {
  stop('Analysis folder not found at: ', analysis_path)
}
name <- 'ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK' 

load(paste0(analysis_path, '/', name, '.RData'))

test_fits
#       A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -2.656 15.077 16.123 15.588  6.110 33.856 19.802     38.148 112.611  -7.196 18.706 118.182 40.967 15.065 166.429 42.435    134.490 1202.775
# 2 -4.075 11.410 11.427 11.419  4.476 24.778 17.958     27.935 126.479  -4.534 37.215 351.483 93.460 33.327 489.590  4.865    357.719 2079.268
# 3 -3.612  7.977  9.643  8.757  3.430 19.140 13.700     21.485  86.373 -14.271 10.174 248.706 33.908 11.001 344.800 24.525    210.750 4067.683
# 4 -9.975  9.982  9.993  9.988  3.915 21.672 24.042     24.433 270.813 -12.413 24.149 161.634 53.974 19.768 227.102 11.682    180.759 2801.849
# 5 -3.465  1.972  8.401  3.735  1.407 12.140  7.893     10.904  45.399  -6.858 17.669 168.749 44.536 15.866 235.014 23.021    171.268 1506.739
# 6 -1.963  1.935  1.939  1.937  0.759  4.203 10.173      4.739  10.333  -3.243 26.890  80.293 44.227 16.952 121.496 22.870    119.396  451.760

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
for (ii in 1:length(test_summary)){
  fit_plot(traces=test_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


