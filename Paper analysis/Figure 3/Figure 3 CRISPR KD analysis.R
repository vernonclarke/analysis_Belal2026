# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 3'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
setwd(analysis_path)

name <- 'CRISPR delta KD'

CRISPR_delta_KD_data <- load_data2(wd=xlsx_path, name=name)[[1]]
CRISPR_delta_KD_data[1:10,]

#      24523004   24523009   24524007   24524011   24530003   24530007   24531005  24613007   24614003   24614006
# 1   1.4967854 -5.0213621  1.8707275 -4.1111245 -6.4617917  1.2576293 -6.2473548 -1.520996 -3.0432127 -0.6980387
# 2  -1.7584228 -3.5972085  1.5655517 -2.2800699  1.6762288  4.9197385 -0.3472900 -6.098633 -1.6190592  4.1847736
# 3  -5.0136309 -2.5799559 -5.1483152 -4.7214760  0.4555257 -3.9303587 -0.7541910 -2.436523  1.0257975  0.5226644
# 4   1.2933349 -1.1558024 -4.2327879  3.0096434 -4.8341876 -0.8786010 -0.9576416 -4.877929 -2.6363117 -1.3083902
# 5   4.7519936 -5.4282631  3.0914305  1.5854898 -1.1720784 -9.1183467 -1.9748941 -4.572754  4.2810057 -1.5118407
# 6   0.2760823 -0.3420003 -1.4862060 -1.8731689  0.6589762 -0.2682495  2.0941161 -8.234863  1.0257975 -1.5118407
# 7   1.9036864 -1.1558024  8.5845943 -8.3835852 -0.7651774 -3.9303587 -2.1783446 -2.436523  1.4326985  2.1502685
# 8   3.7347410 -1.9696044 -4.5379637 -2.8904214 -0.9686279  3.6990355  4.7389728  2.141113  3.2637531  3.9813231
# 9  -0.5377197  0.6752522  0.3448486  0.7716878  3.3038329  3.3938597  0.4665120 -7.014160 -0.1949056 -1.3083902
# 10 -0.7411702 -2.1730549  5.2276609 -0.4490153 -5.4445391  0.3421020 -0.9576416 -7.014160  5.9086097  0.7261149

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
# analyse_PSC(response=CRISPR_delta_KD_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=CRISPR_delta_KD_data[,1], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=355, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# fast -89.478 4.162  30.244  9.572  3.485  42.377  4.538     33.021 3713.756
# slow -46.720 8.136 106.397 22.649  7.845 147.713 23.169    100.525 6150.049

out2  <- analyse_PSC(response=CRISPR_delta_KD_data[,2], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=250, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1 τrise  τdecay tpeak r20_80  d80_20 delay half_width    area1
# fast -102.475 1.129  18.499 3.362  1.139  25.659 8.671     16.726 2273.531
# slow  -38.202 1.354 160.435 6.522  1.726 222.411 4.466    118.197 6383.152

out3  <- analyse_PSC(response=CRISPR_delta_KD_data[,3], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=210, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1 τrise τdecay tpeak r20_80 d80_20 delay half_width    area1
# fast -111.502 1.021  9.396 2.542  0.909 13.093 8.859      9.625 1373.204
# slow -129.494 1.099 58.502 4.452  1.313 81.102 4.437     45.422 8174.684

out4  <- analyse_PSC(response=CRISPR_delta_KD_data[,4], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=380, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area1
# fast  -82.918  3.561  3.566  3.563  1.397   7.732  5.037      8.717   803.165
# slow -117.639 11.302 88.266 26.641  9.646 123.425 10.692     94.433 14041.880

out5  <- analyse_PSC(response=CRISPR_delta_KD_data[,5], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=235, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#           A1 τrise τdecay tpeak r20_80 d80_20 delay half_width    area1
# fast -67.201 2.567 10.365 4.762   1.80 15.056 9.372     13.718 1102.728
# slow -60.572 1.563 65.992 5.992   1.82 91.485 5.062     52.352 4377.214

out6  <- analyse_PSC(response=CRISPR_delta_KD_data[,6], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=240, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80 d80_20 delay half_width     area1
# fast -330.918  1.483 10.825  3.416  1.243 15.165 9.869     11.803  4911.133
# slow -170.532 13.154 65.507 26.424  9.869 93.463 4.316     80.502 16721.646

out7  <- analyse_PSC(response=CRISPR_delta_KD_data[,7], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=280, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width     area1
# fast -169.987  1.013 13.276  2.822  0.977  18.431 9.240     12.538  2791.281
# slow  -99.571 12.316 80.909 27.346 10.029 113.760 3.892     90.988 11295.804

out8  <- analyse_PSC(response=CRISPR_delta_KD_data[,8], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=180, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit) 

#            A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area1
# fast -107.180  5.780  5.796  5.788  2.269 12.559  4.318     14.159 1686.208
# slow  -23.843 10.163 36.275 17.964  6.831 53.454 20.461     50.285 1419.167

out9  <- analyse_PSC(response=CRISPR_delta_KD_data[,9], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=390, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=c(30,500)) 

#            A1  τrise τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast  -30.233  2.186  2.191  2.188  0.858   4.748 9.400      5.353   179.82
# slow -233.229 14.865 93.223 32.470 11.947 131.321 4.964    106.328 30801.34

out10  <- analyse_PSC(response=CRISPR_delta_KD_data[,10], dt=dt, n=n, func=func, method=method, weight_method=weight_method, interval=interval, downsample=downsample,
  MLEsettings=MLEsettings, fit.limits=255, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, latency.limit=20) 

#           A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width    area1
# fast -75.489 5.465  11.706  7.808  3.028  19.022 4.856     20.020 1721.790
# slow -41.947 4.939 102.363 15.731  5.198 141.928 3.436     88.924 5007.153

# organise outputs
CRISPR_delta_KD_summary <- list(out1, out2, out3, out4, out5, out6, out7, out8, out9, out10)

names(CRISPR_delta_KD_summary) <- 1:length(CRISPR_delta_KD_summary)

CRISPR_delta_KD_fits <- t(sapply(1:length(CRISPR_delta_KD_summary), function(ii){
  X <- CRISPR_delta_KD_summary[[ii]]$output
  X <- if (dim(X)[1] == 1) c(rep(NA, dim(X)[2]), as.vector(t(X))) else as.vector(t(X))
  as.vector(t(X))
  })
)
# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(CRISPR_delta_KD_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(CRISPR_delta_KD_fits) <- new_colnames
rownames(CRISPR_delta_KD_fits) <- 1:length(CRISPR_delta_KD_summary)


CRISPR_delta_KD_peaks <- sapply(1:length(CRISPR_delta_KD_data), function(ii){
  y <- CRISPR_delta_KD_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y[x<900], dt=dt, stimulation_time=stimulation_time, baseline=baseline)
  })


CRISPR_delta_KD_areas <- sapply(1:length(CRISPR_delta_KD_data), function(ii){
  y <- CRISPR_delta_KD_data[,ii]
  x = seq(y)*dt - dt
  
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  charge_transfer_fun(y=y, x=x, dt=dt, baseline=baseline, tmax=900)
  })


CRISPR_delta_KD_fits[, "A1"][is.na(CRISPR_delta_KD_fits[, "A1"])] <- 0
CRISPR_delta_KD_fits[, "area1"][is.na(CRISPR_delta_KD_fits[, "area1"])] <- 0
# using n <- 100 and method <- 'BF.LM':

CRISPR_delta_KD_fits
#          A1 τrise τdecay tpeak r20_80 d80_20 delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -89.478 4.162 30.244 9.572  3.485 42.377 4.538     33.021 3713.756  -46.720  8.136 106.397 22.649  7.845 147.713 23.169    100.525  6150.049
# 2  -102.475 1.129 18.499 3.362  1.139 25.659 8.671     16.726 2273.531  -38.202  1.354 160.435  6.522  1.726 222.411  4.466    118.197  6383.152
# 3  -111.502 1.021  9.396 2.542  0.909 13.093 8.859      9.625 1373.204 -129.494  1.099  58.502  4.452  1.313  81.102  4.437     45.422  8174.684
# 4   -82.918 3.560  3.567 3.563  1.397  7.732 5.037      8.717  803.162 -117.639 11.302  88.266 26.640  9.646 123.425 10.692     94.433 14041.891
# 5   -67.201 2.567 10.365 4.762  1.800 15.056 9.372     13.718 1102.728  -60.572  1.563  65.992  5.992  1.820  91.485  5.062     52.352  4377.214
# 6  -330.918 1.483 10.825 3.416  1.243 15.165 9.869     11.803 4911.133 -170.532 13.154  65.507 26.424  9.869  93.463  4.316     80.502 16721.646
# 7  -169.987 1.013 13.276 2.822  0.977 18.431 9.240     12.538 2791.281  -99.571 12.316  80.909 27.346 10.029 113.760  3.892     90.988 11295.804
# 8  -107.180 5.780  5.796 5.788  2.269 12.559 4.318     14.159 1686.208  -23.843 10.163  36.275 17.964  6.831  53.454 20.461     50.285  1419.167
# 9   -30.233 2.186  2.191 2.188  0.858  4.748 9.400      5.353  179.820 -233.229 14.865  93.223 32.470 11.947 131.321  4.964    106.328 30801.338
# 10  -75.489 5.465 11.706 7.808  3.028 19.022 4.856     20.020 1721.790  -41.947  4.939 102.363 15.731  5.198 141.928  3.436     88.924  5007.153

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
data <- cbind(-CRISPR_delta_KD_fits[,1], -CRISPR_delta_KD_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps = 150
setwd(svg_path)
DBSCAN_analyse(data, eps=150, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# List of objects to keep
keep_objects <- c('CRISPR_delta_KD_peaks', 'CRISPR_delta_KD_areas', 'CRISPR_delta_KD_fits', 
  'CRISPR_delta_KD_summary', 'CRISPR_delta_KD_data', 'name', 'analysis_path')

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

# Paths
repo_root <- normalizePath(
  '~/Documents/Repositories/analysis_Belal2026',
  mustWork = TRUE
)

identifier <- 'Figure 3'
analysis_path <- file.path(repo_root, 'Paper analysis', identifier)

name <- 'CRISPR delta KD'

load(paste0(analysis_path, '/', name, '.RData'))

CRISPR_delta_KD_fits
#          A1 τrise τdecay tpeak r20_80 d80_20 delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -89.478 4.162 30.244 9.572  3.485 42.377 4.538     33.021 3713.756  -46.720  8.136 106.397 22.649  7.845 147.713 23.169    100.525  6150.049
# 2  -102.475 1.129 18.499 3.362  1.139 25.659 8.671     16.726 2273.531  -38.202  1.354 160.435  6.522  1.726 222.411  4.466    118.197  6383.152
# 3  -111.502 1.021  9.396 2.542  0.909 13.093 8.859      9.625 1373.204 -129.494  1.099  58.502  4.452  1.313  81.102  4.437     45.422  8174.684
# 4   -82.918 3.560  3.567 3.563  1.397  7.732 5.037      8.717  803.162 -117.639 11.302  88.266 26.640  9.646 123.425 10.692     94.433 14041.891
# 5   -67.201 2.567 10.365 4.762  1.800 15.056 9.372     13.718 1102.728  -60.572  1.563  65.992  5.992  1.820  91.485  5.062     52.352  4377.214
# 6  -330.918 1.483 10.825 3.416  1.243 15.165 9.869     11.803 4911.133 -170.532 13.154  65.507 26.424  9.869  93.463  4.316     80.502 16721.646
# 7  -169.987 1.013 13.276 2.822  0.977 18.431 9.240     12.538 2791.281  -99.571 12.316  80.909 27.346 10.029 113.760  3.892     90.988 11295.804
# 8  -107.180 5.780  5.796 5.788  2.269 12.559 4.318     14.159 1686.208  -23.843 10.163  36.275 17.964  6.831  53.454 20.461     50.285  1419.167
# 9   -30.233 2.186  2.191 2.188  0.858  4.748 9.400      5.353  179.820 -233.229 14.865  93.223 32.470 11.947 131.321  4.964    106.328 30801.338
# 10  -75.489 5.465 11.706 7.808  3.028 19.022 4.856     20.020 1721.790  -41.947  4.939 102.363 15.731  5.198 141.928  3.436     88.924  5007.153

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
  for (ii in 1:length(CRISPR_delta_KD_summary)){
  traces <- CRISPR_delta_KD_summary[[ii]]$traces
  func <- if (dim(traces)[2]==4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


