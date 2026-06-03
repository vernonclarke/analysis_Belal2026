# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

# Settings
identifier <- 'Figure S5'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)


data_source_NWB <- TRUE

name <- if (data_source_NWB) 'ctrl' else 'ctrl axg'
ctrl_data <- load_data2(wd=xlsx_path, name=name)[[1]]
dt_range <- ifelse(grepl('samet|^sk', names(ctrl_data), ignore.case = TRUE), 0.04, 0.05)

old_width <- getOption("width")
options(width = 500)
tmp <- ctrl_data[1:10,]
colnames(tmp) <- abbreviate(colnames(tmp), minlength = 8)
print(tmp)


#    s010423_001  s010423_006 s010423_007 e030723_001 e030723_009 e030923_006 e030923_007 e041223_  e041823_ e050223_  e122822_ e123022_003 e123022_004 e123022_007 e123022_011 e123022_013 e123022_014 s221230_005 s221230_008 s122822_010 s122822_013 s122822_00
# 1    0.4693255  0.003964844   1.4076380   1.3584663  2.73106514  -2.5525428  0.66753983  -3.0579  1.055843 1.863288 -3.339204  1.31708626  0.82329556  -1.8999766 -0.33418732  -2.3831674   -2.851821   -3.147856   1.4941929  -1.0806274    2.582153  1.5278015
# 2    0.1423514  0.199277344   1.7467222   0.9323294  0.23106493  -2.1548154 -1.18028580  -2.2829 -1.037907 3.063288 -3.367614  1.00458644 -0.62556819  -1.3907173  0.09081271  -0.7354406   -2.308343   -3.812062   1.2871094  -0.9789022    2.625750  1.3370667
# 3   -2.1028704 -0.008242187   1.8823559   0.2789206 -0.02461628   0.6292758 -0.09332919  -1.1579  1.212093 2.388288 -3.026704  0.47149796 -0.45511395  -2.1546062 -0.40918705  -3.4058947   -1.656169   -4.601929   1.0800258  -0.7245890    2.832833  1.4515076
# 4   -2.9312047 -0.496523437   1.8710531   1.4436932 -1.30302516  -0.5070880  0.83058272  -1.7329  2.805843 1.963288 -2.771023 -0.09835401 -0.45511426  -0.4416436 -1.05918738  -1.6729400   -3.259430   -4.996862   0.9928327  -0.3482056    2.985421  2.4624023
# 5   -1.9066860 -1.070253906   0.4242938   0.3925568 -0.05302569  -1.9275427  0.28710562  -1.2579  2.180843 2.413288 -2.430113  0.89429272  0.11306792  -2.7101618 -2.15918748  -3.8320311   -2.145299   -4.287777   1.5813860  -0.4194132    2.909127  2.1572266
# 6   -1.7976946 -0.838320313  -0.6833813   1.7846025 -1.44507094   0.2599574 -0.55528636  -1.4579  1.680843 2.163288 -2.458523 -0.22703111 -0.39829576  -1.2055321 -2.08418747  -3.4058947   -1.547473   -2.896535   1.5050921  -0.3176880    2.909127  1.1844788
# 7   -2.7132220 -0.630800781  -0.4008111   2.6936935  1.25379197  -0.3082244 -0.58246027  -2.2079  1.649593 3.638288 -5.242614  0.19576243  0.08465851  -1.9231247 -1.05918738  -3.3206674   -3.368125   -2.420819   1.4614955  -0.9890747    3.007220  1.4896545
# 8   -2.0374756 -1.009218750   0.2208433   0.9607386 -0.84848030  -1.1320883  0.55884477  -1.9829  2.993342 3.138288 -2.941477  0.78399798 -0.68238669  -1.9231252 -0.88418737  -3.1786220   -3.368125   -2.950389   2.1808384  -1.0704549    3.029018  0.9174500
# 9   -1.8412912 -1.021425781   1.2041875   1.5005117 -0.62120787   0.1463210 -0.47376340  -1.4829  2.493343 3.013288 -2.856250 -0.11673718 -1.19375006  -0.9971992  0.81581277  -0.7922581   -1.764864   -3.973626   1.7230748  -0.3380330    1.873710  0.4978333
# 10  -0.8603690 -1.106875000   1.4076380   3.2334659 -1.98484372  -0.3650426 -1.12593797  -1.7829  3.399593 2.513288 -2.714205 -0.24541366 -0.88124972  -1.8536808 -1.55918743  -4.3149857   -3.938777   -3.740256   0.9928327  -0.5007935    2.538557  1.3561401

options(width = old_width)

# ==============================================
# METADATA FOR ANAYSIS
# ==============================================
n <- 100
stimulation_time <- 100
baseline <- 100
# dt <- 0.05
method <- 'BF.LM'
weight_method <- 'none'
downsample <- 1
interval=c(0.2, 0.8) # % rise and decay
# MLE as initial estimate then MCMC and the Metropolis–Hastings algorithm to obtain posterior
MLEsettings <- list(iter=1000, metropolis.scale=1.5, fit.attempts=10, RWm=FALSE)
func <- product2N
# fitting to 10% of the peak to avoid any slow components (see English et al., 2012)
rel.decay.fit.limit <- 0.1

# ==============================================
# FITTING DATA
# ==============================================

# to analyse any trace:
# ii=1; analyse_PSC(response=ctrl_data[[ii]], dt=dt_range[ii], n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=ctrl_data[[1]], dt=dt_range[1], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=800, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)  

out1a <- analyse_PSC(response=ctrl_data[[1]], dt=dt_range[1], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=800, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, latency.limit=15)  

out2  <- analyse_PSC(response=ctrl_data[[2]], dt=dt_range[2], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=610, stimulation_time=stimulation_time+12, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)  

out2a <- analyse_PSC(response=ctrl_data[[2]], dt=dt_range[2], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=610, stimulation_time=stimulation_time+12, baseline=baseline, downsample=downsample, latency.limit=25)  

out3  <- analyse_PSC(response=ctrl_data[[3]], dt=dt_range[3], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=930, stimulation_time=stimulation_time+12, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out3a <- analyse_PSC(response=ctrl_data[[3]], dt=dt_range[3], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=800, 
  MLEsettings=MLEsettings, fit.limits=930, stimulation_time=stimulation_time+12, baseline=baseline, downsample=downsample, latency.limit=25)

out4  <- analyse_PSC(response=ctrl_data[[4]], dt=dt_range[4], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=360, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out4a  <- analyse_PSC(response=ctrl_data[[4]], dt=dt_range[4], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=360, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out5  <- analyse_PSC(response=ctrl_data[[5]], dt=dt_range[5], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=290, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out5a  <- analyse_PSC(response=ctrl_data[[5]], dt=dt_range[5], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=290, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out6  <- analyse_PSC(response=ctrl_data[[6]], dt=dt_range[6], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out6a  <- analyse_PSC(response=ctrl_data[[6]], dt=dt_range[6], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out7  <- analyse_PSC(response=ctrl_data[[7]], dt=dt_range[7], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out7a  <- analyse_PSC(response=ctrl_data[[7]], dt=dt_range[7], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out8  <- analyse_PSC(response=ctrl_data[[8]], dt=dt_range[8], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out8a  <- analyse_PSC(response=ctrl_data[[8]], dt=dt_range[8], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out9  <- analyse_PSC(response=ctrl_data[[9]], dt=dt_range[9], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=255, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out9a  <- analyse_PSC(response=ctrl_data[[9]], dt=dt_range[9], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=265, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out10  <- analyse_PSC(response=ctrl_data[[10]], dt=dt_range[10], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=315, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out10a  <- analyse_PSC(response=ctrl_data[[10]], dt=dt_range[10], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=315, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out11  <- analyse_PSC(response=ctrl_data[[11]], dt=dt_range[11], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=280, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), latency.limit=15, fast.constraint=TRUE)

out11a  <- analyse_PSC(response=ctrl_data[[11]], dt=dt_range[11], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=280, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out12  <- analyse_PSC(response=ctrl_data[[12]], dt=dt_range[12], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=305, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out12a  <- analyse_PSC(response=ctrl_data[[12]], dt=dt_range[12], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=305, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out13  <- analyse_PSC(response=ctrl_data[[13]], dt=dt_range[13], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=340, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(30, 500), fast.constraint=TRUE)

out13a  <- analyse_PSC(response=ctrl_data[[13]], dt=dt_range[13], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=340, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out14  <- analyse_PSC(response=ctrl_data[[14]], dt=dt_range[14], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=1000, 
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 1000), fast.constraint=TRUE)

out14a  <- analyse_PSC(response=ctrl_data[[14]], dt=dt_range[14], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out15  <- analyse_PSC(response=ctrl_data[[15]], dt=dt_range[15], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=340, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out15a  <- analyse_PSC(response=ctrl_data[[15]], dt=dt_range[15], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=340, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out16  <- analyse_PSC(response=ctrl_data[[16]], dt=dt_range[16], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=295, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out16a  <- analyse_PSC(response=ctrl_data[[16]], dt=dt_range[16], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=295, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out17  <- analyse_PSC(response=ctrl_data[[17]], dt=dt_range[17], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=245, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out17a  <- analyse_PSC(response=ctrl_data[[17]], dt=dt_range[17], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=245, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out18  <- analyse_PSC(response=ctrl_data[[18]], dt=dt_range[18], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=390, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out18a  <- analyse_PSC(response=ctrl_data[[18]], dt=dt_range[18], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=390, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out19  <- analyse_PSC(response=ctrl_data[[19]], dt=dt_range[19], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=560, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out19a  <- analyse_PSC(response=ctrl_data[[19]], dt=dt_range[19], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=560, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out20  <- analyse_PSC(response=ctrl_data[[20]], dt=dt_range[20], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out20a  <- analyse_PSC(response=ctrl_data[[20]], dt=dt_range[20], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out21  <- analyse_PSC(response=ctrl_data[[21]], dt=dt_range[21], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out21a  <- analyse_PSC(response=ctrl_data[[21]], dt=dt_range[21], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out22  <- analyse_PSC(response=ctrl_data[[22]], dt=dt_range[22], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), fast.constraint=TRUE)

out22a  <- analyse_PSC(response=ctrl_data[[22]], dt=dt_range[22], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample)

out1$BIC  < out1a$BIC
out2$BIC  < out2a$BIC
out3$BIC  < out3a$BIC
out4$BIC  < out4a$BIC
out5$BIC  < out5a$BIC
out6$BIC  < out6a$BIC
out7$BIC  < out7a$BIC
out8$BIC  < out8a$BIC
out9$BIC  < out9a$BIC
out10$BIC < out10a$BIC
out11$BIC < out11a$BIC
out12$BIC < out12a$BIC
out13$BIC < out13a$BIC
out14$BIC < out14a$BIC
out15$BIC < out15a$BIC
out16$BIC < out16a$BIC
out17$BIC < out17a$BIC
out18$BIC < out18a$BIC
out19$BIC < out19a$BIC
out20$BIC < out20a$BIC
out21$BIC < out21a$BIC
out22$BIC < out22a$BIC

# organise outputs
ctrl_summary <- list(
  out1,  out2,  out3,  out4,  out5,  out6,  out7,  out8,
  out9,  out10, out11, out12, out13, out14, out15, out16,
  out17, out18, out19, out20, out21, out22
)

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


stimulation_range <- c(
  105, 107.5, 107.5, 100, 100, 100, 100, 100, 100, 100, 100, 100,
  100, 100, 100, 100, 100, 100, 100, 100, 100, 100
)


ctrl_peaks <- sapply(1:length(ctrl_data), function(ii){
  y <- ctrl_data[,ii]
  x = seq(y)*dt_range[ii] - dt_range[ii]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt_range[ii], stimulation_time=stimulation_range[ii], baseline=baseline)
  })

ctrl_areas <- sapply(1:length(ctrl_data), function(ii){
  y <- ctrl_data[,ii]
  x = seq(y)*dt_range[ii] - dt_range[ii]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]

  keep <- x >= stimulation_range[ii]
  x <- x[keep] - stimulation_range[ii]
  y <- y[keep]

  charge_transfer_fun(y=y, x=x, tmax=600, dt=dt_range[ii], baseline=baseline)
  })


ctrl_fits
# using n <- 100 and method <- 'BF.LM':

#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -220.164 16.024 38.401 24.034  9.287 60.705  1.730     62.565 15808.673 -116.954 150.694 216.575 179.670 70.256 400.291  0.000    491.657  58064.83
# 2  -209.979 10.303 44.102 19.548  7.364 63.698  0.000     57.143 14425.347 -164.028  19.247 214.875  51.005 17.926 298.679  0.000    210.130  44688.42
# 3  -153.456 16.487 21.373 18.719  7.329 41.172  0.000     46.041  7874.352  -80.148  17.072 444.969  57.885 18.642 616.882 24.795    373.687  40617.82
# 4  -608.629  1.763  6.772  3.208  1.216  9.891 11.578      9.131  6619.007 -940.449  12.568  91.881  28.964 10.540 128.716  9.762    100.141 118431.13
# 5  -329.530  1.941  8.654  3.741  1.406 12.454  9.936     11.052  4393.885 -311.263  12.565  80.539  27.659 10.161 113.351  3.938     91.250  35341.33
# 6  -212.335  1.003 21.455  3.223  1.062 29.748 10.265     18.547  5294.290 -119.149  13.879 138.767  35.507 12.601 193.133  1.974    139.266  21355.19
# 7  -333.218  0.752 16.981  2.454  0.803 23.543  9.822     14.559  6538.070 -265.299   8.818 104.141  23.785  8.317 144.684  3.021    100.565  34717.14
# 8  -150.467  2.876  2.877  2.876  1.128  6.241 13.433      7.037  1176.462 -217.347   9.073  61.668  20.387  7.459  86.598 11.030     68.650  18654.95
# 9  -270.866  2.367 11.037  4.639  1.739 15.825 13.274     13.874  4551.671 -166.315  22.162  56.258  34.065 13.138  87.810  6.373     89.449  17142.96
# 10  -75.272  5.621  5.739  5.680  2.226 12.325 10.300     13.895  1162.183 -126.549  10.868  71.323  24.123  8.847 100.286  5.082     80.233  12658.27
# 11 -256.262  1.273 10.077  3.014  1.090 14.085 10.123     10.740  3482.612 -173.385  13.947  77.059  29.108 10.801 109.231  5.240     91.484  19493.13
# 12 -361.732  1.359 39.082  4.729  1.506 54.180  9.513     32.394 15955.727 -129.936  14.823 120.452  35.414 12.783 168.268  3.385    127.480  21000.57
# 13 -139.589 20.377 20.389 20.383  7.990 44.230  5.867     49.866  7734.325 -316.957  16.438  69.854  31.101 11.722 100.963 29.696     90.749  34558.50
# 14  -75.658  2.138 38.211  6.529  2.192 52.991 11.322     34.018  3429.660  -40.238   5.414 234.149  20.876  6.323 324.600  4.683    185.307  10300.30
# 15 -139.108  9.176 40.393 17.596  6.619 58.193 13.286     51.815  8686.588 -128.696  16.741  94.550  35.219 13.050 133.863  4.411    111.455  17660.29
# 16 -286.536  1.326  5.111  2.416  0.915  7.462  9.621      6.883  2349.554 -295.767   6.831  74.922  18.002  6.337 104.162  6.983     73.567  28178.04
# 17 -195.485  1.177 13.549  3.149  1.104 18.828 10.228     13.158  3341.661  -74.968  15.546  73.877  30.687 11.493 105.766  6.127     92.256   8390.45
# 18 -182.153  0.931 27.935  3.276  1.038 38.726 10.700     23.029  5721.460  -85.369  25.288 138.186  52.565 19.518 196.002  0.269    164.644  17257.00
# 19 -319.336  3.805  3.823  3.814  1.495  8.276  9.881      9.331  3310.742 -419.952  11.176 117.487  29.057 10.267 163.413 16.423    116.506  63182.96
# 20  -48.782 12.026 34.055 19.352  7.433 52.016 15.465     51.744  2932.431  -48.029  27.361 148.057  56.671 21.056 210.122  6.654    176.971  10427.11
# 21 -106.675  2.691 42.824  7.947  2.700 59.403 11.893     38.928  5499.626  -37.134  11.076 235.267  35.519 11.706 326.195  0.000    203.589  10160.05
# 22 -182.692 10.688 37.403 18.744  7.135 55.258 11.747     52.243 11279.026 -124.471  30.808 121.050  56.552 21.402 176.356  6.960    161.835  24039.49


# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-ctrl_fits[,1], -ctrl_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=200
setwd(svg_path)
DBSCAN_analyse(data=data, eps=200, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
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

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure S5'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- 'ctrl'

load(paste0(analysis_path, '/', name, '.RData'))

ctrl_fits

#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -220.164 16.024 38.401 24.034  9.287 60.705  1.730     62.565 15808.673 -116.954 150.694 216.575 179.670 70.256 400.291  0.000    491.657  58064.83
# 2  -209.979 10.303 44.102 19.548  7.364 63.698  0.000     57.143 14425.347 -164.028  19.247 214.875  51.005 17.926 298.679  0.000    210.130  44688.42
# 3  -153.456 16.487 21.373 18.719  7.329 41.172  0.000     46.041  7874.352  -80.148  17.072 444.969  57.885 18.642 616.882 24.795    373.687  40617.82
# 4  -608.629  1.763  6.772  3.208  1.216  9.891 11.578      9.131  6619.007 -940.449  12.568  91.881  28.964 10.540 128.716  9.762    100.141 118431.13
# 5  -329.530  1.941  8.654  3.741  1.406 12.454  9.936     11.052  4393.885 -311.263  12.565  80.539  27.659 10.161 113.351  3.938     91.250  35341.33
# 6  -212.335  1.003 21.455  3.223  1.062 29.748 10.265     18.547  5294.290 -119.149  13.879 138.767  35.507 12.601 193.133  1.974    139.266  21355.19
# 7  -333.218  0.752 16.981  2.454  0.803 23.543  9.822     14.559  6538.070 -265.299   8.818 104.141  23.785  8.317 144.684  3.021    100.565  34717.14
# 8  -150.467  2.876  2.877  2.876  1.128  6.241 13.433      7.037  1176.462 -217.347   9.073  61.668  20.387  7.459  86.598 11.030     68.650  18654.95
# 9  -270.866  2.367 11.037  4.639  1.739 15.825 13.274     13.874  4551.671 -166.315  22.162  56.258  34.065 13.138  87.810  6.373     89.449  17142.96
# 10  -40.380  4.581  4.585  4.583  1.797  9.945  4.577     11.212   503.063 -164.359   2.585  71.530   8.904  2.848  99.164 10.161     59.587  13314.97
# 11 -256.262  1.273 10.077  3.014  1.090 14.085 10.123     10.740  3482.612 -173.385  13.947  77.059  29.108 10.801 109.231  5.240     91.484  19493.13
# 12 -361.732  1.359 39.082  4.729  1.506 54.180  9.513     32.394 15955.727 -129.936  14.823 120.452  35.414 12.783 168.268  3.385    127.480  21000.57
# 13 -139.589 20.316 20.450 20.383  7.990 44.230  5.867     49.865  7734.252 -316.958  16.438  69.854  31.101 11.722 100.963 29.696     90.749  34558.54
# 14  -75.658  2.138 38.211  6.529  2.192 52.991 11.322     34.018  3429.660  -40.238   5.414 234.149  20.876  6.323 324.600  4.683    185.307  10300.30
# 15 -139.108  9.176 40.393 17.596  6.619 58.193 13.286     51.815  8686.588 -128.696  16.741  94.550  35.219 13.050 133.863  4.411    111.455  17660.29
# 16 -286.536  1.326  5.111  2.416  0.915  7.462  9.621      6.883  2349.554 -295.767   6.831  74.922  18.002  6.337 104.162  6.983     73.567  28178.04
# 17 -195.485  1.177 13.549  3.149  1.104 18.828 10.228     13.158  3341.661  -74.968  15.546  73.877  30.687 11.493 105.766  6.127     92.256   8390.45
# 18 -182.153  0.931 27.935  3.276  1.038 38.726 10.700     23.029  5721.460  -85.369  25.288 138.186  52.565 19.518 196.002  0.269    164.644  17257.00
# 19 -319.336  3.805  3.823  3.814  1.495  8.276  9.881      9.331  3310.742 -419.952  11.176 117.487  29.057 10.267 163.413 16.423    116.506  63182.96
# 20  -48.782 12.026 34.055 19.352  7.433 52.016 15.465     51.744  2932.431  -48.029  27.361 148.057  56.671 21.056 210.122  6.654    176.971  10427.11
# 21 -106.675  2.691 42.824  7.947  2.700 59.403 11.893     38.928  5499.626  -37.134  11.076 235.267  35.519 11.706 326.195  0.000    203.589  10160.05
# 22 -182.692 10.688 37.403 18.744  7.135 55.258 11.747     52.243 11279.026 -124.471  30.808 121.050  56.552 21.402 176.356  6.960    161.835  24039.49

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
for (ii in 1:length(ctrl_summary)){
  fit_plot(traces=ctrl_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


