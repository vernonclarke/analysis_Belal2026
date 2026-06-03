# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

# Settings
identifier <- 'Figure S5'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)

data_source_NWB <- TRUE
name <- if (data_source_NWB) '6OHDA' else '6OHDA axg'
OHDA_data <- load_data2(wd=xlsx_path, name=name)[[1]]
dt_range <- ifelse(grepl('samet|^sk', names(OHDA_data), ignore.case = TRUE), 0.04, 0.05)

old_width <- getOption("width")
options(width = 500)
tmp <- OHDA_data[1:10,]
colnames(tmp) <- abbreviate(colnames(tmp), minlength = 8)
print(tmp)

#    s010423_000 s010423_009 e030723_003 e030723_008 e030723_01   e030923_   e041223_ e041823_003 e041823_005   e050223_ e122822_001 e122822_004 e123022_001 e123022_006 e123022_01 s221230_01 s221230_004 s221230_006 s221230_007 s122822_000 s122822_001 s122822_006 s122822_009
# 1   -0.5033140    2.470772 -0.85825470    4.481100  -9.111440  1.6456535  0.4536148    1.624038 -0.12518090 -0.1504577   2.6384374   -4.869520  0.47190348   -4.737946  -2.397738 -1.2268066  -1.2727475    6.326751  0.37821452  -0.3135986    4.393389   -1.684262    6.432359
# 2   -0.2659550    3.510277 -0.59783801    2.280014  -9.356976  0.4524719 -0.9451922    4.532692 -0.20670264 -0.1876598   1.6562945   -4.557020 -0.06786866   -5.303422  -2.427500 -1.9219293  -2.3209600    6.501137  0.41890462  -1.6399395    4.709464   -2.922916    6.678847
# 3    0.6947835    3.662865  0.44382932    5.024579  -7.861439 -0.5702554  0.2155228    3.763462 -0.32445638 -0.5224813   1.2098655   -3.667597 -2.02809673   -4.708184  -1.058453 -1.1589898  -3.3824410    6.658085  1.49719238  -2.8723802    5.101833   -2.527983    5.775057
# 4    1.1695014    4.273216  0.02716147    1.790884 -10.517690  0.1115628  0.7512345    3.547115 -0.01648554 -0.4034336   1.0908179   -3.763751  0.92644897   -3.606994  -3.588215 -0.2265082  -3.7804963    7.242278  1.95495605  -2.9545429    4.415187   -1.432940    5.505094
# 5    1.5764025    4.912178  0.36570374    3.584361 -10.897154  2.7251992 -0.2606684    2.537500 -0.11612333 -0.7159338   1.1801036   -3.499328  0.50031257   -5.749851  -1.713214 -1.9388835  -2.5730617    7.957261  2.09737142  -3.8583327    3.009199   -1.594504    5.728107
# 6    1.9041838    4.988472 -0.07700463    4.345231  -9.758761  0.5945171  0.6917114    2.585577 -0.13423907 -0.3141480   1.9836751   -3.307020  0.47190348   -6.166517  -3.498929 -1.4302572  -0.9543032    7.747998  1.89392090  -4.4686843    2.562334   -2.258710    4.929955
# 7    1.3616491    4.072945 -1.43117141    2.660448  -8.620368  0.3956534  0.3643315    2.681731 -0.06177540  0.3034115   0.2872468   -3.162788 -0.38036932   -3.696280  -3.469167 -0.7012261  -0.1714610    8.053174 -0.04903158  -3.4944693    3.150888   -1.774019    4.953430
# 8    1.2938323    2.661507 -0.62387968    3.665883  -8.040011  1.1342899 -0.7368575    3.354808 -0.53278993  0.1843637   0.5551039   -3.523365  0.47190285   -4.559375  -3.677500 -2.4136013   0.4388905    8.340911 -0.22196452  -2.8254301    3.968323   -2.115098    4.530879
# 9    0.7286919    1.822273  0.10528647    2.035448  -9.915011  1.7024720 -0.3797166    4.388462  0.07409428  0.3182924   0.8229614   -2.129135 -1.14741547   -4.946280  -2.070357 -2.8544108   0.2531314    7.582331  0.83597819  -3.8700702    5.036438   -2.007389    3.720990
# 10   1.1921070    1.898567  0.52195317    1.410448  -8.352511  0.8786080 -0.2011427    3.667308 -0.17047116 -0.2323028   1.0908182   -2.393558 -1.03377783   -5.481994  -1.177500 -1.8710666   0.1204462    8.593771  1.12080892  -4.1869835    5.472403   -2.258710    3.767940

options(width = old_width)

# ==============================================
# METADATA FOR ANAYSIS
# ==============================================
n <- 100
stimulation_time <- 100
baseline <- 100
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
# ii=16; analyse_PSC(response=OHDA_data[[ii]], dt=dt_range[ii], n=n, func=func, method=method, weight_method=weight_method, rel.decay.fit.limit=rel.decay.fit.limit, stimulation_time=stimulation_time, baseline=baseline, fast.constraint=TRUE, fast.decay.limit=fast.decay.limit, interval=interval, downsample=downsample, return.output=FALSE) 

out1  <- analyse_PSC(response=OHDA_data[[1]], dt=dt_range[1], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=910, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out1a <- analyse_PSC(response=OHDA_data[[1]], dt=dt_range[1], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=910, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, latency.limit=25)

out2  <- analyse_PSC(response=OHDA_data[[2]], dt=dt_range[2], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1040, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out2a <- analyse_PSC(response=OHDA_data[[2]], dt=dt_range[2], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=1040, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, latency.limit=25)

out3  <- analyse_PSC(response=OHDA_data[[3]], dt=dt_range[3], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=480, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=20, fast.constraint=TRUE)

out3a <- analyse_PSC(response=OHDA_data[[3]], dt=dt_range[3], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=480, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=20)

out4  <- analyse_PSC(response=OHDA_data[[4]], dt=dt_range[4], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=20, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out4a  <- analyse_PSC(response=OHDA_data[[4]], dt=dt_range[4], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out5  <- analyse_PSC(response=OHDA_data[[5]], dt=dt_range[5], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=270, stimulation_time=stimulation_time+20, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out5a  <- analyse_PSC(response=OHDA_data[[5]], dt=dt_range[5], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=270, stimulation_time=stimulation_time+20, baseline=baseline, downsample=downsample, latency.limit=15)

out6  <- analyse_PSC(response=OHDA_data[[6]], dt=dt_range[6], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out6a  <- analyse_PSC(response=OHDA_data[[6]], dt=dt_range[6], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=450, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out7  <- analyse_PSC(response=OHDA_data[[7]], dt=dt_range[7], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=20, fast.constraint=TRUE)

out7a  <- analyse_PSC(response=OHDA_data[[7]], dt=dt_range[7], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=20)

out8  <- analyse_PSC(response=OHDA_data[[8]], dt=dt_range[8], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out8a  <- analyse_PSC(response=OHDA_data[[8]], dt=dt_range[8], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out9  <- analyse_PSC(response=OHDA_data[[9]], dt=dt_range[9], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=225, stimulation_time=stimulation_time+5, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out9a  <- analyse_PSC(response=OHDA_data[[9]], dt=dt_range[9], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=225, stimulation_time=stimulation_time+5, baseline=baseline, downsample=downsample, latency.limit=15)

out10  <- analyse_PSC(response=OHDA_data[[10]], dt=dt_range[10], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=255, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out10a  <- analyse_PSC(response=OHDA_data[[10]], dt=dt_range[10], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=255, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out11  <- analyse_PSC(response=OHDA_data[[11]], dt=dt_range[11], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out11a  <- analyse_PSC(response=OHDA_data[[11]], dt=dt_range[11], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out12  <- analyse_PSC(response=OHDA_data[[12]], dt=dt_range[12], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time+20, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out12a  <- analyse_PSC(response=OHDA_data[[12]], dt=dt_range[12], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=350, stimulation_time=stimulation_time+20, baseline=baseline, downsample=downsample, latency.limit=15)

out13  <- analyse_PSC(response=OHDA_data[[13]], dt=dt_range[13], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out13a  <- analyse_PSC(response=OHDA_data[[13]], dt=dt_range[13], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out14  <- analyse_PSC(response=OHDA_data[[14]], dt=dt_range[14], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=1000,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=10, fast.decay.limit=c(50, 1000), latency.limit=25, fast.constraint=TRUE)

out14a  <- analyse_PSC(response=OHDA_data[[14]], dt=dt_range[14], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=1000,
  MLEsettings=MLEsettings, fit.limits=300, stimulation_time=stimulation_time, baseline=baseline, downsample=10, latency.limit=25)

out15  <- analyse_PSC(response=OHDA_data[[15]], dt=dt_range[15], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=15, fast.constraint=TRUE)

out15a  <- analyse_PSC(response=OHDA_data[[15]], dt=dt_range[15], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=15)

out16  <- analyse_PSC(response=OHDA_data[[16]], dt=dt_range[16], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out16a  <- analyse_PSC(response=OHDA_data[[16]], dt=dt_range[16], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=500, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out17  <- analyse_PSC(response=OHDA_data[[17]], dt=dt_range[17], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=730, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out17a  <- analyse_PSC(response=OHDA_data[[17]], dt=dt_range[17], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=730, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out18  <- analyse_PSC(response=OHDA_data[[18]], dt=dt_range[18], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=360, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out18a  <- analyse_PSC(response=OHDA_data[[18]], dt=dt_range[18], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=360, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out19  <- analyse_PSC(response=OHDA_data[[19]], dt=dt_range[19], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=520, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out19a  <- analyse_PSC(response=OHDA_data[[19]], dt=dt_range[19], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=520, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out20  <- analyse_PSC(response=OHDA_data[[20]], dt=dt_range[20], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=480, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out20a  <- analyse_PSC(response=OHDA_data[[20]], dt=dt_range[20], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=480, stimulation_time=stimulation_time, baseline=baseline, downsample=downsample, latency.limit=25)

out21  <- analyse_PSC(response=OHDA_data[[21]], dt=dt_range[21], n=n, func=func, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=1000,
  MLEsettings=MLEsettings, fit.limits=560, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out21a  <- analyse_PSC(response=OHDA_data[[21]], dt=dt_range[21], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval, half_width_fit_limit=1000,
  MLEsettings=MLEsettings, fit.limits=560, stimulation_time=stimulation_time+8, baseline=baseline, downsample=downsample, latency.limit=25)

out22  <- analyse_PSC(response=OHDA_data[[22]], dt=dt_range[22], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time+16, baseline=baseline, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out22a  <- analyse_PSC(response=OHDA_data[[22]], dt=dt_range[22], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=400, stimulation_time=stimulation_time+16, baseline=baseline, downsample=downsample)

out23  <- analyse_PSC(response=OHDA_data[[23]], dt=dt_range[23], n=n, func=func, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=250, stimulation_time=stimulation_time, baseline=100, downsample=downsample, fast.decay.limit=c(50, 500), latency.limit=25, fast.constraint=TRUE)

out23a  <- analyse_PSC(response=OHDA_data[[23]], dt=dt_range[23], n=n, func=product1N, method=method, weight_method=weight_method, interval=interval,
  MLEsettings=MLEsettings, fit.limits=250, stimulation_time=stimulation_time, baseline=100, downsample=downsample, latency.limit=25)

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
out23$BIC < out23a$BIC

# organise outputs
OHDA_summary <- list(
  out1,  out2,  out3a,  out4a,  out5,  out6,  out7,  out8,
  out9,  out10, out11a, out12, out13, out14, out15, out16,
  out17a, out18, out19, out20, out21a, out22, out23
)

names(OHDA_summary) <- 1:length(OHDA_summary)

OHDA_fits <- t(sapply(1:length(OHDA_summary), function(ii){
  X <- OHDA_summary[[ii]]$output
  if (dim(X)[1] == 1){
    X <- if (X[3] > 50) c(rep(NA, dim(X)[2]), as.vector(t(X))) else c(as.vector(t(X)), rep(NA, dim(X)[2]))
  }
  as.vector(t(X))
  })
)
# Create new column names by appending 1 and 2 to the original names
new_colnames <- rep(colnames(OHDA_summary[[1]]$output),2)
new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')

colnames(OHDA_fits) <- new_colnames
rownames(OHDA_fits) <- 1:length(OHDA_summary)

cols <- which(colnames(OHDA_fits) %in% c("A1", "area1"))
OHDA_fits[, cols][is.na(OHDA_fits[, cols])] <- 0


stimulation_range <- c(
  105, 105, 100, 100, 120, 100, 100, 100, 105, 100,
  100, 120, 100, 100, 100, 100, 100, 100, 100, 100,
  105, 110, 100
)


OHDA_peaks <- sapply(1:length(OHDA_data), function(ii){
  y <- OHDA_data[,ii]
  x = seq(y)*dt_range[ii] - dt_range[ii]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  peak.fun(y=y, dt=dt_range[ii], stimulation_time=stimulation_range[ii], baseline=baseline)
  })

OHDA_areas <- sapply(1:length(OHDA_data), function(ii){
  y <- OHDA_data[,ii]
  x = seq(y)*dt_range[ii] - dt_range[ii]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]

  keep <- x >= stimulation_range[ii]
  x <- x[keep] - stimulation_range[ii]
  y <- y[keep]

  charge_transfer_fun(y=y, x=x, tmax=600, dt=dt_range[ii], baseline=baseline)
  })


OHDA_fits
# using n <- 100 and method <- 'BF.LM':

#          A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -40.324 44.139 45.554 44.839 17.576  97.317 18.701    109.703  4915.577  -92.461 30.138 302.045 77.160 27.376 420.366  0.000    302.953 36055.783
# 2   -84.684  3.479 27.945  8.280  2.991  39.049 18.425     29.669  3182.557 -212.180 14.423 397.579 49.635 15.882 551.176 11.409    331.375 95575.627
# 3     0.000     NA     NA     NA     NA      NA     NA         NA     0.000  -75.259 23.130  95.886 43.349 16.359 138.920  3.220    125.692 11341.010
# 4    -5.501  2.896  2.897  2.896  1.135   6.285  4.650      7.086    43.313    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 5  -831.364  9.203  9.250  9.226  3.617  20.020 11.758     22.571 20849.887 -506.381 20.359  56.369 32.456 12.478  86.461  7.131     86.444 50766.120
# 6   -44.356  3.731  3.733  3.732  1.463   8.099  4.863      9.131   450.016  -46.596 16.404 131.338 38.994 14.092 183.540  9.200    139.559  8235.392
# 7   -79.585  3.578 10.546  5.853  2.245  15.994  4.978     15.764  1462.065 -136.464 25.866  80.136 43.190 16.529 120.522 18.610    117.413 18746.101
# 8  -111.791  3.809  3.810  3.809  1.493   8.266  9.433      9.319  1157.552 -116.479 44.385  44.435 44.410 17.408  96.366  1.004    108.644 14061.301
# 9   -25.520  5.794  5.835  5.815  2.279  12.617  0.000     14.225   403.364  -41.088 17.105  62.083 30.436 11.564  91.289  5.918     85.512  4164.856
# 10  -18.417  6.678  6.725  6.701  2.627  14.542  4.433     16.394   335.484  -14.573 27.120  47.375 35.384 13.790  81.715  0.800     88.753  1456.979
# 11  -67.863 46.560 46.570 46.565 18.253 101.042  1.300    113.916  8589.927    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 12 -133.463  2.580 18.606  5.918  2.156  26.076 12.574     20.358  3413.101  -30.936 16.900 151.819 41.749 14.949 211.651  0.000    156.484  6183.354
# 13  -73.778  3.159  3.159  3.159  1.238   6.855  5.184      7.729   633.600 -105.299 11.634  79.918 26.239  9.592 112.183 10.317     88.689 11685.812
# 14  -28.310  9.000 30.132 15.507  5.915  44.788  5.121     42.817  1427.167  -22.239  7.095 177.310 23.788  7.698 245.816 21.671    149.780  4509.275
# 15  -33.940  4.944  4.953  4.949  1.940  10.738  4.505     12.106   456.553  -19.233 13.968 157.468 37.130 13.038 218.861 11.000    153.650  3833.854
# 16   -2.862  0.062 42.096  0.405  0.084  58.357  0.073     29.603   121.639  -17.141  3.828 468.389 18.552  4.888 649.326  9.156    344.540  8352.819
# 17    0.000     NA     NA     NA     NA      NA     NA         NA     0.000 -109.537  9.561 194.917 30.313 10.035 270.261 21.749    169.772 24943.173
# 18  -49.564  5.469 20.191  9.797  3.720  29.628  7.016     27.633  1625.776 -112.524  4.090  63.897 12.011  4.088  88.641 24.185     58.283  8676.800
# 19  -50.950  5.148  5.153  5.151  2.019  11.176  7.869     12.600   713.342 -266.232 10.961 112.875 28.309 10.020 157.038 19.423    112.479 38617.205
# 20  -66.809  1.588 11.944  3.695  1.342  16.719 23.830     12.913  1087.303  -41.555 23.958 132.212 49.979 18.546 187.424  7.033    157.022  8018.106
# 21    0.000     NA     NA     NA     NA      NA     NA         NA     0.000  -15.893  6.380 505.789 28.256  7.902 701.173  0.000    381.152  8500.350
# 22 -429.180  5.622 45.166 13.379  4.834  63.113  4.140     47.948 26067.565 -427.582 18.291 101.405 38.222 14.179 143.714  9.653    120.254 63208.194
# 23   -8.659  3.030  3.035  3.032  1.189   6.580 12.376      7.419    71.382  -34.446 14.560  65.250 28.113 10.564  93.859  6.320     83.174  3458.102

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-OHDA_fits[,1], -OHDA_fits[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=200
setwd(svg_path)
DBSCAN_analyse(data=data, eps=200, filename=paste0('dbscan_', name, '.svg'), save=TRUE)
setwd(analysis_path)

# 5 & 22 outliers

# List of objects to keep
keep_objects <- c('OHDA_peaks', 'OHDA_areas', 'OHDA_fits', 
  'OHDA_summary', 'OHDA_data', 'name', 'analysis_path')

# Remove all objects except the ones in the keep_objects list
rm(list = setdiff(ls(), keep_objects))

# Save environment
save.image(file = file.path(analysis_path, paste0(name, '.RData')))

# ==============================================
# RELOAD FITS FROM '.RDATA'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure S5'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path


name <- '6OHDA'

load(paste0(analysis_path, '/', name, '.RData'))

OHDA_fits

#          A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -40.324 44.139 45.554 44.839 17.576  97.317 18.701    109.703  4915.577  -92.461 30.138 302.045 77.160 27.376 420.366  0.000    302.953 36055.783
# 2   -84.684  3.479 27.945  8.280  2.991  39.049 18.425     29.669  3182.557 -212.180 14.423 397.579 49.635 15.882 551.176 11.409    331.375 95575.627
# 3     0.000     NA     NA     NA     NA      NA     NA         NA     0.000  -75.259 23.130  95.886 43.349 16.359 138.920  3.220    125.692 11341.010
# 4    -5.501  2.896  2.897  2.896  1.135   6.285  4.650      7.086    43.313    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 5  -831.364  9.203  9.250  9.226  3.617  20.020 11.758     22.571 20849.887 -506.381 20.359  56.369 32.456 12.478  86.461  7.131     86.444 50766.120
# 6   -44.356  3.731  3.733  3.732  1.463   8.099  4.863      9.131   450.016  -46.596 16.404 131.338 38.994 14.092 183.540  9.200    139.559  8235.392
# 7   -79.585  3.578 10.546  5.853  2.245  15.994  4.978     15.764  1462.065 -136.464 25.866  80.136 43.190 16.529 120.522 18.610    117.413 18746.101
# 8  -111.791  3.809  3.810  3.809  1.493   8.266  9.433      9.319  1157.552 -116.479 44.385  44.435 44.410 17.408  96.366  1.004    108.644 14061.301
# 9   -25.520  5.794  5.835  5.815  2.279  12.617  0.000     14.225   403.364  -41.088 17.105  62.083 30.436 11.564  91.289  5.918     85.512  4164.856
# 10  -18.417  6.678  6.725  6.701  2.627  14.542  4.433     16.394   335.484  -14.573 27.120  47.375 35.384 13.790  81.715  0.800     88.753  1456.979
# 11  -67.863 46.560 46.570 46.565 18.253 101.042  1.300    113.916  8589.927    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 12 -133.463  2.580 18.606  5.918  2.156  26.076 12.574     20.358  3413.101  -30.936 16.900 151.819 41.749 14.949 211.651  0.000    156.484  6183.354
# 13  -73.778  3.159  3.159  3.159  1.238   6.855  5.184      7.729   633.600 -105.299 11.634  79.918 26.239  9.592 112.183 10.317     88.689 11685.812
# 14  -28.310  9.000 30.132 15.507  5.915  44.788  5.121     42.817  1427.167  -22.239  7.095 177.310 23.788  7.698 245.816 21.671    149.780  4509.275
# 15  -33.940  4.944  4.953  4.949  1.940  10.738  4.505     12.106   456.553  -19.233 13.968 157.468 37.130 13.038 218.861 11.000    153.650  3833.854
# 16   -2.862  0.062 42.096  0.405  0.084  58.357  0.073     29.603   121.639  -17.141  3.828 468.389 18.552  4.888 649.326  9.156    344.540  8352.819
# 17    0.000     NA     NA     NA     NA      NA     NA         NA     0.000 -109.537  9.561 194.917 30.313 10.035 270.261 21.749    169.772 24943.173
# 18  -49.564  5.469 20.191  9.797  3.720  29.628  7.016     27.633  1625.776 -112.524  4.090  63.897 12.011  4.088  88.641 24.185     58.283  8676.800
# 19  -50.950  5.148  5.153  5.151  2.019  11.176  7.869     12.600   713.342 -266.232 10.961 112.875 28.309 10.020 157.038 19.423    112.479 38617.205
# 20  -66.809  1.588 11.944  3.695  1.342  16.719 23.830     12.913  1087.303  -41.555 23.958 132.212 49.979 18.546 187.424  7.033    157.022  8018.106
# 21    0.000     NA     NA     NA     NA      NA     NA         NA     0.000  -15.893  6.380 505.789 28.256  7.902 701.173  0.000    381.152  8500.350
# 22 -429.180  5.622 45.166 13.379  4.834  63.113  4.140     47.948 26067.565 -427.582 18.291 101.405 38.222 14.179 143.714  9.653    120.254 63208.194
# 23   -8.659  3.030  3.035  3.032  1.189   6.580 12.376      7.419    71.382  -34.446 14.560  65.250 28.113 10.564  93.859  6.320     83.174  3458.102
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
#   write.csv(OHDA_fits, file = 'OHDA_fits.csv', row.names = FALSE)
#   write.csv(OHDA_peaks, file = 'OHDA_peaks.csv', row.names = FALSE)
#   write.csv(OHDA_areas, file = 'OHDA_areas.csv', row.names = FALSE)
#   }

setwd(svg_path)
for (ii in 1:length(OHDA_summary)){
  fit_plot(traces=OHDA_summary[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), 
    xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=save)
}


