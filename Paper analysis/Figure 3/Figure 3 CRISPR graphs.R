# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- FALSE

UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 3'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
# load RDATA
datasets <- c('CRISPR control',    # ctrl
              'CRISPR delta KD')   # test

datasets2 <- gsub(" ", "_", datasets)

# load *.Rdata datasets
loaded_objects <- list()

# Loop through each dataset
for (ii in 1:length(datasets)) {
  # Get the list of objects in the environment before loading
  objects_before <- ls()
  
  # Load the .RData file
  load(file = file.path(analysis_path, paste0(datasets[ii], '.RData')))
  
  # Get the list of objects in the environment after loading
  objects_after <- ls()
  
  # Identify the new objects that were added
  new_objects <- setdiff(objects_after, objects_before)
  
  # Store the names of the new objects in the list
  loaded_objects[[datasets2[ii]]] <- new_objects
}

fits_list <- datasets2list(datasets2, id='_fits')
peaks_list <- datasets2list(datasets2, id='_peaks')
areas_list <- datasets2list(datasets2, id='_areas')
data_list <- datasets2list(datasets2, id='_data')
output_list <- datasets2list(datasets2, id='_summary')

# Remove the loaded objects
remove_loaded_objects(loaded_objects)

# Display the list of fits objects
fits_list
# $CRISPR_control
#         A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -141.120  1.566  5.357  2.722  1.037  7.939 10.140      7.549  1256.568 -211.089 13.003 100.599 30.553 11.071 140.705  8.138    107.909 28771.24
# 2  -78.429  1.852  5.747  3.094  1.184  8.641 10.657      8.414   772.207  -99.447 15.454 100.255 34.163 12.540 141.031  5.007    113.173 14018.05
# 3 -143.396  2.038  4.137  2.844  1.104  6.821  8.973      7.245  1179.683 -298.921 12.815  83.320 28.351 10.405 117.197  7.356     93.993 35000.99
# 4 -208.987  1.297  4.164  2.198  0.840  6.228 11.636      6.016  1475.159 -241.229 14.317  70.111 28.581 10.685 100.151  4.613     86.655 25424.88
# 5  -74.823  1.993  1.994  1.993  0.781  4.325  8.224      4.877   405.432 -143.314  9.403  85.980 23.366  8.354 119.827  2.946     88.223 16169.95
# 6 -296.552 22.098 22.139 22.119  8.670 47.995  7.403     54.111 17830.022 -257.915 78.696  78.783 78.739 30.865 170.858 23.757    192.626 55202.95
# 7  -99.433  2.688  2.690  2.689  1.054  5.834 11.060      6.578   726.744 -176.352 11.752 105.063 28.985 10.383 146.481  7.180    108.425 24414.29

# $CRISPR_delta_KD
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

setwd(svg_path)

############################################################## AMPLITUDE ############################################################## 

amplitude_ctrl <- create_df(matrix(peaks_list$CRISPR_control, ncol=1), levels = list(drug = 'AP5+NBQX+CGP55845A', condition = 'ctrl'))
amplitude_test <- create_df(matrix(peaks_list$CRISPR_delta_KD,  ncol=1), levels = list(drug = 'AP5+NBQX+CGP55845A', condition = 'test'), 
  start_id = length(unique(amplitude_ctrl$s)) + 1)
amplitude_SPN <- rbind(amplitude_ctrl, amplitude_test)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ condition, df=amplitude_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
log_y <- FALSE
yrange <- if (log_y) c(10, 500) else c(-500, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ condition, data=amplitude_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=log_y, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(10, 1000) else c(-500, 0)
BoxPlot3(formula=amplitude ~ condition, data=transform(amplitude_SPN[,3:4], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=log_y, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_semilog_SPN.svg', width=width, height=height, bg='transparent')

########################################################### CHARGE TRANSFER ###########################################################

area_ctrl <- create_df(matrix(areas_list$CRISPR_control, ncol=1), levels = list(drug = 'AP5+NBQX+CGP55845A', condition = 'ctrl'), var_name = 'charge_transfer')
area_test <- create_df(matrix(areas_list$CRISPR_delta_KD,  ncol=1), levels = list(drug = 'AP5+NBQX+CGP55845A', condition = 'test'), var_name = 'charge_transfer',
  start_id = length(unique(area_ctrl$s)) + 1)
area_SPN <- rbind(area_ctrl, area_test)
area_SPN$charge_transfer <- -area_SPN$charge_transfer

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ condition, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- FALSE
yrange <- if (log_y) c(1, 100) else c(0, 100)
y_tick_interval <- 25

BoxPlot3(formula=charge_transfer ~ condition, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 100) else c(0, 100)

BoxPlot3(formula=charge_transfer ~ condition, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=log_y, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_semilog_SPN.svg', width=width, height=height, bg='transparent')

######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_ctrl <- create2condition_df(fits_list$CRISPR_control, levels = list(kinetics = c('fast', 'slow'), condition = 'ctrl'))
amplitude_2component_test <- create2condition_df(fits_list$CRISPR_delta_KD,  levels = list(kinetics = c('fast', 'slow'), condition = 'test'),
  start_id = dim(fits_list$CRISPR_control)[1] + 1)
amplitude_2component_SPN <- rbind(amplitude_2component_ctrl, amplitude_2component_test)

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ kinetics*condition + Error(s), df=amplitude_2component_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
log_y <- FALSE
yrange <- if (log_y) c(10, 350) else c(-350, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 50
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ kinetics*condition, data=amplitude_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(10, 350) else c(-350, 0)

BoxPlot3(formula=amplitude ~ kinetics*condition, data=transform(amplitude_2component_SPN, amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_semilog_2component_SPN.svg', width=width, height=height, bg='transparent')

# fast and slow plots

log_y <- FALSE
yrange <- if (log_y) c(10, 350) else c(-350, 0)

BoxPlot3(formula=amplitude ~ condition, data=amplitude_2component_SPN[amplitude_2component_SPN$kinetics == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary1[1,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_fast_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ condition, data=amplitude_2component_SPN[amplitude_2component_SPN$kinetics == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary1[2,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_slow_SPN.svg', width=width, height=height, bg='transparent')




# scatter plots
df <- amplitude_2component_SPN
df1 <- df[df$condition=='ctrl',]
df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

df1$x <- ifelse(df1$kinetics == 'fast', 1, 2)
df1$condition <- NULL
names(df1)[names(df1) == 'amplitude'] <- 'y'
df1 <- df1[, c('s', 'x', 'y')]
df1$s <- as.integer(as.character(df1$s))

width=2.875
height=3.5
xlim=c(0, 300)
ylim=c(0, 300)
ScatterPlot(df1, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter_ctrl.svg', lwd=lwd, save=plotsave)
ScatterPlot(df1, sign=-1, xlim=c(10, 300), ylim=c(10, 300), height=height, width=width, filename='scatter_logxy_ctrl.svg', lwd=lwd, log_xy=TRUE, save=plotsave)


# scatter plots
df2 <- df[df$condition=='test',]
df2 <- df2[order(df2$kinetics != 'fast', df2$s), ]

df2$x <- ifelse(df2$kinetics == 'fast', 1, 2)
df2$condition <- NULL
names(df2)[names(df2) == 'amplitude'] <- 'y'
df2 <- df2[, c('s', 'x', 'y')]
df2$s <- as.integer(as.character(df2$s))

ScatterPlot(df2, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter_test.svg', lwd=lwd, save=plotsave, open_symbols=TRUE)
ScatterPlot(df2, sign=-1, xlim=c(10, 300), ylim=c(10, 300), height=height, width=width, filename='scatter_logxy_test.svg', lwd=lwd, log_xy=TRUE, save=plotsave, open_symbols=TRUE)

##################################################### CHARGE TRANSFER FAST vs SLOW ####################################################
area_2component_ctrl <- create2condition_df(fits_list$CRISPR_control,  cols = c(9, 18), levels = list(kinetics = c('fast', 'slow'),   condition = 'ctrl'), var_name='charge_transfer')
area_2component_test <- create2condition_df(fits_list$CRISPR_delta_KD, cols = c(9, 18), levels = list(kinetics = c('fast', 'slow'), condition = 'test'), var_name='charge_transfer',
  start_id = dim(fits_list$CRISPR_control)[1] + 1)
area_2component_SPN <- rbind(area_2component_ctrl, area_2component_test)
area_2component_SPN$charge_transfer <- area_2component_SPN$charge_transfer/1e3

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=charge_transfer ~ kinetics*condition + Error(s), df=area_2component_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- FALSE
yrange <- if (log_y) c(0.1, 100) else c(0, 100)
y_tick_interval <- 20
width <- 3
height <- 3.5

# # if log_ y needs to omit any zeros
# if (log_y){
#   area_2component_SPN$charge_transfer[area_2component_SPN$charge_transfer == 0] <- NA
# }

BoxPlot3(formula=charge_transfer ~ kinetics*condition, data=area_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(0.1, 100) else c(0, 100)
BoxPlot3(formula=charge_transfer ~ kinetics*condition, data=area_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,], log_y=TRUE)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_semilog_2component_SPN.svg', width=width, height=height, bg='transparent')

# scatter plots
df <- area_2component_SPN
df1 <- df[df$condition=='ctrl',]
df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

df1$x <- ifelse(df1$kinetics == 'fast', 1, 2)
df1$condition <- NULL
names(df1)[names(df1) == 'charge_transfer'] <- 'y'
df1 <- df1[, c('s', 'x', 'y')]
df1$s <- as.integer(as.character(df1$s))

width <- 2.875
height <- 3.5
xlim <- c(0, 60)
ylim <- c(0, 60)
x_tick_interval <- 20 
y_tick_interval <- 20
ScatterPlot(df1, xlim=xlim, ylim=ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter_charge_transfer_ctrl.svg', lwd=lwd, save=plotsave)
ScatterPlot(df1, xlim=c(0.1, 60), ylim=c(0.1, 60), x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter_charge_transfer_logxy_ctrl.svg', lwd=lwd, log_xy=TRUE, save=plotsave)

# scatter plots
df <- area_2component_SPN
df1 <- df[df$condition=='test',]
df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

df1$x <- ifelse(df1$kinetics == 'fast', 1, 2)
df1$condition <- NULL
names(df1)[names(df1) == 'charge_transfer'] <- 'y'
df1 <- df1[, c('s', 'x', 'y')]
df1$s <- as.integer(as.character(df1$s))

width <- 2.875
height <- 3.5
xlim <- c(0, 60)
ylim <- c(0, 60)
x_tick_interval <- 20 
y_tick_interval <- 20
ScatterPlot(df1, xlim=xlim, ylim=ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter_charge_transfer_ctrl.svg', lwd=lwd, open_symbols=TRUE, save=plotsave)
ScatterPlot(df1, xlim=c(0.1, 60), ylim=c(0.1, 60), x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter_charge_transfer_logxy_ctrl.svg', lwd=lwd, open_symbols=TRUE, log_xy=TRUE, save=plotsave)

################################################### DECAY TIME CONSTANT FAST vs SLOW ##################################################

tau_2component_ctrl <- create2condition_df(fits_list$CRISPR_control, cols = c(3, 12), levels = list(kinetics = c('fast', 'slow'), condition = 'ctrl'), var_name='tau')
tau_2component_test <- create2condition_df(fits_list$CRISPR_delta_KD,  cols = c(3, 12), levels = list(kinetics = c('fast', 'slow'), condition = 'test'), 
  start_id = dim(fits_list$CRISPR_control)[1] + 1, var_name='tau')
tau_2component_SPN <- rbind(tau_2component_ctrl, tau_2component_test)

# create output for stats
stats_summary <-  rbind(stats_summary,  MCwilcox(formula=tau ~ kinetics*condition + Error(s), df=tau_2component_SPN)) # remove paired as one tau is NA because 1product fit was better
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ kinetics*condition + Error(s), df=tau_2component_SPN))

# update graph properties
log_y <- FALSE
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- if (log_y) c(1, 125) else c(0, 125)
y_tick_interval <- 25
width <- 3
height <- 3.5

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_ctrl.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 200) else c(0, 200)
BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_semilog_ctrl.svg', width=width, height=height, bg='transparent')

ylab <- expression(tau[decay * ',' * fast] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 30) else c(0, 30)
y_tick_interval <- 5

BoxPlot3(formula=tau ~ condition, data=tau_2component_SPN[tau_2component_SPN$kinetics == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[3,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 30) else c(0, 30)
BoxPlot3(formula=tau ~ condition, data=tau_2component_SPN[tau_2component_SPN$kinetics == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[3,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_semilog_SPN.svg', width=width, height=height, bg='transparent')


ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 125) else c(0, 125)
y_tick_interval <- 25

BoxPlot3(formula=tau ~ condition, data=tau_2component_SPN[tau_2component_SPN$kinetics == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[4,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 125) else c(0, 125)
BoxPlot3(formula=tau ~ condition, data=tau_2component_SPN[tau_2component_SPN$kinetics == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[5,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_semilog_SPN.svg', width=width, height=height, bg='transparent')

########################################################## SINGLE EXAMPLES ###########################################################
# single examples
name <- datasets[1]
for (ii in 1:length(output_list$CRISPR_control)){
  fit_plot(traces=output_list$CRISPR_control[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}

name <- datasets[2]
for (ii in 1:length(output_list$CRISPR_delta_KD)){
  traces <- output_list$CRISPR_delta_KD[[ii]]$traces
  func <- if (dim(traces)[2] == 4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}
graphics.off()

# single egs
dx <- 0.1
xlab='time (ms)'
ylab='PSC amplitude (pA)'
lwd=1.0
filter='off'
xlim <- c(0, 800)
ylim <- c(-300, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 7
ctrls_eg <- output_list$CRISPR_control[[ii]]$traces
single_fit_egs(traces=ctrls_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_ctrl.svg', save=plotsave)
single_fit_egs(traces=ctrls_eg, xlim=c(100, 300), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_ctrl.svg', log_y=TRUE, save=plotsave)

ii = 4
tests_eg <- output_list$CRISPR_delta_KD[[ii]]$traces
single_fit_egs(traces=tests_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_test.svg', save=plotsave)
single_fit_egs(traces=tests_eg, xlim=c(100, 300), ylim=NULL, lwd=lwd, colors=colors, 
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_test.svg', log_y=TRUE, save=plotsave)

######################################################### RENAME AND SAVE ##########################################################

# rename output correctly 
for (i in seq_along(fits_list)) {
  # Get the column names of the current element
  cols <- colnames(fits_list[[i]])
  
  # Rename the first A1 to Afast and the second A1 to Aslow
  first_A1 <- which(cols == 'A1')[1]  # Find the first occurrence of A1
  second_A1 <- which(cols == 'A1')[2]  # Find the second occurrence of A1

  if (!is.na(first_A1)) cols[first_A1] <- 'Afast'
  if (!is.na(second_A1)) cols[second_A1] <- 'Aslow'
  
  # Rename all occurrences of area1 to area
  cols <- gsub('^area1$', 'area', cols)
  
  # Update column names in the current element
  colnames(fits_list[[i]]) <- cols
}

fits_list
# $CRISPR_control
# #      Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width      area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area
# # 1 -141.120  1.566  5.357  2.722  1.037  7.939 10.140      7.549  1256.568 -211.089 13.003 100.599 30.553 11.071 140.705  8.138    107.909 28771.24
# # 2  -78.429  1.852  5.747  3.094  1.184  8.641 10.657      8.414   772.207  -99.447 15.454 100.255 34.163 12.540 141.031  5.007    113.173 14018.05
# # 3 -143.396  2.038  4.137  2.844  1.104  6.821  8.973      7.245  1179.683 -298.921 12.815  83.320 28.351 10.405 117.197  7.356     93.993 35000.99
# # 4 -208.987  1.297  4.164  2.198  0.840  6.228 11.636      6.016  1475.159 -241.229 14.317  70.111 28.581 10.685 100.151  4.613     86.655 25424.88
# # 5  -74.823  1.993  1.994  1.993  0.781  4.326  8.224      4.877   405.435 -143.314  9.403  85.980 23.366  8.354 119.827  2.946     88.223 16169.94
# # 6 -296.552 22.099 22.138 22.119  8.670 47.995  7.403     54.111 17830.019 -257.915 78.690  78.788 78.739 30.865 170.858 23.757    192.626 55202.95
# # 7  -99.433  2.688  2.690  2.689  1.054  5.834 11.060      6.578   726.744 -176.352 11.752 105.063 28.985 10.383 146.481  7.180    108.425 24414.29

# # $CRISPR_delta_KD
# #       Afast τrise τdecay tpeak r20_80 d80_20 delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area
# # 1   -89.478 4.162 30.244 9.572  3.485 42.377 4.538     33.021 3713.756  -46.720  8.136 106.397 22.649  7.845 147.713 23.169    100.525  6150.049
# # 2  -102.475 1.129 18.499 3.362  1.139 25.659 8.671     16.726 2273.531  -38.202  1.354 160.435  6.522  1.726 222.411  4.466    118.197  6383.152
# # 3  -111.502 1.021  9.396 2.542  0.909 13.093 8.859      9.625 1373.204 -129.494  1.099  58.502  4.452  1.313  81.102  4.437     45.422  8174.684
# # 4   -82.918 3.561  3.566 3.563  1.397  7.732 5.037      8.717  803.165 -117.639 11.302  88.266 26.641  9.646 123.425 10.692     94.433 14041.880
# # 5   -67.201 2.567 10.365 4.762  1.800 15.056 9.372     13.718 1102.728  -60.572  1.563  65.992  5.992  1.820  91.485  5.062     52.352  4377.214
# # 6  -330.918 1.483 10.825 3.416  1.243 15.165 9.869     11.803 4911.133 -170.532 13.154  65.507 26.424  9.869  93.463  4.316     80.502 16721.646
# # 7  -169.987 1.013 13.276 2.822  0.977 18.431 9.240     12.538 2791.281  -99.571 12.316  80.909 27.346 10.029 113.760  3.892     90.988 11295.804
# # 8  -107.180 5.780  5.796 5.788  2.269 12.559 4.318     14.159 1686.208  -23.843 10.163  36.275 17.964  6.831  53.454 20.461     50.285  1419.167
# # 9   -30.233 2.186  2.191 2.188  0.858  4.748 9.400      5.353  179.820 -233.229 14.865  93.223 32.470 11.947 131.321  4.964    106.328 30801.338
# # 10  -75.489 5.465 11.706 7.808  3.028 19.022 4.856     20.020 1721.790  -41.947  4.939 102.363 15.731  5.198 141.928  3.436     88.924  5007.153

# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'charge transfer 2 components' = area_2component_SPN,
        'tau 2 components' = tau_2component_SPN,
        'ctrl single examples' = ctrls_eg,
        'CRISPR KD single examples' = tests_eg,
        'statistics' = stats_summary
      )
    )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
}

