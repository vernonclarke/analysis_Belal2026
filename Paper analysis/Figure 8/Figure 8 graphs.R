# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- TRUE

UserName <- Sys.getenv('USER')
root_dir <- file.path('/Users', UserName, 'Documents', 'Repositories', 'analysis_Belal2026')

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 8'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
# load RDATA
datasets <- c('Control for MCI-Park',  # ctrl
              'ChAT-Flp X Ndufs2 fl-fl X DAT-Cre-MCI-PARK')   # test

datasets2 <- c('ctrl', 'test')

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
# $ctrl
#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -100.262 2.407  7.228 3.968  1.521 10.927 17.439     10.724 1254.885 -148.460 23.715 181.256 55.492 20.126 253.599 14.358    195.086 36547.86
# 2  -89.222 2.283 16.275 5.217  1.902 22.818 17.147     17.868 2000.751 -224.536 38.027 259.234 85.538 31.287 363.991 14.334    288.332 80961.67
# 3  -44.309 6.049  6.053 6.051  2.372 13.130 20.968     14.803  728.790  -98.728 20.298 109.626 42.012 15.612 155.599 24.370    131.118 15877.83
# 4  -71.638 4.324  9.477 6.241  2.418 15.310 15.765     16.048 1311.611 -144.206 23.122 237.756 59.689 21.131 330.783 12.739    237.006 44070.23
# 5  -30.310 6.983  7.024 7.004  2.745 15.197 22.551     17.134  577.041  -60.592 24.346 241.547 62.129 22.063 336.214 18.367    242.868 18928.71
# 6 -112.088 5.167  5.174 5.171  2.027 11.220 16.033     12.649 1575.402 -255.573 25.873 187.838 59.483 21.657 263.198 18.679    205.132 65891.19
# 7  -54.341 2.729  2.730 2.729  1.070  5.922 12.671      6.677  403.156  -88.457 12.766 164.094 35.350 12.264 227.840 15.420    155.596 18004.56
# 8  -42.510 9.455  9.536 9.495  3.722 20.604 13.302     23.229 1097.229 -108.880 18.795 198.672 48.951 17.288 276.316 23.355    196.756 27675.17

# $test
#       A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width    area1
# 1 -2.656 15.077 16.123 15.588  6.110 33.856 19.802     38.148 112.611  -7.196 18.706 118.182 40.967 15.065 166.429 42.435    134.490 1202.775
# 2 -4.075 11.410 11.427 11.419  4.476 24.778 17.958     27.935 126.479  -4.534 37.215 351.483 93.460 33.327 489.590  4.865    357.719 2079.268
# 3 -3.612  7.977  9.643  8.757  3.430 19.140 13.700     21.485  86.373 -14.271 10.174 248.706 33.908 11.001 344.800 24.525    210.750 4067.683
# 4 -9.975  9.982  9.993  9.988  3.915 21.672 24.042     24.433 270.813 -12.413 24.149 161.634 53.974 19.768 227.102 11.682    180.759 2801.849
# 5 -3.465  1.972  8.401  3.735  1.407 12.140  7.893     10.904  45.399  -6.858 17.669 168.749 44.536 15.866 235.014 23.021    171.268 1506.739
# 6 -1.963  1.935  1.938  1.937  0.759  4.203 10.173      4.738  10.333  -3.243 26.890  80.292 44.227 16.952 121.495 22.870    119.396  451.759

setwd(svg_path)

############################################################## AMPLITUDE ############################################################## 

amplitude_ctrl <- create_df(matrix(peaks_list$ctrl, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'ctrl'))
amplitude_test <- create_df(matrix(peaks_list$test,  ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'test'), 
  start_id = length(unique(amplitude_ctrl$s)) + 1)
amplitude_SPN <- rbind(amplitude_ctrl, amplitude_test)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ cell_type, df=amplitude_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-300, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot(formula=amplitude ~ cell_type, data=amplitude_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot(formula=amplitude ~ cell_type, data=transform(amplitude_SPN[,3:4], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,300),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_amplitude_SPN.svg', width=width, height=height, bg='transparent')

########################################################### CHARGE TRANSFER ###########################################################

area_ctrl <- create_df(matrix(areas_list$ctrl, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'ctrl'), var_name = 'charge_transfer')
area_test <- create_df(matrix(areas_list$test,  ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'test'), var_name = 'charge_transfer',
  start_id = length(unique(area_ctrl$s)) + 1)
area_SPN <- rbind(area_ctrl, area_test)
area_SPN$charge_transfer <- -area_SPN$charge_transfer

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 100)
y_tick_interval <- 25

BoxPlot(formula=charge_transfer ~ cell_type, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot(formula=charge_transfer ~ cell_type, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,100),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_ctrl <- create2condition_df(fits_list$ctrl, levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'))
amplitude_2component_test <- create2condition_df(fits_list$test,  levels = list(condition = c('fast', 'slow'), cell_type = 'test'),
  start_id = dim(fits_list$ctrl)[1] + 1)
amplitude_2component_SPN <- rbind(amplitude_2component_ctrl, amplitude_2component_test)

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ condition*cell_type + Error(s), df=amplitude_2component_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-300, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot(formula=amplitude ~ condition*cell_type + Error(s), data=amplitude_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,])

if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

# paired scatter plot Afast vs Aslow
wide_df <- reshape(
  amplitude_2component_SPN,
  idvar = c('s', 'cell_type'),
  timevar = 'condition',
  direction = 'wide'
)
names(wide_df)[names(wide_df) == 'amplitude.fast'] <- 'Afast'
names(wide_df)[names(wide_df) == 'amplitude.slow'] <- 'Aslow'
wide_df$cell_type <- ifelse(wide_df$cell_type == 'ctrl', 1, 2)
names(wide_df)[names(wide_df) == 'cell_type'] <- 'level'
wide_df$Afast <- -wide_df$Afast
wide_df$Aslow <- -wide_df$Aslow
names(wide_df)[names(wide_df) == 'Afast'] <- 'x'
names(wide_df)[names(wide_df) == 'Aslow'] <- 'y'

width=2.875
height=3.5
xlim=c(0, 300)
ylim=c(0, 300)
scatter_plot(wide_df, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter.svg', lwd=lwd, save=plotsave)

##################################################### CHARGE TRANSFER FAST vs SLOW ####################################################

area_2component_ctrl <- create2condition_df(fits_list$ctrl, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'), var_name='area')
area_2component_test <- create2condition_df(fits_list$test,  cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'test'),
  start_id = dim(fits_list$ctrl)[1] + 1, var_name='area')
area_2component_SPN <- rbind(area_2component_ctrl, area_2component_test)
area_2component_SPN$area <- area_2component_SPN$area/1e3

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=area ~ condition*cell_type + Error(s), df=area_2component_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- TRUE
yrange <- if (log_y) c(0.001, 200) else c(0, 100)
y_tick_interval <- 20
width <- 3
height <- 3.5

# if log_ y needs to omit any zeros

if (log_y){
  area_2component_SPN$area[area_2component_SPN$area == 0] <- NA
}

BoxPlot(formula=area ~ condition*cell_type + Error(s), data=area_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,], log_y=log_y)

if (plotsave) save_graph(svg_path=svg_path, filename='area_2component_SPN.svg', width=width, height=height, bg='transparent')

# paired scatter plot Areafast vs Areaslow
wide_df <- reshape(
  area_2component_SPN,
  idvar = c('s', 'cell_type'),
  timevar = 'condition',
  direction = 'wide'
)
names(wide_df)[names(wide_df) == 'area.fast'] <- 'Areafast'
names(wide_df)[names(wide_df) == 'area.slow'] <- 'Areaslow'
wide_df$cell_type <- ifelse(wide_df$cell_type == 'ctrl', 1, 2)
names(wide_df)[names(wide_df) == 'cell_type'] <- 'level'
names(wide_df)[names(wide_df) == 'Areafast'] <- 'x'
names(wide_df)[names(wide_df) == 'Areaslow'] <- 'y'

width=2.875
height=3.5
xlim=c(0, 5)
ylim=c(0, 60)
scatter_plot(wide_df, xlim = xlim, ylim = ylim, height = height, width = width, 
  xlab = expression(charge~transfer[fast]*' '(pC)), ylab = expression(charge~transfer[slow]*' '(pC)),
  x_tick_interval=1, y_tick_interval=10,
  filename = 'scatter_charge_transfer.svg', lwd = lwd, save = plotsave)

################################################### DECAY TIME CONSTANT FAST vs SLOW ##################################################

tau_2component_ctrl <- create2condition_df(fits_list$ctrl, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'), var_name='tau')
tau_2component_test <- create2condition_df(fits_list$test,  cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'test'),
  start_id = dim(fits_list$ctrl)[1] + 1, var_name='tau')
tau_2component_SPN <- rbind(tau_2component_ctrl, tau_2component_test)

# create output for stats
stats_summary <-  rbind(stats_summary,  MCwilcox(formula=tau ~ condition*cell_type, df=tau_2component_SPN)[1:2,]) # remove paired as one tau is NA because 1product fit was better
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ condition*cell_type , df=tau_2component_SPN))

# update graph properties
log_y <- TRUE
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- if (log_y) c(1, 500) else c(0, 500)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot(formula=tau ~ condition*cell_type + Error(s), data=tau_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_ctrl.svg', width=width, height=height, bg='transparent')

ylab <- expression(tau[decay * ',' * fast] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 30) else c(0, 30)
y_tick_interval <- 5

BoxPlot(formula=tau ~ cell_type, data=tau_2component_SPN[tau_2component_SPN$condition == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[3,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_SPN.svg', width=width, height=height, bg='transparent')


ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 500) else c(0, 500)
y_tick_interval <- 100

BoxPlot(formula=tau ~ cell_type, data=tau_2component_SPN[tau_2component_SPN$condition == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[4,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_SPN.svg', width=width, height=height, bg='transparent')


########################################################## SINGLE EXAMPLES ###########################################################
# single examples

name <- datasets[1]
for (ii in 1:length(output_list$ctrl)){
  fit_plot(traces=output_list$ctrl[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}

name <- datasets[2]
for (ii in 1:length(output_list$test)){
  traces <- output_list$test[[ii]]$traces
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
xlim <- c(0, 1000)
ylim <- c(-200, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 1
ctrls_eg <- output_list$ctrl[[ii]]$traces
single_fit_egs(traces=ctrls_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_ctrl.svg', save=plotsave)
single_fit_egs(traces=ctrls_eg, xlim=c(100, 500), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_ctrl.svg', log_y=TRUE, save=plotsave)

ii = 2
tests_eg <- output_list$test[[ii]]$traces
single_fit_egs(traces=tests_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_test.svg', save=plotsave)
single_fit_egs(traces=tests_eg, xlim=c(100, 500), ylim=NULL, lwd=lwd, colors=colors, 
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_test.svg', log_y=TRUE, save=plotsave)

######################################################### SIMULATED EXAMPLE ##########################################################

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
# $ctrl
#      Afast τrise τdecay tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area
# 1 -100.262 2.407  7.228 3.968  1.521 10.927 17.439     10.724 1254.885 -148.460 23.715 181.256 55.492 20.126 253.599 14.358    195.086 36547.86
# 2  -89.222 2.283 16.275 5.217  1.902 22.818 17.147     17.868 2000.751 -224.536 38.027 259.234 85.538 31.287 363.991 14.334    288.332 80961.67
# 3  -44.309 6.049  6.053 6.051  2.372 13.130 20.968     14.803  728.790  -98.728 20.298 109.626 42.012 15.612 155.599 24.370    131.118 15877.83
# 4  -71.638 4.324  9.477 6.241  2.418 15.310 15.765     16.048 1311.611 -144.206 23.122 237.756 59.689 21.131 330.783 12.739    237.006 44070.23
# 5  -30.310 6.983  7.024 7.004  2.745 15.197 22.551     17.134  577.041  -60.592 24.346 241.547 62.129 22.063 336.214 18.367    242.868 18928.71
# 6 -112.088 5.167  5.174 5.171  2.027 11.220 16.033     12.649 1575.402 -255.573 25.873 187.838 59.483 21.657 263.198 18.679    205.132 65891.19
# 7  -54.341 2.729  2.730 2.729  1.070  5.922 12.671      6.677  403.156  -88.457 12.766 164.094 35.350 12.264 227.840 15.420    155.596 18004.56
# 8  -42.510 9.455  9.536 9.495  3.722 20.604 13.302     23.229 1097.229 -108.880 18.795 198.672 48.951 17.288 276.316 23.355    196.756 27675.17

# $test
#    Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area   Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area
# 1 -2.656 15.077 16.123 15.588  6.110 33.856 19.802     38.148 112.611  -7.196 18.706 118.182 40.967 15.065 166.429 42.435    134.490 1202.775
# 2 -4.075 11.410 11.427 11.419  4.476 24.778 17.958     27.935 126.479  -4.534 37.215 351.483 93.460 33.327 489.590  4.865    357.719 2079.268
# 3 -3.612  7.977  9.643  8.757  3.430 19.140 13.700     21.485  86.373 -14.271 10.174 248.706 33.908 11.001 344.800 24.525    210.750 4067.683
# 4 -9.975  9.982  9.993  9.988  3.915 21.672 24.042     24.433 270.813 -12.413 24.149 161.634 53.974 19.768 227.102 11.682    180.759 2801.849
# 5 -3.465  1.972  8.401  3.735  1.407 12.140  7.893     10.904  45.399  -6.858 17.669 168.749 44.536 15.866 235.014 23.021    171.268 1506.739
# 6 -1.963  1.935  1.938  1.937  0.759  4.203 10.173      4.738  10.333  -3.243 26.890  80.292 44.227 16.952 121.495 22.870    119.396  451.759

# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'amplitude 2 components' = amplitude_2component_SPN,
        'tau 2 components' = tau_2component_SPN,
        'ctrl single examples' = ctrls_eg,
        'test single examples' = tests_eg,
        'statistics' = stats_summary
      )
    )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.csv'), wd=xlsx_path)
}
