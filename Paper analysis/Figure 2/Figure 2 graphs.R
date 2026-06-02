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

identifier <- 'Figure 2'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}
# load RDATA
datasets <- c('ChAT-Cre X tdTomato',  # D1 dSPN
            'ChAT-Cre X De eGFP')     # D2 iSPN

# Remove _nwb suffix if present
datasets_clean <- sub('_nwb$', '', datasets)
datasets2 <- gsub('(?<!\\bX)\\s+|\\s+(?!X\\b)|-', '_', datasets_clean, perl = TRUE)


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
# $ChAT_Cre_X_tdTomato
#          A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -88.305 3.935  3.942 3.938  1.544  8.546 11.360      9.635  945.380 -117.851 11.549  97.056 27.905 10.046 135.483 13.294    101.797 15248.338
# 2   -88.606 6.407  7.100 6.742  2.642 14.660 25.720     16.507 1625.942  -70.643  7.041 259.004 26.094  8.065 359.056 31.497    208.470 20236.213
# 3   -79.712 1.442 31.399 4.655  1.530 43.534 14.273     27.066 2902.916 -210.653 21.857 157.532 50.125 18.260 220.782 27.329    172.395 45616.726
# 4   -69.261 2.956 13.383 5.730  2.152 19.235 10.288     16.999 1422.304  -41.612  3.047 186.091 12.738  3.689 257.977 11.467    142.867  8292.307
# 5  -132.419 2.720  2.723 2.721  1.067  5.905 14.385      6.658  979.584 -177.149 10.615  90.084 25.732  9.257 125.724 16.978     94.239 21234.360
# 6   -51.147 6.279  6.280 6.279  2.461 13.626 22.778     15.362  873.039  -45.042  3.075 117.152 11.494  3.537 162.408 18.579     93.934  5820.659
# 7  -101.555 1.565  7.671 3.125  1.168 10.958 10.127      9.478 1170.867  -86.024  8.587  75.015 21.018  7.543 104.633  7.961     77.880  8539.756
# 8   -23.844 7.794  8.711 8.236  3.228 17.915 20.381     20.168  534.621  -10.815  6.774 595.574 30.671  8.458 825.641 18.800    445.916  6781.617
# 9    -3.217 0.778  2.969 1.412  0.535  4.340 24.680      4.013   15.368  -10.282 18.944 230.012 51.541 17.977 319.489 10.651    220.773  2958.884
# 10 -127.557 3.980  3.991 3.986  1.562  8.648 15.773      9.750 1381.940 -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.689

# $ChAT_Cre_X_De_eGFP
#        A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -65.476 9.646 10.032 9.836  3.856 21.350  8.937     24.066 1751.023 -221.043 23.193 168.395 53.323 19.414 235.953 12.500    183.894 51088.809
# 2 -73.551 4.967  4.971 4.969  1.948 10.782 17.850     12.156  993.478 -205.394 19.079 172.792 47.259 16.911 240.854 20.759    177.732 46654.255
# 3 -26.094 6.899  6.910 6.905  2.707 14.983  5.336     16.892  489.755  -19.980  2.996 140.667 11.783  3.532 195.006 11.173    110.451  3056.172
# 4 -81.963 1.274  6.060 2.516  0.942  8.675 12.862      7.565  752.274 -162.985  4.713  71.321 13.711  4.682  98.950  6.609     65.450 14088.337
# 5 -80.942 4.700  4.712 4.706  1.845 10.212 17.195     11.513 1035.417 -140.761 18.779 141.793 43.759 15.885 198.452 22.078    153.133 27174.772
# 6 -88.302 5.124  5.128 5.126  2.009 11.123 18.266     12.540 1230.388 -177.971 27.222 232.987 66.176 23.789 325.109  5.281    243.187 55085.135
setwd(svg_path)

############################################################## AMPLITUDE ############################################################## 

amplitude_dSPN <- create_df(matrix(peaks_list$ChAT_Cre_X_tdTomato, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'dSPN'))
amplitude_iSPN <- create_df(matrix(peaks_list$ChAT_Cre_X_De_eGFP,  ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'iSPN'), 
  start_id = length(unique(amplitude_dSPN$s)) + 1)
amplitude_SPN <- rbind(amplitude_dSPN, amplitude_iSPN)

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

BoxPlot3(formula=amplitude ~ cell_type, data=amplitude_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')


BoxPlot3(formula=amplitude ~ cell_type, data=transform(amplitude_SPN[,3:4], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,300),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1:2,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_semilog_SPN.svg', width=width, height=height, bg='transparent')


########################################################### CHARGE TRANSFER ###########################################################

area_dSPN <- create_df(matrix(areas_list$ChAT_Cre_X_tdTomato, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'dSPN'), var_name = 'charge_transfer')
area_iSPN <- create_df(matrix(areas_list$ChAT_Cre_X_De_eGFP,  ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'iSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN$s)) + 1)
area_SPN <- rbind(area_dSPN, area_iSPN)
area_SPN$charge_transfer <- -area_SPN$charge_transfer

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 80)
y_tick_interval <- 20

BoxPlot3(formula=charge_transfer ~ cell_type, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')


BoxPlot3(formula=charge_transfer ~ cell_type, data=area_SPN[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,100),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1:2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_semilog_SPN.svg', width=width, height=height, bg='transparent')



######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'))
amplitude_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP,  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1)
amplitude_2component_SPN <- rbind(amplitude_2component_dSPN, amplitude_2component_iSPN)
names(amplitude_2component_SPN)[names(amplitude_2component_SPN) == 'condition'] <- 'kinetics'

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ kinetics*cell_type + Error(s), df=amplitude_2component_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-400, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ kinetics*cell_type, data=amplitude_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ kinetics*cell_type, data=transform(amplitude_2component_SPN, amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,400),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[1:4,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_semilog_SPN.svg', width=width, height=height, bg='transparent')


# paired scatter plot Afast vs Aslow
wide_df <- reshape(
  amplitude_2component_SPN,
  idvar = c('s', 'cell_type'),
  timevar = 'kinetics',
  direction = 'wide'
)
names(wide_df)[names(wide_df) == 'amplitude.fast'] <- 'Afast'
names(wide_df)[names(wide_df) == 'amplitude.slow'] <- 'Aslow'
wide_df$cell_type <- ifelse(wide_df$cell_type == 'dSPN', 1, 2)
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

area_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='area')
area_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP,  cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1, var_name='area')
area_2component_SPN <- rbind(area_2component_dSPN, area_2component_iSPN)
area_2component_SPN$area <- area_2component_SPN$area/1e3

names(area_2component_SPN)[names(area_2component_SPN) == 'condition'] <- 'kinetics'


# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=area ~ kinetics*cell_type + Error(s), df=area_2component_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 60)
y_tick_interval <- 20
width <- 3
height <- 3.5

BoxPlot3(formula=area ~ kinetics*cell_type, data=area_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='area_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=area ~ kinetics*cell_type, data=area_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.01,60),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_semilog_SPN.svg', width=width, height=height, bg='transparent')


# paired scatter plot Areafast vs Areaslow
wide_df <- reshape(
  area_2component_SPN,
  idvar = c('s', 'cell_type'),
  timevar = 'kinetics',
  direction = 'wide'
)
names(wide_df)[names(wide_df) == 'area.fast'] <- 'Areafast'
names(wide_df)[names(wide_df) == 'area.slow'] <- 'Areaslow'
wide_df$cell_type <- ifelse(wide_df$cell_type == 'dSPN', 1, 2)
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

tau_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau')
tau_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP,  cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1, var_name='tau')
tau_2component_SPN <- rbind(tau_2component_dSPN, tau_2component_iSPN)

# create output for stats
stats_summary <-  rbind(stats_summary,  MCwilcox(formula=tau ~ condition*cell_type + Error(s), df=tau_2component_SPN))
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ condition*cell_type + Error(s), df=tau_2component_SPN))

# update graph properties
log_y <- TRUE
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- if (log_y) c(1, 600) else c(0, 600)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot3(formula=tau ~ condition*cell_type, data=tau_2component_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[3:6,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_dSPN.svg', width=width, height=height, bg='transparent')

ylab <- expression(tau[decay * ',' * fast] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 30) else c(0, 30)
y_tick_interval <- 5

BoxPlot3(formula=tau ~ cell_type, data=tau_2component_SPN[tau_2component_SPN$condition == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[3,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_SPN.svg', width=width, height=height, bg='transparent')


ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 600) else c(0, 600)
y_tick_interval <- 100

BoxPlot3(formula=tau ~ cell_type, data=tau_2component_SPN[tau_2component_SPN$condition == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[4,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_SPN.svg', width=width, height=height, bg='transparent')


########################################################## SINGLE EXAMPLES ###########################################################
# single examples

name <- datasets[1]
for (ii in 1:length(output_list$ChAT_Cre_X_tdTomato)){
  fit_plot(traces=output_list$ChAT_Cre_X_tdTomato[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}

name <- datasets[2]
for (ii in 1:length(output_list$ChAT_Cre_X_De_eGFP)){
  fit_plot(traces=output_list$ChAT_Cre_X_De_eGFP[[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}
graphics.off()

# single egs
dx <- 0.1
xlab='time (ms)'
ylab='PSC amplitude (pA)'
lwd=1.0
filter='off'
xlim <- c(50, 400)
ylim <- c(-200, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 7
dSPNs_eg <- output_list$ChAT_Cre_X_tdTomato[[ii]]$traces
single_fit_egs(traces=dSPNs_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_dSPN.svg', save=plotsave)
single_fit_egs(traces=dSPNs_eg, xlim=c(108, 348), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_dSPN.svg', log_y=TRUE, save=plotsave)

ylim <- c(-250, 10)
ii = 4
iSPNs_eg <- output_list$ChAT_Cre_X_De_eGFP[[ii]]$traces
single_fit_egs(traces=iSPNs_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_iSPN.svg', save=plotsave)
single_fit_egs(traces=iSPNs_eg, xlim=c(100, 340), ylim=NULL, lwd=lwd, colors=colors, 
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_iSPN.svg', log_y=TRUE, save=plotsave)

############################################################ SAVE OUTPUTS #############################################################

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
# $ChAT_Cre_X_tdTomato
#       Afast τrise τdecay tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area
# 1   -88.305 3.935  3.942 3.938  1.544  8.546 11.360      9.635  945.380 -117.851 11.549  97.056 27.905 10.046 135.483 13.294    101.797 15248.338
# 2   -88.606 6.407  7.100 6.742  2.642 14.660 25.720     16.507 1625.942  -70.643  7.041 259.004 26.094  8.065 359.056 31.497    208.470 20236.213
# 3   -79.712 1.442 31.399 4.655  1.530 43.534 14.273     27.066 2902.916 -210.653 21.857 157.532 50.125 18.260 220.782 27.329    172.395 45616.726
# 4   -69.261 2.956 13.383 5.730  2.152 19.235 10.288     16.999 1422.304  -41.612  3.047 186.091 12.738  3.689 257.977 11.467    142.867  8292.307
# 5  -132.419 2.720  2.723 2.721  1.067  5.905 14.385      6.658  979.584 -177.149 10.615  90.084 25.732  9.257 125.724 16.978     94.239 21234.360
# 6   -51.147 6.279  6.280 6.279  2.461 13.626 22.778     15.362  873.039  -45.042  3.075 117.152 11.494  3.537 162.408 18.579     93.934  5820.659
# 7  -101.555 1.565  7.671 3.125  1.168 10.958 10.127      9.478 1170.867  -86.024  8.587  75.015 21.018  7.543 104.633  7.961     77.880  8539.756
# 8   -23.844 7.794  8.711 8.236  3.228 17.915 20.381     20.168  534.621  -10.815  6.774 595.574 30.671  8.458 825.641 18.800    445.916  6781.617
# 9    -3.217 0.778  2.969 1.412  0.535  4.340 24.680      4.013   15.368  -10.282 18.944 230.012 51.541 17.977 319.489 10.651    220.773  2958.884
# 10 -127.557 3.980  3.991 3.986  1.562  8.648 15.773      9.750 1381.940 -272.359 19.267 125.422 42.643 15.649 176.410 12.011    141.436 47992.689

# $ChAT_Cre_X_De_eGFP
#     Afast τrise τdecay tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area
# 1 -65.476 9.646 10.032 9.836  3.856 21.350  8.937     24.066 1751.023 -221.043 23.193 168.395 53.323 19.414 235.953 12.500    183.894 51088.809
# 2 -73.551 4.967  4.971 4.969  1.948 10.782 17.850     12.156  993.478 -205.394 19.079 172.792 47.259 16.911 240.854 20.759    177.732 46654.255
# 3 -26.094 6.899  6.910 6.905  2.707 14.983  5.336     16.892  489.755  -19.980  2.996 140.667 11.783  3.532 195.006 11.173    110.451  3056.172
# 4 -81.963 1.274  6.060 2.516  0.942  8.675 12.862      7.565  752.274 -162.985  4.713  71.321 13.711  4.682  98.950  6.609     65.450 14088.337
# 5 -80.942 4.700  4.712 4.706  1.845 10.212 17.195     11.513 1035.417 -140.761 18.779 141.793 43.759 15.885 198.452 22.078    153.133 27174.772
# 6 -88.302 5.124  5.128 5.126  2.009 11.123 18.266     12.540 1230.388 -177.971 27.222 232.987 66.176 23.789 325.109  5.281    243.187 55085.135

# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'amplitude 2 components' = amplitude_2component_SPN,
        'tau 2 components' = tau_2component_SPN,
        'dSPN single examples' = dSPNs_eg,
        'iSPN single examples' = iSPNs_eg,
        'statistics' = stats_summary
      )
    )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
}



