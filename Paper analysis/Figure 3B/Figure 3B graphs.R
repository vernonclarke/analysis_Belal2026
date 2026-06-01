# =========================================# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- TRUE

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

identifier <- 'Figure 3B'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
# load RDATA
datasets <- c('NDNF GABA PSCs')  
datasets2 <- gsub('(?<!\\bX)\\s+|\\s+(?!X\\b)|-', '_', datasets, perl = TRUE)

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
# $NDNF_GABA_PSCs
#         A1 τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area1
# 1    0.000    NA     NA     NA     NA     NA     NA         NA     0.000  -83.338 27.970 220.487  66.139 23.931 308.240 22.490    235.288  24802.818
# 2    0.000    NA     NA     NA     NA     NA     NA         NA     0.000 -132.380 16.370 294.826  50.107 16.809 408.856 20.454    262.131  46259.139
# 3 -199.125 7.096 30.000 13.400  5.052 43.382 16.204     39.047  9337.491 -350.525 41.871 384.282 104.169 37.235 535.528 12.749    394.165 176642.619
# 4 -171.350 2.722 26.226  6.880  2.449 36.520 15.170     26.559  5841.865 -241.536 39.880 302.065  93.029 33.763 422.730 17.000    325.936  99273.913
# 5 -536.365 1.412 30.000  4.529  1.492 41.595 14.628     25.960 18712.888 -696.715 31.644 425.388  88.833 30.688 590.484 12.905    399.624 365202.220
# 6  -23.617 8.850  8.866  8.858  3.472 19.221 20.328     21.670   568.675  -47.292 20.882 135.323  46.145 16.939 190.369 16.019    152.810   9000.183
# 7  -81.882 7.095  8.363  7.694  3.014 16.786 19.314     18.863  1718.308 -162.242 30.230 313.420  78.245 27.680 436.008 18.000    311.818  65269.794
# 8 -153.808 5.357  5.363  5.360  2.101 11.631 19.129     13.113  2241.005 -146.833 26.691 259.408  67.659 24.069 361.181 16.042    262.146  49440.443

setwd(svg_path)

##

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
data <- cbind(-fits_list$NDNF_GABA_PSCs[,1], -fits_list$NDNF_GABA_PSCs[,10])
colnames(data) <- c('Afast', 'Aslow')

# DBSCAN_analyse(data) # eps=400
DBSCAN_analyse(data=data, eps=400, filename=paste0('dbscan.svg'), save=plotsave)

fits_list$NDNF_GABA_PSCs <- fits_list$NDNF_GABA_PSCs[-5,]
peaks_list$NDNF_GABA_PSCs <- peaks_list$NDNF_GABA_PSCs[-5]
areas_list$NDNF_GABA_PSCs <- areas_list$NDNF_GABA_PSCs[-5]
data_list$NDNF_GABA_PSCs <- data_list$NDNF_GABA_PSCs[-5]
output_list$NDNF_GABA_PSCs <- output_list$NDNF_GABA_PSCs[-5]

############################################################ AMPLITUDE ############################################################## 
amplitude_NDNF_GABA_PSCs <- create_df(matrix(peaks_list$NDNF_GABA_PSCs, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'SPN'))

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

BoxPlot3(formula=amplitude ~ cell_type, data=amplitude_NDNF_GABA_PSCs[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=NULL)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(10, 400) else c(-400, 0)
BoxPlot3(formula=amplitude ~ cell_type, data=transform(amplitude_NDNF_GABA_PSCs[,3:4], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_semilog_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

########################################################### CHARGE TRANSFER ###########################################################

area_NDNF_GABA_PSCs <- create_df(matrix(areas_list$NDNF_GABA_PSCs, ncol=1), levels = list(condition = 'GABAzine sensitive', cell_type = 'SPN'), var_name = 'charge_transfer')
area_NDNF_GABA_PSCs$charge_transfer <- -area_NDNF_GABA_PSCs$charge_transfer

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 200)
y_tick_interval <- 100

BoxPlot3(formula=charge_transfer ~ cell_type, data=area_NDNF_GABA_PSCs[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=NULL)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 200) else c(0, 200)

BoxPlot3(formula=charge_transfer ~ cell_type, data=area_NDNF_GABA_PSCs[,3:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_semilog_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')


######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_NDNF_GABA_PSCs <- create2condition_df(fits_list$NDNF_GABA_PSCs, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'))


# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ condition + Error(s), df=amplitude_2component_NDNF_GABA_PSCs)

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
y_tick_interval <- 200
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ condition, data=amplitude_2component_NDNF_GABA_PSCs, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

# paired scatter plot Afast vs Aslow
wide_df <- reshape(
  amplitude_2component_NDNF_GABA_PSCs,
  idvar = c('s', 'cell_type'),
  timevar = 'condition',
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
xlim=c(0, 400)
ylim=c(0, 400)
x_tick_interval <- 100
y_tick_interval <- 100

scatter_plot(wide_df, xlim=xlim, ylim=ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter.svg', lwd=lwd, save=plotsave)

##################################################### CHARGE TRANSFER FAST vs SLOW ####################################################

area_2component_NDNF_GABA_PSCs <- create2condition_df(fits_list$NDNF_GABA_PSCs, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='area')
area_2component_NDNF_GABA_PSCs$area <- area_2component_NDNF_GABA_PSCs$area/1e3

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=area ~ condition + Error(s), df=area_2component_NDNF_GABA_PSCs))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 200)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot3(formula=area ~ condition*cell_type, data=area_2component_NDNF_GABA_PSCs, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='area_2component_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

# paired scatter plot Areafast vs Areaslow
wide_df <- reshape(
  area_2component_NDNF_GABA_PSCs,
  idvar = c('s', 'cell_type'),
  timevar = 'condition',
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
xlim=c(0, 30)
ylim=c(0, 200)
x_tick_interval <- 10
y_tick_interval <- 50

scatter_plot(wide_df, xlim = xlim, ylim = ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval,
  height = height, width = width, xlab = expression(charge~transfer[fast]*' '(pC)), ylab = expression(charge~transfer[slow]*' '(pC)),
  filename = 'scatter_charge_transfer.svg', lwd = lwd, save = plotsave)

################################################### DECAY TIME CONSTANT FAST vs SLOW ##################################################

tau_2component_NDNF_GABA_PSCs <- create2condition_df(fits_list$NDNF_GABA_PSCs, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'SPN'), var_name='tau')

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ condition, df=tau_2component_NDNF_GABA_PSCs))

# update graph properties
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- c(0, 500)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot3(formula=tau ~ condition*cell_type, data=tau_2component_NDNF_GABA_PSCs, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[3,])
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')


BoxPlot3(formula=tau ~ condition*cell_type, data=tau_2component_NDNF_GABA_PSCs, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,500),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[3,], log_y=TRUE)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

ylab <- expression(tau[decay * ',' * fast] * ' ' * (ms))
yrange <- c(0, 30)
y_tick_interval <- 5

BoxPlot3(formula=tau ~ cell_type, data=tau_2component_NDNF_GABA_PSCs[tau_2component_NDNF_GABA_PSCs$condition == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=NULL)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')

ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
yrange <- c(0, 500)
y_tick_interval <- 100

BoxPlot3(formula=tau ~ cell_type, data=tau_2component_NDNF_GABA_PSCs[tau_2component_NDNF_GABA_PSCs$condition == 'slow',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=NULL)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_NDNF_GABA_PSCs.svg', width=width, height=height, bg='transparent')


########################################################## SINGLE EXAMPLES ###########################################################
# single examples

name <- datasets[1]
for (ii in 1:length(output_list$NDNF_GABA_PSCs)){
  traces <- output_list$NDNF_GABA_PSCs[[ii]]$traces
  func <- if (dim(traces)[2]==4) product1N else product2N
    fit_plot(traces=output_list$NDNF_GABA_PSCs[[ii]]$traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,1200), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}

# single egs
dx <- 0.1
xlab='time (ms)'
ylab='PSC amplitude (pA)'
lwd=1.0
filter='off'
xlim <- c(0, 1200)
ylim <- c(-250, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 4
SPNs_eg <- output_list$NDNF_GABA_PSCs[[ii]]$traces
single_fit_egs(traces=SPNs_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_NDNF_GABA_PSCs.svg', save=plotsave)
single_fit_egs(traces=SPNs_eg, xlim=c(100, 800), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_NDNF_GABA_PSCs.svg', log_y=TRUE, save=plotsave)

########################################################## RENAME OUTPUTS ############################################################

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
# $NDNF_GABA_PSCs
#      Afast τrise τdecay  tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay   tpeak r20_80  d80_20  delay half_width       area
# 1    0.000    NA     NA     NA     NA     NA     NA         NA    0.000  -83.338 27.970 220.487  66.139 23.931 308.240 22.490    235.288  24802.818
# 2    0.000    NA     NA     NA     NA     NA     NA         NA    0.000 -132.380 16.370 294.826  50.107 16.809 408.856 20.454    262.131  46259.139
# 3 -199.125 7.096 30.000 13.400  5.052 43.382 16.204     39.047 9337.491 -350.525 41.871 384.282 104.169 37.235 535.528 12.749    394.165 176642.619
# 4 -171.350 2.722 26.226  6.880  2.449 36.520 15.170     26.559 5841.865 -241.536 39.880 302.065  93.029 33.763 422.730 17.000    325.936  99273.913
# 6  -23.617 8.849  8.868  8.858  3.472 19.221 20.328     21.670  568.676  -47.292 20.882 135.323  46.145 16.940 190.369 16.019    152.810   9000.177
# 7  -81.882 7.095  8.363  7.694  3.014 16.786 19.314     18.863 1718.308 -162.242 30.230 313.420  78.245 27.680 436.008 18.000    311.818  65269.794
# 8 -153.808 5.357  5.363  5.360  2.101 11.631 19.129     13.113 2241.004 -146.833 26.691 259.408  67.659 24.069 361.181 16.042    262.146  49440.441

# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_NDNF_GABA_PSCs,
        'charge transfer'  = area_NDNF_GABA_PSCs,
        'amplitude 2 components' = amplitude_2component_NDNF_GABA_PSCs,
        'tau 2 components' = tau_2component_NDNF_GABA_PSCs,
        'SPN single examples' = SPNs_eg
      )
    )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  graphics.off()
}


