# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- TRUE

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

# settings
identifier <- 'Figure 9'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# path where all graphs are stored
svg_path <- paste0(analysis_path, '/svg')

# load RDATA

# names of datasets to be loaded
datasets <- c('ChAT-Cre X tdTomato',        # D1 dSPN
              'ChAT-Cre X tdTomato 6OHDA',  # D1 dSPN
              'ChAT-Cre X De eGFP',         # D2 iSPN
              'ChAT-Cre X De eGFP 6OHDA')   # D2 iSPN

datasets2 <- gsub('(?<!\\bX)\\s+|\\s+(?!X\\b)|-', '_', datasets, perl = TRUE)
datasets2 <- sub('_all_traces$', '', datasets2)


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
#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -56.038  3.829  3.832  3.830  1.501  8.311 18.404      9.370   583.438 -210.479 24.793 226.340  61.573 22.019 315.451 11.586    232.336  62533.64
# 2  -255.239  3.819 16.579  7.285  2.743 23.916 15.642     21.374  6566.786 -433.480 26.791 136.782  54.317 20.255 194.827 17.900    166.685  88198.38
# 3   -35.507  2.041  2.046  2.043  0.801  4.434 22.969      4.998   197.202 -145.341 24.325 152.419  53.117 19.545 214.718 11.392    173.890  31388.74
# 4  -291.374 19.742 23.672 21.588  8.457 47.155  0.000     52.952 17169.174 -532.137 49.593 256.724 101.061 37.653 365.337  0.000    311.411 202512.35
# 5  -336.909  1.399 12.237  3.426  1.229 17.068 13.833     12.700  5454.644 -508.551 20.904 141.027  46.850 17.149 198.091 10.827    157.342  99980.05
# 6  -180.028  2.521  6.184  3.819  1.475  9.724 12.428      9.976  2064.510 -238.322 14.853 129.534  36.334 13.042 180.683  9.087    134.542  40866.58
# 7   -80.479  2.686  8.855  4.600  1.756 13.192 17.954     12.659  1198.019  -83.755 17.660 138.627  41.700 15.093 193.821 14.254    148.104  15685.41
# 8   -11.571  8.829  9.152  8.989  3.523 19.510  7.458     21.992   282.773  -56.619 12.564 339.811  43.021 13.796 471.092 17.429    283.945  21836.57
# 9    -8.834  4.909  5.082  4.995  1.958 10.840  7.265     12.220   119.949  -44.599  4.020 208.204  16.180  4.789 288.632 17.469    162.036  10035.98
# 10  -27.649  4.184 10.879  6.497  2.503 16.896  6.483     17.124   546.552   -8.092 15.383  95.482  33.478 12.326 134.563  0.499    109.250   1097.13

# $ChAT_Cre_X_tdTomato_6OHDA
#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1 -12.480  6.041  6.078  6.059  2.375 13.149 10.384     14.824 205.556  -53.077  26.912 276.899  69.487 24.598 385.238 19.517    275.984 18888.979
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000 -104.272 145.881 145.949 145.915 57.198 316.623  1.332    374.162 41358.180
# 3 -13.210  6.387  6.394  6.391  2.505 13.867  5.154     15.634 229.482  -27.725  18.756 133.055  42.778 15.602 186.571 12.826    146.276  5087.855
# 4  -5.708  1.087 10.934  2.786  0.988 15.216  3.815     10.956  80.526   -2.411   4.501 423.278  20.672  5.648 586.788 13.059    315.665  1071.680
# 5 -13.978  3.410  3.414  3.412  1.338  7.404  5.992      8.348 129.651  -29.445 144.230 144.512 144.371 56.593 313.274  0.000    369.327 11555.395
# 6 -22.285  6.343  6.344  6.344  2.487 13.765  4.363     15.519 384.272  -28.608  19.438 266.530  54.897 18.928 369.934 24.772    249.379  9368.700
# 7 -11.396  9.582  9.815  9.698  3.801 21.046  5.664     23.725 300.440   -5.356  62.251 113.101  82.675 32.193 192.663  3.600    208.138  1258.278
# 8  -6.091  3.959  4.612  4.269  1.673  9.307  3.736     10.463  70.883   -3.474  28.442  36.264  32.037 12.545  70.347 10.500     78.745   304.736
# 9  -0.972 20.145 20.207 20.176  7.909 43.780  4.530     49.358  53.287   -1.444 114.863 115.927 115.393 45.233 250.399 18.086    284.434   453.001

# $ChAT_Cre_X_De_eGFP
#         A1 τrise τdecay tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1 -163.623 4.117 25.722 8.981  3.305 36.240 14.935     29.372 5967.360 -502.946 30.551 223.181 70.389 25.616 312.661 17.000    243.298 153869.15
# 2  -84.247 6.614  7.668 7.115  2.788 15.507 10.754     17.436 1633.832 -365.648 19.393 255.620 54.116 18.730 354.866 13.300    241.114 115504.78
# 3  -11.687 1.341  5.154 2.440  0.925  7.527 15.678      6.947   96.709  -69.100 34.474 245.407 78.721 28.703 344.074  4.687    269.522  23371.08
# 4  -72.115 4.205  4.208 4.206  1.649  9.127 20.066     10.290  824.548 -160.548 19.992 306.514 58.385 19.911 425.235 15.842    280.580  59536.03
# 5  -73.373 3.637  3.642 3.640  1.427  7.898 12.663      8.904  725.950 -241.665 22.674 104.786 44.290 16.614 150.345 15.409    132.131  38643.95
# 6  -55.194 3.369  3.372 3.370  1.321  7.313 14.062      8.245  505.645 -143.395 10.360 165.569 30.628 10.402 229.669 16.154    150.386  28566.18

# $ChAT_Cre_X_De_eGFP_6OHDA
#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width   area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# 1  -6.366  1.327 20.132  3.863  1.319 27.930  4.541     18.465 155.260   0.000      NA      NA      NA     NA      NA     NA         NA    0.000
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000  -8.239   1.974 239.117   9.546  2.519 331.487 24.628    175.974 2050.384
# 3 -12.837  2.219  2.223  2.221  0.871  4.820 18.845      5.434  77.506 -18.058  24.719 108.719  47.388 17.826 156.642 15.309    139.508 3035.804
# 4 -19.727  4.849  4.850  4.849  1.901 10.523  4.213     11.864 260.039 -21.597  20.611 219.211  53.786 18.985 304.861 12.858    216.784 6050.871
# 5  -7.150 22.028 22.073 22.050  8.644 47.848 21.671     53.944 428.592 -14.631 125.706 125.874 125.790 49.309 272.955  4.712    313.522 5002.983

setwd(svg_path)

############################################################## AMPLITUDE ############################################################## 

amplitude_dSPN_control <- create_df(matrix(peaks_list$ChAT_Cre_X_tdTomato, ncol=1), levels = list(condition = 'control', cell_type = 'dSPN'))
amplitude_iSPN_control <- create_df(matrix(peaks_list$ChAT_Cre_X_De_eGFP,  ncol=1), levels = list(condition = 'control', cell_type = 'iSPN'), 
  start_id = length(unique(amplitude_dSPN_control$s)) + 1)
amplitude_dSPN_6OHDA <- create_df(matrix(peaks_list$ChAT_Cre_X_tdTomato_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'dSPN'),
  start_id = length(unique(amplitude_dSPN_control$s)) + length(unique(amplitude_iSPN_control$s)) + 1)
amplitude_iSPN_6OHDA <- create_df(matrix(peaks_list$ChAT_Cre_X_De_eGFP_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'iSPN'),
  start_id = length(unique(amplitude_dSPN_control$s)) + length(unique(amplitude_iSPN_control$s)) + length(unique(amplitude_dSPN_6OHDA$s)) + 1)

amplitude_SPN <- rbind(amplitude_dSPN_control, amplitude_iSPN_control, amplitude_dSPN_6OHDA, amplitude_iSPN_6OHDA)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ cell_type*condition, df=amplitude_SPN)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-700, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ cell_type*condition, data=amplitude_SPN[,2:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ cell_type*condition, data=transform(amplitude_SPN[,2:4], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,700),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_amplitude_SPN.svg', width=width, height=height, bg='transparent')

# nb sig bars are offscale so needs to be added manually

########################################################### CHARGE TRANSFER ###########################################################
area_dSPN_control <- create_df(matrix(areas_list$ChAT_Cre_X_tdTomato, ncol=1), levels = list(condition = 'control', cell_type = 'dSPN'), var_name = 'charge_transfer')
area_iSPN_control <- create_df(matrix(areas_list$ChAT_Cre_X_De_eGFP,  ncol=1), levels = list(condition = 'control', cell_type = 'iSPN'), var_name = 'charge_transfer', 
  start_id = length(unique(area_dSPN_control$s)) + 1)
area_dSPN_6OHDA <- create_df(matrix(areas_list$ChAT_Cre_X_tdTomato_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'dSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN_control$s)) + length(unique(area_iSPN_control$s)) + 1)
area_iSPN_6OHDA <- create_df(matrix(areas_list$ChAT_Cre_X_De_eGFP_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'iSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN_control$s)) + length(unique(area_iSPN_control$s)) + length(unique(area_dSPN_6OHDA$s)) + 1)

area_SPN <- rbind(area_dSPN_control, area_iSPN_control, area_dSPN_6OHDA, area_iSPN_6OHDA)

area_SPN$charge_transfer <- -area_SPN$charge_transfer

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- FALSE
yrange <- if (log_y) c(0.1, 250) else c(0, 250)
y_tick_interval <- 50

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_SPN[,2:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[5:8,], log_y=FALSE,)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_SPN[,2:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,300),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'))
amplitude_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP,  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1)
amplitude_2component_dSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato_6OHDA, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'),
    start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1]  + dim(fits_list$ChAT_Cre_X_tdTomato_6OHDA)[1] + 1)
amplitude_2component_iSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP_6OHDA,  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + dim(fits_list$ChAT_Cre_X_tdTomato_6OHDA)[1] + dim(fits_list$ChAT_Cre_X_De_eGFP)[1]+ 1)

add_condition_column <- function(df, new_condition, var_name) {
  names(df)[names(df) == 'condition'] <- 'kinetics'
  df$condition <- new_condition
  df <- df[, c('s', 'cell_type', 'condition', 'kinetics', var_name)]
  return(df)
}

amplitude_2component_SPN <- rbind(
  add_condition_column(amplitude_2component_dSPN, 'ctrl', 'amplitude'),
  add_condition_column(amplitude_2component_iSPN, 'ctrl', 'amplitude'),
  add_condition_column(amplitude_2component_dSPN_6OHDA, '6OHDA', 'amplitude'),
  add_condition_column(amplitude_2component_iSPN_6OHDA, '6OHDA', 'amplitude')
  )

# specify order for plots
amplitude_2component_SPN$condition <- factor(amplitude_2component_SPN$condition, levels = c('ctrl', '6OHDA'))
amplitude_2component_SPN$cell_type <- factor(amplitude_2component_SPN$cell_type, levels = c('dSPN', 'iSPN'))
amplitude_2component_SPN$kinetics  <- factor(amplitude_2component_SPN$kinetics, levels = c('fast', 'slow'))

# MCwilcox(formula=amplitude ~ condition*cell_type, df=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='fast',])
# MCwilcox(formula=amplitude ~ kinetics*cell_type + Error(s), df=amplitude_2component_SPN[amplitude_2component_SPN$condition=='ctrl',-3])

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ cell_type*condition, df=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='fast',])
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=amplitude ~ cell_type*condition, df=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='slow',]))

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-600, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 100
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ cell_type*condition, data=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='fast',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,])
if (plotsave) save_graph(svg_path=svg_path, filename='fast_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

# boxplot not feasible as responses are fits and some amplitudes are designated zero by the fitting procedure
# no way to deal with this and have semilog boxplot accurate
# BoxPlot3(formula=amplitude ~ cell_type*condition, data=transform(amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='fast',c(3,2,5)], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=c(0.1,600),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[1:4,])
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_fast_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ cell_type*condition, data=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='slow_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

# BoxPlot3(formula=amplitude ~ cell_type*condition, data=transform(amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='slow',c(3,2,5)], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=c(0.1,600),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[5:8,])
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_slow_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

# scatter plots
df <- amplitude_2component_SPN
df1 <- df[df$cell_type=='dSPN',]
df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

df1$kinetics <- NULL
df1$x <- ifelse(df1$condition == 'ctrl', 1, 2)
df1$condition <- NULL
names(df1)[names(df1) == 'amplitude'] <- 'y'
df1 <- df1[, c('s', 'x', 'y')]
df1$s <- as.integer(as.character(df1$s))

width=2.875
height=3.5
xlim=c(0, 500)
ylim=c(0, 500)
ScatterPlot(df1, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter_dSPN.svg', lwd=lwd, save=plotsave)

# ScatterPlot(df1, sign=-1, xlim=c(1, 500), ylim=c(1, 500), height=height, width=width, filename='scatter_logxy_ctrl.svg', lwd=lwd, log_xy=TRUE, save=plotsave)

df2 <- df[df$cell_type=='iSPN',]
df2 <- df2[order(df2$kinetics != 'fast', df2$s), ]

df2$kinetics <- NULL
df2$x <- ifelse(df2$condition == 'ctrl', 1, 2)
df2$condition <- NULL
names(df2)[names(df2) == 'amplitude'] <- 'y'
df2 <- df2[, c('s', 'x', 'y')]
df2$s <- as.integer(as.character(df2$s))

ScatterPlot(df2, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, open_symbols=TRUE, filename='scatter_iSPN.svg', lwd=lwd, save=plotsave)

# ScatterPlot(df2, sign=-1, xlim=c(1, 500), ylim=c(1, 500), height=height, width=width, open_symbols=TRUE, filename='scatter_logxy_6OHDA.svg', lwd=lwd, log_xy=TRUE, save=plotsave)

#################################################### CHARGE TRANSFER FAST vs SLOW ####################################################
area_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='charge_transfer')
area_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1)
area_2component_dSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato_6OHDA, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1]  + dim(fits_list$ChAT_Cre_X_tdTomato_6OHDA)[1] + 1)
area_2component_iSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP_6OHDA, cols = c(9, 18),  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + dim(fits_list$ChAT_Cre_X_tdTomato_6OHDA)[1] + dim(fits_list$ChAT_Cre_X_De_eGFP)[1]+ 1)

area_2component_SPN <- rbind(
  add_condition_column(area_2component_dSPN, 'ctrl', 'charge_transfer'),
  add_condition_column(area_2component_iSPN, 'ctrl', 'charge_transfer'),
  add_condition_column(area_2component_dSPN_6OHDA, '6OHDA', 'charge_transfer'),
  add_condition_column(area_2component_iSPN_6OHDA, '6OHDA', 'charge_transfer')
  )
area_2component_SPN$charge_transfer <- area_2component_SPN$charge_transfer/1e3

area_2component_SPN$condition <- factor(area_2component_SPN$condition, levels = c('ctrl', '6OHDA'))
area_2component_SPN$cell_type <- factor(area_2component_SPN$cell_type, levels = c('dSPN', 'iSPN'))
area_2component_SPN$kinetics  <- factor(area_2component_SPN$kinetics, levels = c('fast', 'slow'))

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_2component_SPN[area_2component_SPN$kinetics=='fast',]))
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_2component_SPN[area_2component_SPN$kinetics=='slow',]))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 10)
y_tick_interval <- 2
width <- 3
height <- 3.5

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='fast',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[9:12,])
if (plotsave) save_graph(svg_path=svg_path, filename='fast_area_2component_SPN.svg', width=width, height=height, bg='transparent')

# BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='fast',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=c(0.01,10),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[9:12,])
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_fast_area_2component_SPN.svg', width=width, height=height, bg='transparent')

yrange <- c(0, 250)
y_tick_interval <- 50

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[13:16,])
if (plotsave) save_graph(svg_path=svg_path, filename='slow_area_2component_SPN.svg', width=width, height=height, bg='transparent')

# BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=c(0.01,250),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[13:16,])
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_slow_area_2component_SPN.svg', width=width, height=height, bg='transparent')

# # scatter plots
# df <- area_2component_SPN
# df1 <- df[df$condition=='ctrl',]
# df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

# df1$kinetics <- NULL
# df1$x <- ifelse(df1$cell_type == 'dSPN', 1, 2)
# df1$condition <- NULL
# names(df1)[names(df1) == 'charge_transfer'] <- 'y'
# df1 <- df1[, c('s', 'x', 'y')]
# df1$s <- as.integer(as.character(df1$s))

# width=2.875
# height=3.5
# xlim=c(0, 50)
# ylim=c(0, 50)
# x_tick_interval=10 
# y_tick_interval=10

# ScatterPlot(df1, xlim=xlim, ylim=ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, filename='scatter_charge_transfer_ctrl.svg', lwd=lwd, save=plotsave)

# df2 <- df[df$condition=='6OHDA',]
# df2 <- df2[order(df2$kinetics != 'fast', df2$s), ]

# df2$kinetics <- NULL
# df2$x <- ifelse(df2$cell_type == 'dSPN', 1, 2)
# df2$condition <- NULL
# names(df2)[names(df2) == 'charge_transfer'] <- 'y'
# df2 <- df2[, c('s', 'x', 'y')]
# df2$s <- as.integer(as.character(df2$s))
# ScatterPlot(df2, xlim=xlim, ylim=ylim, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, height=height, width=width, open_symbols=TRUE, filename='scatter_charge_transfer_6OHDA.svg', lwd=lwd, save=plotsave)

################################################### DECAY TIME CONSTANT FAST vs SLOW ##################################################
tau_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau')
tau_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='tau',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1)
tau_2component_dSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato_6OHDA, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1]  + dim(fits_list$ChAT_Cre_X_De_eGFP)[1] + 1)
tau_2component_iSPN_6OHDA <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP_6OHDA, cols = c(3, 12),  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='tau',
  start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + dim(fits_list$ChAT_Cre_X_De_eGFP)[1] + dim(fits_list$ChAT_Cre_X_tdTomato_6OHDA)[1]+ 1)

tau_2component_SPN <- rbind(
  add_condition_column(tau_2component_dSPN, 'ctrl', 'tau'),
  add_condition_column(tau_2component_iSPN, 'ctrl', 'tau'),
  add_condition_column(tau_2component_dSPN_6OHDA, '6OHDA', 'tau'),
  add_condition_column(tau_2component_iSPN_6OHDA, '6OHDA', 'tau')
  )

# specify order for plots
tau_2component_SPN$condition <- factor(tau_2component_SPN$condition, levels = c('ctrl', '6OHDA'))
tau_2component_SPN$cell_type <- factor(tau_2component_SPN$cell_type, levels = c('dSPN', 'iSPN'))
tau_2component_SPN$kinetics  <- factor(tau_2component_SPN$kinetics, levels = c('fast', 'slow'))

# tau_2component_dSPN <- create2condition_df(fits_list$ChAT_Cre_X_tdTomato, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau')
# tau_2component_iSPN <- create2condition_df(fits_list$ChAT_Cre_X_De_eGFP,  cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
#   start_id = dim(fits_list$ChAT_Cre_X_tdTomato)[1] + 1, var_name='tau')
# tau_2component_SPN <- rbind(tau_2component_dSPN, tau_2component_iSPN)

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ kinetics*condition , tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')]))
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ kinetics*condition, tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')]))

# update graph properties
log_y <- FALSE
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- if (log_y) c(1,600) else c(0,600)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_dSPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1,600) else c(0,600)

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_dSPN.svg', width=width, height=height, bg='transparent')

log_y <- FALSE
yrange <- if (log_y) c(1, 600) else c(0,600)
BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[21:24,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_iSPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 600) else c(0,600)

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[21:24,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_iSPN.svg', width=width, height=height, bg='transparent')


# # ctrl vs 6OHDA plots
# log_y <- FALSE
# yrange <- if (log_y) c(1,600) else c(0,600)

# BoxPlot3(formula=tau ~ kinetics*cell_type, data=tau_2component_SPN[tau_2component_SPN$condition=='ctrl', setdiff(names(tau_2component_SPN), 'condition')], 
#   wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
# if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_ctrl.svg', width=width, height=height, bg='transparent')

# log_y <- TRUE
# yrange <- if (log_y) c(1,600) else c(0,600)

# BoxPlot3(formula=tau ~ kinetics*cell_type, data=tau_2component_SPN[tau_2component_SPN$condition=='ctrl', setdiff(names(tau_2component_SPN), 'condition')], 
#   wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_ctrl.svg', width=width, height=height, bg='transparent')


# log_y <- FALSE
# yrange <- if (log_y) c(1,600) else c(0,600)

# BoxPlot3(formula=tau ~ kinetics*cell_type, data=tau_2component_SPN[tau_2component_SPN$condition=='6OHDA', setdiff(names(tau_2component_SPN), 'condition')], 
#   wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
# if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_6OHDA.svg', width=width, height=height, bg='transparent')

# log_y <- TRUE
# yrange <- if (log_y) c(1,600) else c(0,600)

# BoxPlot3(formula=tau ~ kinetics*cell_type, data=tau_2component_SPN[tau_2component_SPN$condition=='6OHDA', setdiff(names(tau_2component_SPN), 'condition')], 
#   wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
#   yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#   p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
# if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_6OHDA.svg', width=width, height=height, bg='transparent')


# nb these are unpaired tests as there are missing values

# update graph properties

# fast or slow as separate graphs 
stats_summary <- rbind(stats_summary, 
  MCwilcox(formula=tau ~ cell_type*condition, df=tau_2component_SPN[tau_2component_SPN$kinetics=='fast', setdiff(names(tau_2component_SPN), 'kinetics')])
  )

stats_summary <- rbind(stats_summary, 
  MCwilcox(formula=tau ~ cell_type*condition, df=tau_2component_SPN[tau_2component_SPN$kinetics=='slow', setdiff(names(tau_2component_SPN), 'kinetics')])
)

log_y <- FALSE
ylab <- expression(tau[decay * ',' * fast] * ' ' * (ms))
yrange <- if (log_y) c(1, 30) else c(0,30)
y_tick_interval <- 5

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='fast', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 30) else c(0,30)

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='fast', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, 
  lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_fast_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- FALSE
ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
yrange <- if (log_y) c(10, 600) else c(0,600)
y_tick_interval <- 100

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='slow', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[13:16,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(10, 600) else c(0,600)

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='slow', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, 
  lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_slow_2component_SPN.svg', width=width, height=height, bg='transparent')

# ########################################################## SINGLE EXAMPLES ###########################################################
# # single examples

median(c(peaks_list[[1]], peaks_list[[3]]))
# [1] -190.8226

median(c(peaks_list[[2]], peaks_list[[4]]))
# [1] -21.20044

for (iii in 1:length(output_list)) {
  name <- names(output_list)[iii]
  for (ii in 1:length(output_list[[iii]])){
    traces <- output_list[[iii]][[ii]]$traces
    func <- if (dim(traces)[2]==4) product1N else product2N  
    fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
      ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
  }
}
graphics.off()

# single egs
dx <- 0.1
xlab <- 'time (ms)'
ylab <- 'PSC amplitude (pA)'
lwd <- 1.0
filter='off'
xlim <- c(0, 800)
ylim <- c(-250, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii <- 7
dSPNs_control_egs <- output_list$ChAT_Cre_X_tdTomato[[ii]]$traces
single_fit_egs(traces=dSPNs_control_egs, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_dSPN.svg', save=plotsave)
single_fit_egs(traces=dSPNs_control_egs, xlim=c(100, 500), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_dSPN.svg', log_y=TRUE, save=plotsave)

ii = 7
dSPNs_6OHDA_egs <- output_list$ChAT_Cre_X_tdTomato_6OHDA[[ii]]$traces
single_fit_egs(traces=dSPNs_6OHDA_egs, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_dSPN_6OHDA.svg', save=plotsave)
single_fit_egs(traces=dSPNs_6OHDA_egs, xlim=c(100, 500), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_dSPN_6OHDA.svg', log_y=TRUE, save=plotsave)

graphics.off()

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
#       Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width      area    Aslow  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area
# 1   -56.038  3.829  3.832  3.830  1.501  8.311 18.404      9.370   583.438 -210.479 24.793 226.340  61.573 22.019 315.451 11.586    232.336  62533.64
# 2  -255.239  3.819 16.579  7.285  2.743 23.916 15.642     21.374  6566.786 -433.480 26.791 136.782  54.317 20.255 194.827 17.900    166.685  88198.38
# 3   -35.507  2.041  2.046  2.043  0.801  4.434 22.969      4.998   197.202 -145.341 24.325 152.419  53.117 19.545 214.718 11.392    173.890  31388.74
# 4  -291.374 19.742 23.672 21.588  8.457 47.155  0.000     52.952 17169.174 -532.137 49.593 256.724 101.061 37.653 365.337  0.000    311.411 202512.35
# 5  -336.909  1.399 12.237  3.426  1.229 17.068 13.833     12.700  5454.644 -508.551 20.904 141.027  46.850 17.149 198.091 10.827    157.342  99980.05
# 6  -180.028  2.521  6.184  3.819  1.475  9.724 12.428      9.976  2064.510 -238.322 14.853 129.534  36.334 13.042 180.683  9.087    134.542  40866.58
# 7   -80.479  2.686  8.855  4.600  1.756 13.192 17.954     12.659  1198.019  -83.755 17.660 138.627  41.700 15.093 193.821 14.254    148.104  15685.41
# 8   -11.571  8.829  9.152  8.989  3.523 19.510  7.458     21.992   282.773  -56.619 12.564 339.811  43.021 13.796 471.092 17.429    283.945  21836.57
# 9    -8.834  4.909  5.082  4.995  1.958 10.840  7.265     12.220   119.949  -44.599  4.020 208.204  16.180  4.789 288.632 17.469    162.036  10035.98
# 10  -27.649  4.184 10.879  6.497  2.503 16.896  6.483     17.124   546.552   -8.092 15.383  95.482  33.478 12.326 134.563  0.499    109.250   1097.13

# $ChAT_Cre_X_tdTomato_6OHDA
#     Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area    Aslow   τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area
# 1 -12.480  6.041  6.078  6.059  2.375 13.149 10.384     14.824 205.556  -53.077  26.912 276.899  69.487 24.598 385.238 19.517    275.984 18888.979
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000 -104.272 145.881 145.949 145.915 57.198 316.623  1.332    374.162 41358.180
# 3 -13.210  6.387  6.394  6.391  2.505 13.867  5.154     15.634 229.482  -27.725  18.756 133.055  42.778 15.602 186.571 12.826    146.276  5087.855
# 4  -5.708  1.087 10.934  2.786  0.988 15.216  3.815     10.956  80.526   -2.411   4.501 423.278  20.672  5.648 586.788 13.059    315.665  1071.680
# 5 -13.978  3.410  3.414  3.412  1.338  7.404  5.992      8.348 129.651  -29.445 144.230 144.512 144.371 56.593 313.274  0.000    369.327 11555.395
# 6 -22.285  6.343  6.344  6.344  2.487 13.765  4.363     15.519 384.272  -28.608  19.438 266.530  54.897 18.928 369.934 24.772    249.379  9368.700
# 7 -11.396  9.582  9.815  9.698  3.801 21.046  5.664     23.725 300.440   -5.356  62.251 113.101  82.675 32.193 192.663  3.600    208.138  1258.278
# 8  -6.091  3.959  4.612  4.269  1.673  9.307  3.736     10.463  70.883   -3.474  28.442  36.264  32.037 12.545  70.347 10.500     78.745   304.736
# 9  -0.972 20.145 20.207 20.176  7.909 43.780  4.530     49.358  53.287   -1.444 114.863 115.927 115.393 45.233 250.399 18.086    284.434   453.001

# $ChAT_Cre_X_De_eGFP
#      Afast τrise τdecay tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area
# 1 -163.623 4.117 25.722 8.981  3.305 36.240 14.935     29.372 5967.360 -502.946 30.551 223.181 70.389 25.616 312.661 17.000    243.298 153869.15
# 2  -84.247 6.614  7.668 7.115  2.788 15.507 10.754     17.436 1633.832 -365.648 19.393 255.620 54.116 18.730 354.866 13.300    241.114 115504.78
# 3  -11.687 1.341  5.154 2.440  0.925  7.527 15.678      6.947   96.709  -69.100 34.474 245.407 78.721 28.703 344.074  4.687    269.522  23371.08
# 4  -72.115 4.205  4.208 4.206  1.649  9.127 20.066     10.290  824.548 -160.548 19.992 306.514 58.385 19.911 425.235 15.842    280.580  59536.03
# 5  -73.373 3.637  3.642 3.640  1.427  7.898 12.663      8.904  725.950 -241.665 22.674 104.786 44.290 16.614 150.345 15.409    132.131  38643.95
# 6  -55.194 3.369  3.372 3.370  1.321  7.313 14.062      8.245  505.645 -143.395 10.360 165.569 30.628 10.402 229.669 16.154    150.386  28566.18

# $ChAT_Cre_X_De_eGFP_6OHDA
#     Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area   Aslow   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area
# 1  -6.366  1.327 20.132  3.863  1.319 27.930  4.541     18.465 155.260   0.000      NA      NA      NA     NA      NA     NA         NA    0.000
# 2   0.000     NA     NA     NA     NA     NA     NA         NA   0.000  -8.239   1.974 239.117   9.546  2.519 331.487 24.628    175.974 2050.384
# 3 -12.837  2.219  2.223  2.221  0.871  4.820 18.845      5.434  77.506 -18.058  24.719 108.719  47.388 17.826 156.642 15.309    139.508 3035.804
# 4 -19.727  4.849  4.850  4.849  1.901 10.523  4.213     11.864 260.039 -21.597  20.611 219.211  53.786 18.985 304.861 12.858    216.784 6050.871
# 5  -7.150 22.028 22.073 22.050  8.644 47.848 21.671     53.944 428.592 -14.631 125.706 125.874 125.790 49.309 272.955  4.712    313.522 5002.983


# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'amplitude 2 components' = amplitude_2component_SPN,
        'tau 2 components' = tau_2component_SPN,
        'dSPN single examples control' = dSPNs_control_egs,
        'dSPN single examples 6OHDA' = dSPNs_6OHDA_egs,
        'statistics' = stats_summary,
        'additional statistics' = stats_summary1
      )
    )
  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.csv'), wd=xlsx_path)
}


