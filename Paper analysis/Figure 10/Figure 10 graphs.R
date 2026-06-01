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
identifier <- 'Figure 10'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# path where all graphs are stored
svg_path <- paste0(analysis_path, '/svg')

# names of datasets to be loaded
datasets <- c('NPY Cre X dSPN WT',          # D1 dSPN
              'NPY Cre X dSPN 6OHDA',       # D1 dSPN
              'NPY Cre X iSPN WT',          # D2 iSPN
              'NPY Cre X iSPN 6OHDA')       # D2 iSPN

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

# $NPY_Cre_X_dSPN_WT
#          A1  τrise τdecay  tpeak r20_80 d80_20 delay half_width    area1      A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1   -20.740  7.432  7.451  7.441  2.917 16.147 5.666     18.204  419.512 -20.298   3.010 110.691  11.153  3.447 153.451  8.069     89.097  2485.035
# 2   -40.701  3.945  8.127  5.540  2.151 13.344 4.696     14.138  654.040 -11.121   2.316  85.005   8.577  2.652 117.842  0.000     68.435  1045.684
# 3  -123.390  1.685 11.723  3.817  1.394 16.448 6.436     12.960 2003.107 -32.411   1.911 116.823   7.990  2.314 161.951  4.362     89.680  4054.373
# 4  -179.594  3.031 16.531  6.296  2.338 23.449 6.353     19.708 4344.955 -37.949  54.203  54.352  54.277 21.276 117.777 19.643    132.783  5599.033
# 5  -138.271  2.182  6.505  3.586  1.375  9.845 5.142      9.678 1560.914 -66.797   2.606 134.418  10.478  3.103 186.343  9.084    104.648  9706.684
# 6   -15.969  7.562  7.726  7.643  2.996 16.587 4.150     18.699  331.793  -2.981  20.166  75.856  36.392 13.803 111.060 20.000    103.080   365.297
# 7   -31.382  3.527 11.098  5.926  2.266 16.649 5.082     16.159  594.049 -20.130   3.323 120.873  12.281  3.801 167.566  6.558     97.411  2693.428
# 8   -16.752 11.340 11.773 11.554  4.529 25.078 4.400     28.268  526.204 -10.555  16.807 117.288  38.115 13.918 164.551 19.846    129.560  1713.322
# 9   -20.204 17.458 60.535 30.505 11.616 89.542 7.147     84.853 2024.401  -6.449 173.471 174.316 173.893 68.165 377.335 37.014    466.359  3048.475
# 10  -19.817  5.441  5.448  5.445  2.134 11.815 4.869     13.320  293.298 -86.816  17.612 317.013  53.898 18.082 439.625  6.030    281.885 32622.547
# 11   -6.155  6.563  8.606  7.492  2.933 16.499 6.040     18.436  126.511 -22.931  16.730 228.626  47.201 16.280 317.329  3.767    214.061  6444.854
# 12  -73.848  2.042 10.417  4.138  1.543 14.838 3.622     12.697 1144.483 -44.866   1.292  81.211   5.438  1.569 112.583  5.364     62.211  3895.942
# 13  -76.051  2.489 11.443  4.853  1.821 16.426 3.981     14.456 1329.930 -32.462   1.959  93.919   7.743  2.314 130.200  5.719     73.601  3310.801

# $NPY_Cre_X_dSPN_6OHDA
#        A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width    area1       A1  τrise  τdecay   tpeak r20_80  d80_20 delay half_width     area1
# 1 -36.672  6.523 62.230 16.436  5.856 86.669  6.720     63.177 2971.913  -15.339 20.347 276.321  57.297 19.774 383.542 3.825    259.046  5215.046
# 2 -37.024 27.708 59.162 39.533 15.331 96.217  6.835    101.324 4272.994  -27.455 59.497 241.521 110.603 41.787 350.658 4.932    320.232 10482.424
# 3 -48.507  9.638  9.660  9.649  3.782 20.937  5.275     23.605 1272.256  -92.606 16.090 111.171  36.363 13.288 156.023 8.568    123.163 14278.678
# 4  -9.857  3.183  3.189  3.186  1.249  6.914 13.253      7.795   85.374  -48.987 14.055 131.626  35.199 12.560 183.370 4.138    134.245  8424.807
# 5  -7.445  9.234 10.721  9.940  3.895 21.667  6.066     24.361  201.719   -4.769 34.549  39.853  37.075 14.528  80.780 8.501     90.847   481.862
# 6 -16.864  3.178 62.249  9.963  3.312 86.315  5.942     54.570 1231.948   -5.829 22.011 497.027  71.789 23.499 689.094 0.000    426.097  3347.283
# 7 -58.909  7.084 10.927  8.730  3.410 19.669  8.634     21.679 1430.995 -251.516 16.411 119.770  37.798 13.757 167.795 5.015    130.601 41302.221
# 8 -43.639 19.155 49.441 29.651 11.429 76.900  8.476     78.063 3930.240  -69.250 56.866 183.765  96.592 36.899 274.567 4.844    264.776 21525.759

# $NPY_Cre_X_iSPN_WT
#         A1 τrise τdecay tpeak r20_80  d80_20 delay half_width     area1      A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -24.765 3.283  6.864 4.641  1.801  11.226 5.335     11.864   334.248 -18.260  3.599 129.618  13.266  4.111 179.689  9.946    104.572  2621.938
# 2  -38.914 2.987  7.263 4.508  1.741  11.444 4.673     11.760   525.756 -29.206  1.450 133.800   6.632  1.817 185.487  2.680     99.892  4106.368
# 3  -44.072 3.819  3.824 3.822  1.498   8.293 4.043      9.349   457.844  -2.458 26.678  26.959  26.818 10.512  58.193 17.682     65.607   179.176
# 4  -32.985 2.757 21.545 6.500  2.354  30.127 4.951     23.047   960.950  -8.676 20.492 293.323  58.630 20.131 407.041  0.900    272.099  3107.883
# 5 -191.684 2.162 99.992 8.474  2.545 138.619 6.051     78.625 20861.978 -62.773 64.643 343.276 132.971 49.464 487.726  0.000    424.679 31742.637
# 6  -41.619 3.400 15.686 6.637  2.490  22.510 5.677     19.791   996.712 -22.638  2.196 181.160   9.809  2.728 251.140  4.084    136.170  4329.315

# $NPY_Cre_X_iSPN_6OHDA
#         A1 τrise  τdecay  tpeak r20_80  d80_20 delay half_width     area1       A1  τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area1
# 1  -45.758 3.080  91.172 10.799  3.427 126.393 4.577     75.291  4696.459  -28.192 63.403 326.170 128.905 48.046 464.351  1.962    405.532 13652.173
# 2 -183.839 1.711  10.275  3.681  1.358  14.502 5.179     11.879  2702.604  -96.653  1.130 118.579   5.309  1.429 164.385  3.787     87.898 11985.771
# 3   -9.302 6.473   6.485  6.479  2.540  14.059 4.383     15.850   163.830  -37.810 13.205 139.381  34.376 12.142 193.856  6.032    138.085  6744.163
# 4  -31.284 5.492   5.502  5.497  2.155  11.928 4.981     13.447   467.451  -49.662  3.613 151.573  13.829  4.205 210.125  3.365    120.321  8246.571
# 5 -159.597 5.225   5.242  5.234  2.052  11.357 2.746     12.804  2270.651  -56.981 10.220 138.260  28.746  9.924 191.913 13.640    129.718  9698.934
# 6 -234.010 1.833   8.699  3.617  1.355  12.456 4.161     10.869  3085.194 -113.134  1.137 135.556   5.483  1.450 187.920  5.550     99.837 15968.887
# 7  -42.266 1.971   9.667  3.936  1.471  13.807 4.337     11.941   613.927  -54.131  6.116 158.054  20.689  6.670 219.119  3.464    132.890  9752.172
# 8 -440.901 2.841  76.153  9.705  3.115 105.574 4.635     63.710 38139.375 -150.636 31.571 189.214  67.854 25.043 267.087  6.800    218.872 40796.846
# 9  -47.117 7.552 138.952 23.257  7.784 192.688 7.207    123.091  7739.849  -13.979 42.129 522.430 115.374 40.162 725.550  0.298    502.572  9107.897

setwd(svg_path)

############################################################## AMPLITUDE ############################################################## 

amplitude_dSPN_control <- create_df(matrix(peaks_list$NPY_Cre_X_dSPN_WT, ncol=1), levels = list(condition = 'control', cell_type = 'dSPN'))
amplitude_iSPN_control <- create_df(matrix(peaks_list$NPY_Cre_X_iSPN_WT,  ncol=1), levels = list(condition = 'control', cell_type = 'iSPN'), 
  start_id = length(unique(amplitude_dSPN_control$s)) + 1)
amplitude_dSPN_6OHDA <- create_df(matrix(peaks_list$NPY_Cre_X_dSPN_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'dSPN'),
  start_id = length(unique(amplitude_dSPN_control$s)) + length(unique(amplitude_iSPN_control$s)) + 1)
amplitude_iSPN_6OHDA <- create_df(matrix(peaks_list$NPY_Cre_X_iSPN_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'iSPN'),
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
yrange <- c(-500, 0)
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
  yrange=c(1,500),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_amplitude_SPN.svg', width=width, height=height, bg='transparent')


########################################################### CHARGE TRANSFER ###########################################################
area_dSPN_control <- create_df(matrix(areas_list$NPY_Cre_X_dSPN_WT, ncol=1), levels = list(condition = 'control', cell_type = 'dSPN'), var_name = 'charge_transfer')
area_iSPN_control <- create_df(matrix(areas_list$NPY_Cre_X_iSPN_WT,  ncol=1), levels = list(condition = 'control', cell_type = 'iSPN'), var_name = 'charge_transfer', 
  start_id = length(unique(area_dSPN_control$s)) + 1)
area_dSPN_6OHDA <- create_df(matrix(areas_list$NPY_Cre_X_dSPN_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'dSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN_control$s)) + length(unique(area_iSPN_control$s)) + 1)
area_iSPN_6OHDA <- create_df(matrix(areas_list$NPY_Cre_X_iSPN_6OHDA, ncol=1), levels = list(condition = '6OHDA', cell_type = 'iSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN_control$s)) + length(unique(area_iSPN_control$s)) + length(unique(area_dSPN_6OHDA$s)) + 1)

area_SPN <- rbind(area_dSPN_control, area_iSPN_control, area_dSPN_6OHDA, area_iSPN_6OHDA)

area_SPN$charge_transfer <- -area_SPN$charge_transfer

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- FALSE
yrange <- if (log_y) c(0.1, 100) else c(0, 100)
y_tick_interval <- 25

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_SPN[,2:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[5:8,], log_y=FALSE,)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_SPN[,2:4], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,100),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_dSPN <- create2condition_df(fits_list$NPY_Cre_X_dSPN_WT, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'))
amplitude_2component_iSPN <- create2condition_df(fits_list$NPY_Cre_X_iSPN_WT,  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + 1)
amplitude_2component_dSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_dSPN_6OHDA, levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'),
    start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1]  + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + 1)
amplitude_2component_iSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_iSPN_6OHDA,  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + dim(fits_list$NPY_Cre_X_iSPN_WT)[1]+ 1)

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
yrange <- c(-500, 0)
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

BoxPlot3(formula=amplitude ~ cell_type*condition, data=transform(amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='fast',c(3,2,5)], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,500),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[1:4,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_fast_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ cell_type*condition, data=amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='slow_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=amplitude ~ cell_type*condition, data=transform(amplitude_2component_SPN[amplitude_2component_SPN$kinetics=='slow',c(3,2,5)], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,500),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[5:8,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_slow_amplitude_2component_SPN.svg', width=width, height=height, bg='transparent')

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


# df <- amplitude_2component_SPN
# df1 <- df[df$condition=='ctrl',]
# df1 <- df1[order(df1$kinetics != 'fast', df1$s), ]

# df1$kinetics <- NULL
# df1$x <- ifelse(df1$cell_type == 'dSPN', 1, 2)
# df1$condition <- NULL
# names(df1)[names(df1) == 'amplitude'] <- 'y'
# df1 <- df1[, c('s', 'x', 'y')]
# df1$s <- as.integer(as.character(df1$s))


width=2.875
height=3.5
xlim=c(0, 500)
ylim=c(0, 500)
ScatterPlot(df1, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter_dSPN.svg', lwd=lwd, save=plotsave)

ScatterPlot(df1, sign=-1, xlim=c(1, 500), ylim=c(1, 500), height=height, width=width, filename='scatter_logxy_dSPN.svg', lwd=lwd, log_xy=TRUE, save=plotsave)


df <- amplitude_2component_SPN
df2 <- df[df$cell_type=='iSPN',]
df2 <- df2[order(df2$kinetics != 'fast', df2$s), ]


df2$kinetics <- NULL
df2$x <- ifelse(df2$condition == 'ctrl', 1, 2)
df2$condition <- NULL
names(df2)[names(df2) == 'amplitude'] <- 'y'
df2 <- df2[, c('s', 'x', 'y')]
df2$s <- as.integer(as.character(df2$s))

# df2 <- df[df$condition=='6OHDA',]
# df2 <- df2[order(df2$kinetics != 'fast', df2$s), ]

# df2$kinetics <- NULL
# df2$x <- ifelse(df2$cell_type == 'dSPN', 1, 2)
# df2$condition <- NULL
# names(df2)[names(df2) == 'amplitude'] <- 'y'
# df2 <- df2[, c('s', 'x', 'y')]
# df2$s <- as.integer(as.character(df2$s))


ScatterPlot(df2, sign=-1, xlim=xlim, ylim=ylim, height=height, width=width, open_symbols=TRUE, filename='scatter_iSPN.svg', lwd=lwd, save=plotsave)

ScatterPlot(df2, sign=-1, xlim=c(1, 500), ylim=c(1, 500), height=height, width=width, open_symbols=TRUE, filename='scatter_logxy_iSPN.svg', lwd=lwd, log_xy=TRUE, save=plotsave)

#################################################### CHARGE TRANSFER FAST vs SLOW ####################################################
area_2component_dSPN <- create2condition_df(fits_list$NPY_Cre_X_dSPN_WT, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='charge_transfer')
area_2component_iSPN <- create2condition_df(fits_list$NPY_Cre_X_iSPN_WT, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + 1)
area_2component_dSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_dSPN_6OHDA, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1]  + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + 1)
area_2component_iSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_iSPN_6OHDA, cols = c(9, 18),  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='charge_transfer',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + dim(fits_list$NPY_Cre_X_iSPN_WT)[1]+ 1)

area_2component_SPN <- rbind(
  add_condition_column(area_2component_dSPN, 'ctrl', 'charge_transfer'),
  add_condition_column(area_2component_iSPN, 'ctrl', 'charge_transfer'),
  add_condition_column(area_2component_dSPN_6OHDA, '6OHDA', 'charge_transfer'),
  add_condition_column(area_2component_iSPN_6OHDA, '6OHDA', 'charge_transfer')
  )
area_2component_SPN$charge_transfer <- area_2component_SPN$charge_transfer/1e3

# MCwilcox(formula=amplitude ~ condition*cell_type, df=area_2component_SPN[area_2component_SPN$kinetics=='fast',])
# MCwilcox(formula=amplitude ~ kinetics*cell_type + Error(s), df=area_2component_SPN[area_2component_SPN$condition=='ctrl',-3])

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_2component_SPN[area_2component_SPN$kinetics=='fast',]))
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=charge_transfer ~ cell_type*condition, df=area_2component_SPN[area_2component_SPN$kinetics=='slow',]))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 60)
y_tick_interval <- 20
width <- 3
height <- 3.5

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='fast',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[9:12,])
if (plotsave) save_graph(svg_path=svg_path, filename='fast_area_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='fast',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.01,50),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[9:12,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_fast_area_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[13:16,])
if (plotsave) save_graph(svg_path=svg_path, filename='slow_area_2component_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ cell_type*condition, data=area_2component_SPN[area_2component_SPN$kinetics=='slow',c(3,2,5)], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.01,50),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary1[13:16,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_slow_area_2component_SPN.svg', width=width, height=height, bg='transparent')

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
tau_2component_dSPN <- create2condition_df(fits_list$NPY_Cre_X_dSPN_WT, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau')
tau_2component_iSPN <- create2condition_df(fits_list$NPY_Cre_X_iSPN_WT, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='tau',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + 1)
tau_2component_dSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_dSPN_6OHDA, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1]  + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + 1)
tau_2component_iSPN_6OHDA <- create2condition_df(fits_list$NPY_Cre_X_iSPN_6OHDA, cols = c(3, 12),  levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'), var_name='tau',
  start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + dim(fits_list$NPY_Cre_X_dSPN_6OHDA)[1] + dim(fits_list$NPY_Cre_X_iSPN_WT)[1]+ 1)

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

# tau_2component_dSPN <- create2condition_df(fits_list$NPY_Cre_X_dSPN_WT, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'dSPN'), var_name='tau')
# tau_2component_iSPN <- create2condition_df(fits_list$NPY_Cre_X_iSPN_WT,  cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'iSPN'),
#   start_id = dim(fits_list$NPY_Cre_X_dSPN_WT)[1] + 1, var_name='tau')
# tau_2component_SPN <- rbind(tau_2component_dSPN, tau_2component_iSPN)

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ kinetics*condition + Error(s), tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')]))
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=tau ~ kinetics*condition + Error(s), tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')]))

# update graph properties
log_y <- FALSE
ylab <- expression(tau[decay] * ' ' * (ms))
yrange <- if (log_y) c(1, 500) else c(0,500)
y_tick_interval <- 100
width <- 3
height <- 3.5

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_dSPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 500) else c(0,500)

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='dSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[17:20,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_dSPN.svg', width=width, height=height, bg='transparent')


log_y <- FALSE
yrange <- if (log_y) c(1, 500) else c(0,500)
BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[21:24,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_2component_iSPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 500) else c(0,500)

BoxPlot3(formula=tau ~ kinetics*condition, data=tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[21:24,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_2component_iSPN.svg', width=width, height=height, bg='transparent')



# library(qboxplot)
# for visual quecking
# # Draw boxplots for groups using type 6 quantiles
# qboxplot(tau ~ kinetics*condition, data = tau_2component_SPN[tau_2component_SPN$cell_type=='iSPN', setdiff(names(tau_2component_SPN), 'cell_type')],
#          qtype = 6, range = 1.5, log='y',
#          main = "Boxplot using type 6 quantiles")


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
yrange <- if (log_y) c(1, 100) else c(0,100)
y_tick_interval <- 25

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='fast', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 100) else c(0,100)

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='fast', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, 
  lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_fast_2component_SPN.svg', width=width, height=height, bg='transparent')



log_y <- FALSE
ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
yrange <- if (log_y) c(1, 500) else c(0,500)
y_tick_interval <- 100

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='slow', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[13:16,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_slow_2component_SPN.svg', width=width, height=height, bg='transparent')

log_y <- TRUE
yrange <- if (log_y) c(1, 500) else c(0,500)

BoxPlot3(formula=tau ~ cell_type*condition, data=tau_2component_SPN[tau_2component_SPN$kinetics=='slow', setdiff(names(tau_2component_SPN), c('s', 'kinetics'))], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, 
  lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, test_result=stats_summary[9:12,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_tau_slow_2component_SPN.svg', width=width, height=height, bg='transparent')

# ########################################################## SINGLE EXAMPLES ###########################################################
# # single examples

for (iii in 1:length(output_list)) {
  name <- names(output_list)[iii]
  for (ii in 1:length(output_list[[iii]])){
    fit_plot(traces=output_list[[iii]][[ii]]$traces, func=product2N, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
      ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
  }
}
graphics.off()

# single egs
dx <- 0.1
xlab='time (ms)'
ylab='PSC amplitude (pA)'
lwd=1.0
filter='off'
xlim <- c(0, 600)
ylim <- c(-200, 10)
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 6
iSPNs_control_egs <- output_list$NPY_Cre_X_iSPN_WT[[ii]]$traces
single_fit_egs(traces=iSPNs_control_egs, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_iSPN.svg', save=plotsave)
single_fit_egs(traces=iSPNs_control_egs, xlim=c(100, 345), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_iSPN.svg', log_y=TRUE, save=plotsave)

ii = 7
iSPNs_6OHDA_egs <- output_list$NPY_Cre_X_iSPN_6OHDA[[ii]]$traces
single_fit_egs(traces=iSPNs_6OHDA_egs, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_iSPN_6OHDA.svg', save=plotsave)
single_fit_egs(traces=iSPNs_6OHDA_egs, xlim=c(100, 345), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_iSPN_6OHDA.svg', log_y=TRUE, save=plotsave)

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
# $NPY_Cre_X_dSPN_WT
#       Afast  τrise τdecay  tpeak r20_80 d80_20 delay half_width     area   Aslow   τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area
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

# $NPY_Cre_X_dSPN_6OHDA
#     Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area    Aslow  τrise  τdecay   tpeak r20_80  d80_20 delay half_width      area
# 1 -36.672  6.523 62.230 16.436  5.856 86.669  6.720     63.177 2971.913  -15.339 20.347 276.321  57.297 19.774 383.542 3.825    259.046  5215.046
# 2 -37.024 27.708 59.162 39.533 15.331 96.217  6.835    101.324 4272.994  -27.455 59.497 241.521 110.603 41.787 350.658 4.932    320.232 10482.424
# 3 -48.507  9.638  9.660  9.649  3.782 20.937  5.275     23.605 1272.256  -92.606 16.090 111.171  36.363 13.288 156.023 8.568    123.163 14278.678
# 4  -9.857  3.183  3.189  3.186  1.249  6.914 13.253      7.795   85.374  -48.987 14.055 131.626  35.199 12.560 183.370 4.138    134.245  8424.807
# 5  -7.445  9.234 10.721  9.940  3.895 21.667  6.066     24.361  201.719   -4.769 34.549  39.853  37.075 14.528  80.780 8.501     90.847   481.862
# 6 -16.864  3.178 62.249  9.963  3.312 86.315  5.942     54.570 1231.948   -5.829 22.011 497.027  71.789 23.499 689.094 0.000    426.097  3347.283
# 7 -58.909  7.084 10.927  8.730  3.410 19.669  8.634     21.679 1430.995 -251.516 16.411 119.770  37.798 13.757 167.795 5.015    130.601 41302.221
# 8 -43.639 19.155 49.441 29.651 11.429 76.900  8.476     78.063 3930.240  -69.250 56.866 183.765  96.592 36.899 274.567 4.844    264.776 21525.759

# $NPY_Cre_X_iSPN_WT
#      Afast τrise τdecay tpeak r20_80  d80_20 delay half_width      area   Aslow  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area
# 1  -24.765 3.283  6.864 4.641  1.801  11.226 5.335     11.864   334.248 -18.260  3.599 129.618  13.266  4.111 179.689  9.946    104.572  2621.938
# 2  -38.914 2.987  7.263 4.508  1.741  11.444 4.673     11.760   525.756 -29.206  1.450 133.800   6.632  1.817 185.487  2.680     99.892  4106.368
# 3  -44.072 3.819  3.824 3.822  1.498   8.293 4.043      9.349   457.844  -2.458 26.774  26.862  26.818 10.512  58.193 17.681     65.607   179.176
# 4  -32.985 2.757 21.545 6.500  2.354  30.127 4.951     23.047   960.950  -8.676 20.492 293.323  58.630 20.131 407.041  0.900    272.099  3107.883
# 5 -191.684 2.162 99.992 8.474  2.545 138.619 6.051     78.625 20861.978 -62.773 64.643 343.276 132.971 49.464 487.726  0.000    424.679 31742.637
# 6  -41.619 3.400 15.686 6.637  2.490  22.510 5.677     19.791   996.712 -22.638  2.196 181.160   9.809  2.728 251.140  4.084    136.170  4329.315

# $NPY_Cre_X_iSPN_6OHDA
#      Afast τrise  τdecay  tpeak r20_80  d80_20 delay half_width      area    Aslow  τrise  τdecay   tpeak r20_80  d80_20  delay half_width      area
# 1  -45.758 3.080  91.172 10.799  3.427 126.393 4.577     75.291  4696.459  -28.192 63.403 326.170 128.905 48.046 464.351  1.962    405.532 13652.173
# 2 -183.839 1.711  10.275  3.681  1.358  14.502 5.179     11.879  2702.604  -96.653  1.130 118.579   5.309  1.429 164.385  3.787     87.898 11985.771
# 3   -6.016 6.164   6.171  6.168  2.418  13.383 6.361     15.088   100.854  -38.008 13.531 139.373  34.949 12.371 193.902  4.385    138.876  6807.087
# 4  -31.284 5.492   5.502  5.497  2.155  11.928 4.981     13.447   467.453  -49.662  3.613 151.573  13.829  4.205 210.125  3.365    120.321  8246.571
# 5 -159.597 5.225   5.242  5.234  2.052  11.357 2.746     12.804  2270.651  -56.981 10.220 138.260  28.746  9.924 191.913 13.640    129.718  9698.934
# 6 -234.010 1.833   8.699  3.617  1.355  12.456 4.161     10.869  3085.194 -113.134  1.137 135.556   5.483  1.450 187.920  5.550     99.837 15968.887
# 7  -42.266 1.971   9.667  3.936  1.471  13.807 4.337     11.941   613.927  -54.131  6.116 158.054  20.689  6.670 219.119  3.464    132.890  9752.172
# 8 -440.901 2.841  76.153  9.705  3.115 105.574 4.635     63.710 38139.375 -150.636 31.571 189.214  67.854 25.043 267.087  6.800    218.872 40796.846
# 9  -47.117 7.552 138.952 23.257  7.784 192.688 7.207    123.091  7739.849  -13.979 42.129 522.430 115.374 40.162 725.550  0.298    502.572  9107.897

# save all to single 'xlsx'
if (plotsave){
  data_list <-
    c(fits_list, 
      list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'amplitude 2 components' = amplitude_2component_SPN,
        'tau 2 components' = tau_2component_SPN,
        'iSPN single examples control' = iSPNs_control_egs,
        'iSPN single examples 6OHDA' = iSPNs_6OHDA_egs,
        'statistics' = stats_summary,
        'additional statistics' = stats_summary1
      )
    )
  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.csv'), wd=xlsx_path)
}


