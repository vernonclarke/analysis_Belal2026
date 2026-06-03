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

# settings
identifier <- 'Figure S5'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# path where all graphs are stored
svg_path <- paste0(analysis_path, '/svg')

# load RDATA

# names of datasets to be loaded
datasets <- c('ctrl',    # ctrl
              '6OHDA')   # OHDA

datasets2 <- c('ctrl', 'OHDA')

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

# outliers / cf dbscan results based on Afast and Aslow
# Use malanobis distance
# Peter J. Rousseeuw & Bert C. van Zomeren (1990)
# Unmasking Multivariate Outliers and Leverage Points
# Journal of the American Statistical Association, 85(411), 633–639.

amplitude_2component_ctrl <- create2condition_df(fits_list$ctrl, levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'))
amplitude_2component_OHDA <- create2condition_df(fits_list$OHDA,  levels = list(condition = c('fast', 'slow'), cell_type = 'OHDA'),
  start_id = dim(fits_list$ctrl)[1] + 1)
amplitude_2component_SPN <- rbind(amplitude_2component_ctrl, amplitude_2component_OHDA)

wide_df <- reshape(
  amplitude_2component_SPN,
  idvar = c('s', 'cell_type'),
  timevar = 'condition',
  direction = 'wide'
)
names(wide_df)[names(wide_df) == 'amplitude.fast'] <- 'Afast'
names(wide_df)[names(wide_df) == 'amplitude.slow'] <- 'Aslow'


mat <- -wide_df[, c('Afast', 'Aslow')]

mat <- cbind(
  wide_df["cell_type"],
  -wide_df[, -c(1,2)]
)
cutoff <- sqrt(qchisq(0.99999, df = 2)) # very conservative; makes sense for small samples

setwd(svg_path)
out <- mv_outliers(mat, method='MCD', alpha=0.5, quant=0.9999995, plot=TRUE, type=c('both'),  palette='Roma', 
  xlab=expression(A[fast]~'(pA)'), ylab=expression(A[slow]~'(pA)'), filename='outlier_plot.svg', save=plotsave)


outlier <- out[names(out)]
outliers <- wide_df[outlier, ]

outlier_id <- lapply( split(outliers, outliers$cell_type), function(x) match(x$s, wide_df$s[wide_df$cell_type == x$cell_type[1]]))
# ctrl: 4
# OHDA: 5, 22


fits_list   <- mapply(function(df, idx) df[-idx, ], fits_list,   outlier_id[names(fits_list)],   SIMPLIFY = FALSE)
peaks_list <- mapply(function(v, idx) v[-idx], peaks_list, outlier_id[names(peaks_list)], SIMPLIFY = FALSE)
areas_list <- mapply(function(v, idx) v[-idx], areas_list, outlier_id[names(areas_list)], SIMPLIFY = FALSE)
data_list   <- mapply(function(df, idx) df[-idx, ], data_list,   outlier_id[names(data_list)],   SIMPLIFY = FALSE)
output_list <- mapply(function(lst, idx) lst[!names(lst) %in% as.character(idx)],
                      output_list, outlier_id[names(output_list)], SIMPLIFY = FALSE)

# Remove the loaded objects
remove_loaded_objects(loaded_objects)

# Display the list of fits objects
fits_list

# $ctrl
#          A1  τrise τdecay  tpeak r20_80 d80_20  delay half_width     area1       A1   τrise  τdecay   tpeak r20_80  d80_20  delay half_width    area1
# 1  -220.164 16.024 38.401 24.034  9.287 60.705  1.730     62.565 15808.673 -116.954 150.694 216.575 179.670 70.256 400.291  0.000    491.657 58064.83
# 2  -209.979 10.303 44.102 19.548  7.364 63.698  0.000     57.143 14425.347 -164.028  19.247 214.875  51.005 17.926 298.679  0.000    210.130 44688.42
# 3  -153.456 16.487 21.373 18.719  7.329 41.172  0.000     46.041  7874.352  -80.148  17.072 444.969  57.885 18.642 616.882 24.795    373.687 40617.82
# 5  -329.530  1.941  8.654  3.741  1.406 12.454  9.936     11.052  4393.885 -311.263  12.565  80.539  27.659 10.161 113.351  3.938     91.250 35341.33
# 6  -212.335  1.003 21.455  3.223  1.062 29.748 10.265     18.547  5294.290 -119.149  13.879 138.767  35.507 12.601 193.133  1.974    139.266 21355.19
# 7  -333.218  0.752 16.981  2.454  0.803 23.543  9.822     14.559  6538.070 -265.299   8.818 104.141  23.785  8.317 144.684  3.021    100.565 34717.14
# 8  -150.467  2.876  2.877  2.876  1.128  6.241 13.433      7.037  1176.462 -217.347   9.073  61.668  20.387  7.459  86.598 11.030     68.650 18654.95
# 9  -270.866  2.367 11.037  4.639  1.739 15.825 13.274     13.874  4551.671 -166.315  22.162  56.258  34.065 13.138  87.810  6.373     89.449 17142.96
# 10  -75.272  5.621  5.739  5.680  2.226 12.325 10.300     13.895  1162.183 -126.549  10.868  71.323  24.123  8.847 100.286  5.082     80.233 12658.27
# 11 -256.262  1.273 10.077  3.014  1.090 14.085 10.123     10.740  3482.612 -173.385  13.947  77.059  29.108 10.801 109.231  5.240     91.484 19493.13
# 12 -361.732  1.359 39.082  4.729  1.506 54.180  9.513     32.394 15955.727 -129.936  14.823 120.452  35.414 12.783 168.268  3.385    127.480 21000.57
# 13 -139.589 20.377 20.389 20.383  7.990 44.230  5.867     49.866  7734.325 -316.957  16.438  69.854  31.101 11.722 100.963 29.696     90.749 34558.50
# 14  -75.658  2.138 38.211  6.529  2.192 52.991 11.322     34.018  3429.660  -40.238   5.414 234.149  20.876  6.323 324.600  4.683    185.307 10300.30
# 15 -139.108  9.176 40.393 17.596  6.619 58.193 13.286     51.815  8686.588 -128.696  16.741  94.550  35.219 13.050 133.863  4.411    111.455 17660.29
# 16 -286.536  1.326  5.111  2.416  0.915  7.462  9.621      6.883  2349.554 -295.767   6.831  74.922  18.002  6.337 104.162  6.983     73.567 28178.04
# 17 -195.485  1.177 13.549  3.149  1.104 18.828 10.228     13.158  3341.661  -74.968  15.546  73.877  30.687 11.493 105.766  6.127     92.256  8390.45
# 18 -182.153  0.931 27.935  3.276  1.038 38.726 10.700     23.029  5721.460  -85.369  25.288 138.186  52.565 19.518 196.002  0.269    164.644 17257.00
# 19 -319.336  3.805  3.823  3.814  1.495  8.276  9.881      9.331  3310.742 -419.952  11.176 117.487  29.057 10.267 163.413 16.423    116.506 63182.96
# 20  -48.782 12.026 34.055 19.352  7.433 52.016 15.465     51.744  2932.431  -48.029  27.361 148.057  56.671 21.056 210.122  6.654    176.971 10427.11
# 21 -106.675  2.691 42.824  7.947  2.700 59.403 11.893     38.928  5499.626  -37.134  11.076 235.267  35.519 11.706 326.195  0.000    203.589 10160.05
# 22 -182.692 10.688 37.403 18.744  7.135 55.258 11.747     52.243 11279.026 -124.471  30.808 121.050  56.552 21.402 176.356  6.960    161.835 24039.49

# $OHDA
#          A1  τrise τdecay  tpeak r20_80  d80_20  delay half_width    area1       A1  τrise  τdecay  tpeak r20_80  d80_20  delay half_width     area1
# 1   -40.324 44.139 45.554 44.839 17.576  97.317 18.701    109.703 4915.577  -92.461 30.138 302.045 77.160 27.376 420.366  0.000    302.953 36055.783
# 2   -84.684  3.479 27.945  8.280  2.991  39.049 18.425     29.669 3182.557 -212.180 14.423 397.579 49.635 15.882 551.176 11.409    331.375 95575.627
# 3     0.000     NA     NA     NA     NA      NA     NA         NA    0.000  -75.259 23.130  95.886 43.349 16.359 138.920  3.220    125.692 11341.010
# 4    -5.501  2.896  2.897  2.896  1.135   6.285  4.650      7.086   43.313    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 6   -44.356  3.731  3.733  3.732  1.463   8.099  4.863      9.131  450.016  -46.596 16.404 131.338 38.994 14.092 183.540  9.200    139.559  8235.392
# 7   -79.585  3.578 10.546  5.853  2.245  15.994  4.978     15.764 1462.065 -136.464 25.866  80.136 43.190 16.529 120.522 18.610    117.413 18746.101
# 8  -111.791  3.809  3.810  3.809  1.493   8.266  9.433      9.319 1157.552 -116.479 44.385  44.435 44.410 17.408  96.366  1.004    108.644 14061.301
# 9   -25.520  5.794  5.835  5.815  2.279  12.617  0.000     14.225  403.364  -41.088 17.105  62.083 30.436 11.564  91.289  5.918     85.512  4164.856
# 10  -18.417  6.678  6.725  6.701  2.627  14.542  4.433     16.394  335.484  -14.573 27.120  47.375 35.384 13.790  81.715  0.800     88.753  1456.979
# 11  -67.863 46.560 46.570 46.565 18.253 101.042  1.300    113.916 8589.927    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 12 -133.463  2.580 18.606  5.918  2.156  26.076 12.574     20.358 3413.101  -30.936 16.900 151.819 41.749 14.949 211.651  0.000    156.484  6183.354
# 13  -73.778  3.159  3.159  3.159  1.238   6.855  5.184      7.729  633.600 -105.299 11.634  79.918 26.239  9.592 112.183 10.317     88.689 11685.812
# 14  -28.310  9.000 30.132 15.507  5.915  44.788  5.121     42.817 1427.167  -22.239  7.095 177.310 23.788  7.698 245.816 21.671    149.780  4509.275
# 15  -33.940  4.944  4.953  4.949  1.940  10.738  4.505     12.106  456.553  -19.233 13.968 157.468 37.130 13.038 218.861 11.000    153.650  3833.854
# 16   -2.862  0.062 42.096  0.405  0.084  58.357  0.073     29.603  121.639  -17.141  3.828 468.389 18.552  4.888 649.326  9.156    344.540  8352.819
# 17    0.000     NA     NA     NA     NA      NA     NA         NA    0.000 -109.537  9.561 194.917 30.313 10.035 270.261 21.749    169.772 24943.173
# 18  -49.564  5.469 20.191  9.797  3.720  29.628  7.016     27.633 1625.776 -112.524  4.090  63.897 12.011  4.088  88.641 24.185     58.283  8676.800
# 19  -50.950  5.148  5.153  5.151  2.019  11.176  7.869     12.600  713.342 -266.232 10.961 112.875 28.309 10.020 157.038 19.423    112.479 38617.205
# 20  -66.809  1.588 11.944  3.695  1.342  16.719 23.830     12.913 1087.303  -41.555 23.958 132.212 49.979 18.546 187.424  7.033    157.022  8018.106
# 21    0.000     NA     NA     NA     NA      NA     NA         NA    0.000  -15.893  6.380 505.789 28.256  7.902 701.173  0.000    381.152  8500.350
# 23   -8.659  3.030  3.035  3.032  1.189   6.580 12.376      7.419   71.382  -34.446 14.560  65.250 28.113 10.564  93.859  6.320     83.174  3458.102

ctrl_names <- names(data_list$ctrl)[as.integer(rownames(fits_list$ctrl))]
OHDA_names <- names(data_list$OHDA)[as.integer(rownames(fits_list$OHDA))]

# ctrl_names <- abbreviate(colnames(data_list$ctrl), minlength = 8)
# OHDA_names <- abbreviate(colnames(data_list$OHDA), minlength = 8)

cat(sprintf("ctrl: n = %d cells from %d animals\n",
    length(ctrl_names),
    length(unique(regmatches(ctrl_names, regexpr("\\d{6,8}", ctrl_names))))))

# ctrl: n = 21 cells from 9 animals

cat(sprintf("6OHDA: n = %d cells from %d animals\n",
    length(OHDA_names),
    length(unique(regmatches(OHDA_names, regexpr("\\d{6,8}", OHDA_names))))))

# 6OHDA: n = 21 cells from 9 animals


############################################################## AMPLITUDE ############################################################## 

amplitude_ctrl <- create_df(matrix(peaks_list$ctrl, ncol=1), levels = list(cell_type = 'ctrl'))
amplitude_OHDA <- create_df(matrix(peaks_list$OHDA,  ncol=1), levels = list(cell_type = 'OHDA'), 
  start_id = length(unique(amplitude_ctrl$s)) + 1)
amplitude_SPN <- rbind(amplitude_ctrl, amplitude_OHDA)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ cell_type, df=amplitude_SPN)

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

BoxPlot(formula=amplitude ~ cell_type, data=amplitude_SPN[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot(formula=amplitude ~ cell_type, data=transform(amplitude_SPN[,2:3], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,700),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_amplitude_SPN.svg', width=width, height=height, bg='transparent')

########################################################### CHARGE TRANSFER ###########################################################

area_ctrl <- create_df(matrix(areas_list$ctrl, ncol=1), levels = list(cell_type = 'ctrl'), var_name = 'charge_transfer')
area_OHDA <- create_df(matrix(areas_list$OHDA,  ncol=1), levels = list(cell_type = 'OHDA'), var_name = 'charge_transfer',
  start_id = length(unique(area_ctrl$s)) + 1)
area_SPN <- rbind(area_ctrl, area_OHDA)


# stats OHDAs
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type, df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 150)
y_tick_interval <- 50

BoxPlot(formula=charge_transfer ~ cell_type, data=area_SPN[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot(formula=charge_transfer ~ cell_type, data=area_SPN[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,150),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

######################################################## AMPLITUDE FAST vs SLOW ####################################################### 
amplitude_2component_ctrl <- create2condition_df(fits_list$ctrl, levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'))
amplitude_2component_OHDA <- create2condition_df(fits_list$OHDA,  levels = list(condition = c('fast', 'slow'), cell_type = 'OHDA'),
  start_id = dim(fits_list$ctrl)[1] + 1)
amplitude_2component_SPN <- rbind(amplitude_2component_ctrl, amplitude_2component_OHDA)

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ condition*cell_type + Error(s), df=amplitude_2component_SPN)

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
xlim=c(0, 500)
ylim=c(0, 500)
scatter_plot(wide_df, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter.svg', lwd=lwd, save=plotsave)

##################################################### CHARGE TRANSFER FAST vs SLOW ####################################################

area_2component_ctrl <- create2condition_df(fits_list$ctrl, cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'), var_name='area')
area_2component_OHDA <- create2condition_df(fits_list$OHDA,  cols = c(9, 18), levels = list(condition = c('fast', 'slow'), cell_type = 'OHDA'),
  start_id = dim(fits_list$ctrl)[1] + 1, var_name='area')
area_2component_SPN <- rbind(area_2component_ctrl, area_2component_OHDA)
area_2component_SPN$area <- area_2component_SPN$area/1e3

# create output for stats
stats_summary1 <- rbind(stats_summary1, MCwilcox(formula=area ~ condition*cell_type + Error(s), df=area_2component_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- TRUE
yrange <- if (log_y) c(0.01, 200) else c(0, 100)
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
xlim=c(0, 20)
ylim=c(0, 60)
scatter_plot(wide_df, xlim = xlim, ylim = ylim, height = height, width = width, 
  xlab = expression(charge~transfer[fast]*' '(pC)), ylab = expression(charge~transfer[slow]*' '(pC)),
  x_tick_interval=4, y_tick_interval=10,
  filename = 'scatter_charge_transfer.svg', lwd = lwd, save = plotsave)

################################################### DECAY TIME CONSTANT FAST vs SLOW ##################################################

tau_2component_ctrl <- create2condition_df(fits_list$ctrl, cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'ctrl'), var_name='tau')
tau_2component_OHDA <- create2condition_df(fits_list$OHDA,  cols = c(3, 12), levels = list(condition = c('fast', 'slow'), cell_type = 'OHDA'),
  start_id = dim(fits_list$ctrl)[1] + 1, var_name='tau')
tau_2component_SPN <- rbind(tau_2component_ctrl, tau_2component_OHDA)

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
yrange <- if (log_y) c(1, 50) else c(0, 50)
y_tick_interval <- 10

BoxPlot(formula=tau ~ cell_type, data=tau_2component_SPN[tau_2component_SPN$condition == 'fast',3:4], 
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, 
  test_result=stats_summary[3,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='tau_fast_SPN.svg', width=width, height=height, bg='transparent')


ylab <- expression(tau[decay * ',' * slow] * ' ' * (ms))
log_y <- FALSE
yrange <- if (log_y) c(1, 300) else c(0, 300)
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
  traces <- output_list$ctrl[[ii]]$traces
  func <- if (dim(traces)[2] == 4) product1N else product2N
  fit_plot(traces=traces, func=func, filename=paste0(name, '_', ii, '.svg'), xlab='time (ms)', 
    ylab='PSC amplitude (pA)', xlim=c(0,800), ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=plotsave)
}

name <- datasets[2]
for (ii in 1:length(output_list$OHDA)){
  traces <- output_list$OHDA[[ii]]$traces
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
xlim <- c(50, 400)
ylim <- c(-600, 20)
width <- 3.5
height <- 4.5
ybar <- 100
xbar <- 50
colors <- c('#4C77BB', '#5A9B79', '#F28E2B')

ii = 15
ctrls_eg <- output_list$ctrl[[ii]]$traces
single_fit_egs(traces=ctrls_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_ctrl.svg', save=plotsave)
single_fit_egs(traces=ctrls_eg, xlim=c(100, 400), ylim=NULL, lwd=lwd, colors=colors,
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_ctrl.svg', log_y=TRUE, save=plotsave)

ii = 14
OHDAs_eg <- output_list$OHDA[[ii]]$traces
single_fit_egs(traces=OHDAs_eg, xlim=xlim, ylim=ylim, lwd=lwd, colors=colors, 
  height=height, width=width, xbar=xbar, ybar=ybar, filename='egs_OHDA.svg', save=plotsave)
single_fit_egs(traces=OHDAs_eg, xlim=c(100, 400), ylim=NULL, lwd=lwd, colors=colors, 
  height=height*0.645, width=width*0.61, xbar=xbar, ybar=ybar, filename='semilog_egs_OHDA.svg', log_y=TRUE, save=plotsave)

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
#       Afast  τrise τdecay  tpeak r20_80 d80_20  delay half_width      area    Aslow   τrise  τdecay   tpeak r20_80  d80_20  delay half_width     area
# 1  -220.164 16.024 38.401 24.034  9.287 60.705  1.730     62.565 15808.673 -116.954 150.694 216.575 179.670 70.256 400.291  0.000    491.657 58064.83
# 2  -209.979 10.303 44.102 19.548  7.364 63.698  0.000     57.143 14425.347 -164.028  19.247 214.875  51.005 17.926 298.679  0.000    210.130 44688.42
# 3  -153.456 16.487 21.373 18.719  7.329 41.172  0.000     46.041  7874.352  -80.148  17.072 444.969  57.885 18.642 616.882 24.795    373.687 40617.82
# 5  -329.530  1.941  8.654  3.741  1.406 12.454  9.936     11.052  4393.885 -311.263  12.565  80.539  27.659 10.161 113.351  3.938     91.250 35341.33
# 6  -212.335  1.003 21.455  3.223  1.062 29.748 10.265     18.547  5294.290 -119.149  13.879 138.767  35.507 12.601 193.133  1.974    139.266 21355.19
# 7  -333.218  0.752 16.981  2.454  0.803 23.543  9.822     14.559  6538.070 -265.299   8.818 104.141  23.785  8.317 144.684  3.021    100.565 34717.14
# 8  -150.467  2.876  2.877  2.876  1.128  6.241 13.433      7.037  1176.462 -217.347   9.073  61.668  20.387  7.459  86.598 11.030     68.650 18654.95
# 9  -270.866  2.367 11.037  4.639  1.739 15.825 13.274     13.874  4551.671 -166.315  22.162  56.258  34.065 13.138  87.810  6.373     89.449 17142.96
# 10  -75.272  5.621  5.739  5.680  2.226 12.325 10.300     13.895  1162.183 -126.549  10.868  71.323  24.123  8.847 100.286  5.082     80.233 12658.27
# 11 -256.262  1.273 10.077  3.014  1.090 14.085 10.123     10.740  3482.612 -173.385  13.947  77.059  29.108 10.801 109.231  5.240     91.484 19493.13
# 12 -361.732  1.359 39.082  4.729  1.506 54.180  9.513     32.394 15955.727 -129.936  14.823 120.452  35.414 12.783 168.268  3.385    127.480 21000.57
# 13 -139.589 20.377 20.389 20.383  7.990 44.230  5.867     49.866  7734.325 -316.957  16.438  69.854  31.101 11.722 100.963 29.696     90.749 34558.50
# 14  -75.658  2.138 38.211  6.529  2.192 52.991 11.322     34.018  3429.660  -40.238   5.414 234.149  20.876  6.323 324.600  4.683    185.307 10300.30
# 15 -139.108  9.176 40.393 17.596  6.619 58.193 13.286     51.815  8686.588 -128.696  16.741  94.550  35.219 13.050 133.863  4.411    111.455 17660.29
# 16 -286.536  1.326  5.111  2.416  0.915  7.462  9.621      6.883  2349.554 -295.767   6.831  74.922  18.002  6.337 104.162  6.983     73.567 28178.04
# 17 -195.485  1.177 13.549  3.149  1.104 18.828 10.228     13.158  3341.661  -74.968  15.546  73.877  30.687 11.493 105.766  6.127     92.256  8390.45
# 18 -182.153  0.931 27.935  3.276  1.038 38.726 10.700     23.029  5721.460  -85.369  25.288 138.186  52.565 19.518 196.002  0.269    164.644 17257.00
# 19 -319.336  3.805  3.823  3.814  1.495  8.276  9.881      9.331  3310.742 -419.952  11.176 117.487  29.057 10.267 163.413 16.423    116.506 63182.96
# 20  -48.782 12.026 34.055 19.352  7.433 52.016 15.465     51.744  2932.431  -48.029  27.361 148.057  56.671 21.056 210.122  6.654    176.971 10427.11
# 21 -106.675  2.691 42.824  7.947  2.700 59.403 11.893     38.928  5499.626  -37.134  11.076 235.267  35.519 11.706 326.195  0.000    203.589 10160.05
# 22 -182.692 10.688 37.403 18.744  7.135 55.258 11.747     52.243 11279.026 -124.471  30.808 121.050  56.552 21.402 176.356  6.960    161.835 24039.49

# $OHDA
#       Afast  τrise τdecay  tpeak r20_80  d80_20  delay half_width     area    Aslow  τrise  τdecay  tpeak r20_80  d80_20  delay half_width      area
# 1   -40.324 44.139 45.554 44.839 17.576  97.317 18.701    109.703 4915.577  -92.461 30.138 302.045 77.160 27.376 420.366  0.000    302.953 36055.783
# 2   -84.684  3.479 27.945  8.280  2.991  39.049 18.425     29.669 3182.557 -212.180 14.423 397.579 49.635 15.882 551.176 11.409    331.375 95575.627
# 3     0.000     NA     NA     NA     NA      NA     NA         NA    0.000  -75.259 23.130  95.886 43.349 16.359 138.920  3.220    125.692 11341.010
# 4    -5.501  2.896  2.897  2.896  1.135   6.285  4.650      7.086   43.313    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 6   -44.356  3.731  3.733  3.732  1.463   8.099  4.863      9.131  450.016  -46.596 16.404 131.338 38.994 14.092 183.540  9.200    139.559  8235.392
# 7   -79.585  3.578 10.546  5.853  2.245  15.994  4.978     15.764 1462.065 -136.464 25.866  80.136 43.190 16.529 120.522 18.610    117.413 18746.101
# 8  -111.791  3.809  3.810  3.809  1.493   8.266  9.433      9.319 1157.552 -116.479 44.385  44.435 44.410 17.408  96.366  1.004    108.644 14061.301
# 9   -25.520  5.794  5.835  5.815  2.279  12.617  0.000     14.225  403.364  -41.088 17.105  62.083 30.436 11.564  91.289  5.918     85.512  4164.856
# 10  -18.417  6.678  6.725  6.701  2.627  14.542  4.433     16.394  335.484  -14.573 27.120  47.375 35.384 13.790  81.715  0.800     88.753  1456.979
# 11  -67.863 46.560 46.570 46.565 18.253 101.042  1.300    113.916 8589.927    0.000     NA      NA     NA     NA      NA     NA         NA     0.000
# 12 -133.463  2.580 18.606  5.918  2.156  26.076 12.574     20.358 3413.101  -30.936 16.900 151.819 41.749 14.949 211.651  0.000    156.484  6183.354
# 13  -73.778  3.159  3.159  3.159  1.238   6.855  5.184      7.729  633.600 -105.299 11.634  79.918 26.239  9.592 112.183 10.317     88.689 11685.812
# 14  -28.310  9.000 30.132 15.507  5.915  44.788  5.121     42.817 1427.167  -22.239  7.095 177.310 23.788  7.698 245.816 21.671    149.780  4509.275
# 15  -33.940  4.944  4.953  4.949  1.940  10.738  4.505     12.106  456.553  -19.233 13.968 157.468 37.130 13.038 218.861 11.000    153.650  3833.854
# 16   -2.862  0.062 42.096  0.405  0.084  58.357  0.073     29.603  121.639  -17.141  3.828 468.389 18.552  4.888 649.326  9.156    344.540  8352.819
# 17    0.000     NA     NA     NA     NA      NA     NA         NA    0.000 -109.537  9.561 194.917 30.313 10.035 270.261 21.749    169.772 24943.173
# 18  -49.564  5.469 20.191  9.797  3.720  29.628  7.016     27.633 1625.776 -112.524  4.090  63.897 12.011  4.088  88.641 24.185     58.283  8676.800
# 19  -50.950  5.148  5.153  5.151  2.019  11.176  7.869     12.600  713.342 -266.232 10.961 112.875 28.309 10.020 157.038 19.423    112.479 38617.205
# 20  -66.809  1.588 11.944  3.695  1.342  16.719 23.830     12.913 1087.303  -41.555 23.958 132.212 49.979 18.546 187.424  7.033    157.022  8018.106
# 21    0.000     NA     NA     NA     NA      NA     NA         NA    0.000  -15.893  6.380 505.789 28.256  7.902 701.173  0.000    381.152  8500.350
# 23   -8.659  3.030  3.035  3.032  1.189   6.580 12.376      7.419   71.382  -34.446 14.560  65.250 28.113 10.564  93.859  6.320     83.174  3458.102

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
        'OHDA single examples' = OHDAs_eg,
        'statistics' = stats_summary
      )
    )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.csv'), wd=xlsx_path)
}


############################################ extra plots using separate analysis method  ############################################

amplitude_2component_SPN_alt <- load_data2(wd=xlsx_path, name='alternative_method')[[1]]

amplitude_2component_SPN_alt$amplitude <- amplitude_2component_SPN_alt$amplitude * -1

amplitude_2component_SPN_alt <- amplitude_2component_SPN_alt[
  !(paste(amplitude_2component_SPN_alt$cell_type, amplitude_2component_SPN_alt$s) %in% paste(outliers$cell_type, outliers$s)
  ),
]

rownames(amplitude_2component_SPN_alt) <- 1:dim(amplitude_2component_SPN_alt)[1]

# create output for stats
stats_summary1 <- MCwilcox(formula=amplitude ~ cell_type*condition + Error(s), df=amplitude_2component_SPN_alt)

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

BoxPlot(formula=amplitude ~ cell_type*condition, data=amplitude_2component_SPN_alt, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary1[1:4,])

if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_2component_SPN_alt.svg', width=width, height=height, bg='transparent')

# paired scatter plot Afast vs Aslow
wide_df <- reshape(
  amplitude_2component_SPN_alt,
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
xlim=c(0, 600)
ylim=c(0, 600)
scatter_plot(wide_df, xlim=xlim, ylim=ylim, height=height, width=width, filename='scatter_alt.svg', lwd=lwd, save=plotsave)


if (plotsave){
  data_list2 <-
    c(list(
        'amplitude 2 components' = amplitude_2component_SPN_alt,
        'ctrl single examples' = ctrls_eg,
        'OHDA single examples' = OHDAs_eg,
        'statistics' = stats_summary1
      )
    )

  # save to excel spreadsheet
  list2excel(data_list2, paste0(identifier, '_alt.xlsx'), wd=xlsx_path)
  list2csv(data_list2, paste0(identifier, '_alt.csv'), wd=xlsx_path)
}



