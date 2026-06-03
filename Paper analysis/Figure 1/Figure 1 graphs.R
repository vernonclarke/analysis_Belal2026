# =========================================# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- TRUE

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') Sys.getenv('USERPROFILE') else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))

identifier <- 'Figure 1'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Settings
# load RDATA
load(file.path(analysis_path, paste0(identifier, '.RData')))

setwd(svg_path)

# only GABAzine
dSPN_amplitude <- dSPN_peaks[, c(1,2)]
dSPN_area <- dSPN_areas[, c(1,2)]

iSPN_amplitude <- iSPN_peaks[, c(1,2)]
iSPN_area <- iSPN_areas[, c(1,2)]

############################################################### DBscan ################################################################ 

# DBSCAN_analyse(data=dSPN_area[,c(1,2)], minPts=5) # eps=50 / 100; at least 5 points in cluster
DBSCAN_analyse(data=dSPN_areas[,c(1,2)], minPts=5, eps=30, filename='dbscan_dSPN_svg', save=plotsave)
DBSCAN_analyse(data=iSPN_areas[,c(1,2)], minPts=5, eps=40, filename='dbscan_iSPN_svg', save=plotsave)

############################################################## AMPLITUDE ############################################################## 

amplitude_dSPN <- create_df(dSPN_amplitude, levels = list(condition = c('ctrl', 'GABAzine'), cell_type = 'dSPN'))
amplitude_iSPN <- create_df(iSPN_amplitude, levels = list(condition = c('ctrl', 'GABAzine'), cell_type = 'iSPN'), 
  start_id = length(unique(amplitude_dSPN$s)) + 1)
amplitude_SPN <- rbind(amplitude_dSPN, amplitude_iSPN)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ cell_type*condition + Error(s), df=amplitude_SPN)

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

BoxPlot(formula=amplitude ~ condition*cell_type + Error(s), data=amplitude_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[1:2,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')


BoxPlot(formula=amplitude ~ condition*cell_type + Error(s), data=transform(amplitude_SPN, amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,400),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1:2,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_semilog_SPN.svg', width=width, height=height, bg='transparent')

########################################################### CHARGE TRANSFER ###########################################################

area_dSPN <- create_df(dSPN_area, levels = list(condition = c('ctrl', 'GABAzine'), cell_type = 'dSPN'), var_name = 'charge_transfer')
area_iSPN <- create_df(iSPN_area, levels = list(condition = c('ctrl', 'GABAzine'), cell_type = 'iSPN'), var_name = 'charge_transfer',
  start_id = length(unique(area_dSPN$s)) + 1)
area_SPN <- rbind(area_dSPN, area_iSPN)

# stats tests
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ cell_type*condition + Error(s), df=area_SPN))

# update graph properties
ylab <- expression(charge~transfer~(pC))
yrange <- c(0, 80)
y_tick_interval <- 20

BoxPlot(formula=charge_transfer ~ condition*cell_type + Error(s), data=area_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[5:6,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

BoxPlot(formula=charge_transfer ~ condition*cell_type + Error(s), data=area_SPN, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,80),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1:2,])
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer_semilog_SPN.svg', width=width, height=height, bg='transparent')

########################################################### SINGLE EXAMPLES ###########################################################

# single egs
lwd <- 1.0
width <- 3.5
height <- 4.5
ybar <- 50
xbar <- 50

ii <- 10
dSPNs_egs <- dSPN_data[[ii]][, colnames(dSPN_data[[ii]]) %in% c('time', 'control', 'GABAzine')]
single_egs(traces=dSPNs_egs, xlim=c(50,800), lwd=lwd, width=width, height=height, 
  ybar=ybar, xbar=xbar, filename='egs_dSPN.svg', save=plotsave)
single_egs(traces=dSPNs_egs, xlim=c(100,600), lwd=lwd, height=height*0.645, width=width*0.61, log_y=TRUE, 
  ybar=ybar, xbar=xbar, filename='egs_semilog_dSPN.svg', save=plotsave)

ii <- 4
iSPNs_egs <- iSPN_data[[ii]][, colnames(iSPN_data[[ii]]) %in% c('time', 'control', 'GABAzine')]
single_egs(traces=iSPNs_egs, xlim=c(50,800), lwd=lwd, width=width, height=height, 
  ybar=ybar, xbar=xbar, filename='egs_iSPN.svg', save=plotsave)
single_egs(traces=iSPNs_egs, xlim=c(100,600), lwd=lwd, height=height*0.645, width=width*0.61, log_y=TRUE, 
  ybar=ybar, xbar=xbar, filename='egs_semilog_iSPN.svg', save=plotsave)

############################################################ SAVE OUTPUTS #############################################################

# save to one 'xlsx' file AND as individual 'csv' files
setwd(analysis_path)
# save all to single 'xlsx'
if (plotsave){
  data_list <-
    list(
        'amplitude'  = amplitude_SPN,
        'charge transfer'  = area_SPN,
        'dSPN single eg'  = dSPNs_egs,
        'iSPN single eg'  = iSPNs_egs,
        'statistics' = stats_summary
    )
  # save to excel and csv spreadsheets
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
  list2csv(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
}
