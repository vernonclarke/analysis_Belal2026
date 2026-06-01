rm(list = ls(all = TRUE))
graphics.off()

plotsave <- TRUE

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')

identifier <- 'Figure 11'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

setwd(analysis_path)

load(paste0(analysis_path, '/', identifier, '.RData'))


setwd(svg_path)

############################################################## AMPLITUDE ##############################################################
amplitude_ChI_NGF_control <- create_df(matrix(ChI_NGF_control_peaks, ncol=1), levels = list(condition = 'control'))
amplitude_ChI_NGF_6OHDA   <- create_df(matrix(ChI_NGF_6OHDA_peaks, ncol=1),   levels = list(condition = '6OHDA'),
  start_id = length(unique(amplitude_ChI_NGF_control$s)) + 1)

amplitude <- rbind(amplitude_ChI_NGF_control, amplitude_ChI_NGF_6OHDA)

# create output for stats
stats_summary <- MCwilcox(formula=amplitude ~ condition, df=amplitude)

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(PSC~amplitude~(pA))
xrange <- c(0.75, 4.25)
yrange <- c(-800, 0)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 200
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot3(formula=amplitude ~ condition, data=amplitude[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange,  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude.svg', width=width, height=height, bg='transparent')


BoxPlot3(formula=amplitude ~ condition, data=transform(amplitude[,2:3], amplitude = -1 * amplitude), wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(1,1000),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[1,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_amplitude_SPN.svg', width=width, height=height, bg='transparent')


########################################################### CHARGE TRANSFER ###########################################################
area_ChI_NGF_control <- create_df(matrix(ChI_NGF_control_areas, ncol=1), levels = list(condition = 'control'), var_name='charge_transfer')
area_ChI_NGF_6OHDA   <- create_df(matrix(ChI_NGF_6OHDA_areas, ncol=1),   levels = list(condition = '6OHDA'),   var_name='charge_transfer',
  start_id = length(unique(area_ChI_NGF_control$s)) + 1)

area <- rbind(area_ChI_NGF_control, area_ChI_NGF_6OHDA)

# create output for stats
stats_summary <- rbind(stats_summary, MCwilcox(formula=charge_transfer ~ condition, df=area))


# update graph properties
ylab <- expression(charge~transfer~(pC))
log_y <- FALSE
yrange <- if (log_y) c(1, 60) else c(0, 60)
y_tick_interval <- 20

BoxPlot3(formula=charge_transfer ~ condition, data=area[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, test_result=stats_summary[2,], log_y=log_y)
if (plotsave) save_graph(svg_path=svg_path, filename='charge_transfer.svg', width=width, height=height, bg='transparent')

BoxPlot3(formula=charge_transfer ~ condition, data=area[,2:3], wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=c(0.1,100),  xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width, log_y=TRUE, test_result=stats_summary[2,])
if (plotsave) save_graph(svg_path=svg_path, filename='semilog_charge_transfer_SPN.svg', width=width, height=height, bg='transparent')

# ########################################################## SINGLE EXAMPLES ###########################################################
# single examples

datasets <- c('ChI_NGF_control', 'ChI_NGF_6OHDA')
name <- datasets[1]

# single egs
lwd <- 1.0
width <- 3.5
height <- 4.5
ylim <- c(-200, 20)
ybar <- 50
xbar <- 50


for (ii in 2:dim(ChI_NGF_control_data)[2]){
  single_egs(traces=ChI_NGF_control_data[c(1,ii)], xlim=c(0,1000), ylim=NULL, lwd=lwd, width=width, height=height, 
    ybar=ybar, xbar=xbar, filename=paste0('ChI_NGF_control_', ii, '.svg'), save=plotsave)
}

for (ii in 2:dim(ChI_NGF_6OHDA_data)[2]){
  single_egs(traces=ChI_NGF_6OHDA_data[c(1,ii)], xlim=c(0,1000), ylim=NULL, lwd=lwd, width=width, height=height, 
    ybar=ybar, xbar=xbar, filename=paste0('ChI_NGF_6OHDA_', ii, '.svg'), save=plotsave)
} 



ii <- 3
ctrl_egs <- ChI_NGF_control_data[c(1,ii)]
single_egs(traces=ctrl_egs, xlim=c(0,1000), ylim=ylim, lwd=lwd, width=width, height=height, 
  ybar=ybar, xbar=xbar, filename='egs_ctrl.svg', save=plotsave)

single_egs(traces=ctrl_egs, xlim=c(100,300), ylim=NULL, lwd=lwd, height=height*0.645, width=width*0.61, 
  log_y=TRUE, ybar=ybar, xbar=xbar, filename='egs_semilog_ctrl.svg', save=plotsave)


ii <- 7
OHDA_egs <- ChI_NGF_6OHDA_data[c(1,ii)]
single_egs(traces=OHDA_egs, xlim=c(0,1000), ylim=ylim, lwd=lwd, width=width, height=height, 
  ybar=ybar, xbar=xbar, filename='egs_OHDA.svg', save=plotsave)

single_egs(traces=OHDA_egs, xlim=c(100,300), ylim=NULL, lwd=lwd, height=height*0.645, width=width*0.61, 
  log_y=TRUE, ybar=ybar, xbar=xbar, filename='egs_semilog_OHDA.svg', save=plotsave)


graphics.off()

# save all to single 'xlsx'
if (plotsave){
  data_list <-
      list(
        'amplitude'  = amplitude,
        'charge transfer'  = area,
        'ChI-NGF single examples' = ctrl_egs,
        'ChI-NGF 6OHDA single examples' = OHDA_egs,
        'statistics' = stats_summary
      )

  # save to excel spreadsheet
  list2excel(data_list, paste0(identifier, '.xlsx'), wd=xlsx_path)
}



############################################################### Vcmd ################################################################
file_path3 <- paste0(repo_root, '/Paper analysis/Raw ABF data summaries/', identifier)
setwd(file_path3)
expt_id <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)

# initial raw summaries; good for approx amplitudes and exact command potential
summary2 <- setNames(lapply(expt_id, function(expt) {
  folder_base <- file.path(file_path3, expt, 'xlsx')
  folders <- list.dirs(path = folder_base, full.names = FALSE, recursive = FALSE)
  
  if (length(folders) == 0) return(list())
  
  out_summary <- lapply(folders, function(folder) {
    folder_path <- file.path(folder_base, folder)
    out1 <- load_data2(wd = folder_path, name = folder, header = TRUE)
    out1$summary
  })
  names(out_summary) <- folders
  out_summary
}), expt_id)

Vcmd <-  setNames(lapply(1:length(summary2), function(iii){
  sapply(1:length(summary2[[iii]]), function(ii) mean(summary2[[iii]][[ii]][, 'holding.potential.(mV)']))
}), expt_id)


Vcmd_NDNF <- do.call(rbind, lapply(names(Vcmd), function(name) {
  data.frame(
    condition = if (grepl('control', name)) 'control' else '6OHDA',
    Vcmd = Vcmd[[name]]
  )
}))

Vcmd_NDNF <- Vcmd_NDNF[order(Vcmd_NDNF$condition, decreasing = TRUE), ]
rownames(Vcmd_NDNF) <- 1:dim(Vcmd_NDNF)[1]

Vcmd_NDNF$condition <- factor(Vcmd_NDNF$condition, levels = c("control", "6OHDA"))

width <- 3
height <- 3.5

# update graph properties
ylab <- expression(V[cmd]~(mV))
yrange <- c(-100, -60)
y_tick_interval <- 10


BoxPlot3(formula= Vcmd~condition, data=Vcmd_NDNF, wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, 
  yrange=yrange, xlabel_angle=45, tick_length=tick_length, y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
  p.cex=p.cex, height=height, width=width)
if (plotsave) save_graph(svg_path=svg_path, filename='Vcmd_SPN.svg', width=width, height=height, bg='transparent')
