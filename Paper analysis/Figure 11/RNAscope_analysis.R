# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

root_dir <- Sys.getenv('ANALYSIS_ROOT', unset = file.path(if (.Platform$OS.type == 'windows') { user_profile <- Sys.getenv('USERPROFILE'); if (nzchar(user_profile)) user_profile else file.path('C:/Users', Sys.getenv('USERNAME')) } else path.expand('~'), 'Documents', 'Repositories', 'analysis_Belal2026'))

source(file.path(root_dir, 'R functions', 'setup.R'))
load_required_packages(c('bayesplot', 'brms', 'lme4', 'parallel', 'robustlmm', 'sciplot'))

identifier <- 'Figure 11'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# choose cell type: 'NDNF' or 'TH'
cell_type <- 'NDNF'
cell_type <- 'TH'

if (cell_type == 'NDNF') {
  cell_marker <- 'NDNF+'
  df_name <- 'df_NDNF'
  ylab <- expression('CHRNB2 Particles/NDNF Cell')
  yrange <- c(0, 60)
  line_ylim <- c(16, 26)
  dens_ylim <- c(0, 0.6)
  pred_ylim <- c(0, 0.08)
} else if (cell_type == 'TH') {
  cell_marker <- 'TH+'
  df_name <- 'df_TH'
  ylab <- expression('CHRNB2 Particles/TH Cell')
  yrange <- c(0, 20)
  line_ylim <- c(4, 12)
  dens_ylim <- c(0, 1.2)
  pred_ylim <- c(0, 0.16)
} else {
  stop("cell_type must be 'NDNF' or 'TH'")
}

nboot <- 1999
plotsave <- TRUE

# Load data
reanalysis = FALSE

if (reanalysis){
  RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_reanalysed')
}else{
  RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_experimenter')
}

RNAscope_data  <- RNAscope_data[order(RNAscope_data$condition,
                                      RNAscope_data$cell_type,
                                      RNAscope_data$hemisphere,
                                      RNAscope_data$slice_id), ]

rownames(RNAscope_data) = 1:dim(RNAscope_data)[1]

cell_data <- subset(RNAscope_data, cell_type == cell_marker)
stats_summary <- MCwilcox(formula=count ~ condition, df=cell_data)

df <- subset(RNAscope_data, cell_type == cell_marker)
df$Animal <- sub('^(L\\d+).*', '\\1', as.character(df$slice_id))
df$Group <- ifelse(df$hemisphere == 'UL', 'Control', '6OHDA')

df$Animal <- factor(df$Animal)
df$SliceID <- factor(df$slice_id)
df$Group <- factor(df$Group, levels = c('Control', '6OHDA'))
df$var <- df$count

options(contrasts=c('contr.sum', 'contr.poly'))

# Fit mixed-effects model
model <- lmerTest::lmer(var ~ Group + (1 | Animal/field), data = df)
summary(model)

# extract lmer fixed effects summary in csv-friendly format
lmer_sum <- summary(model)
fe <- as.data.frame(coef(summary(model)))
fe$parameter <- rownames(fe)
fe$model <- 'lmer: var ~ Group + (1 | Animal/field)'
fe$data <- df_name
rownames(fe) <- NULL
names(fe) <- c('estimate', 'se', 'df', 't_value', 'p_value', 'parameter', 'model', 'data')
fe <- fe[, c('model', 'data', 'parameter', 'estimate', 'se', 'df', 't_value', 'p_value')]

# Bootstrap
name <- paste0('RNAscope_', cell_type)
RData_path <- file.path(analysis_path, paste0('/', name, '.RData'))

if (file.exists(RData_path)) {
  load(RData_path)

} else {

  df$Animal  <- factor(df$Animal)
  df$field   <- factor(df$field)
  df$Group   <- factor(df$Group, levels = c('Control', '6OHDA'))

  rmod <- rlmer(var ~ Group + (1 | Animal/field), data = df)

  set.seed(42)
  animal_list  <- split(df, df$Animal)
  n_animals    <- length(animal_list)
  coef_name    <- grep('^Group', names(fixef(rmod)), value = TRUE)[1]
  group_levels <- levels(df$Group)

  est <- numeric(nboot)
  failed <- logical(nboot)
  dropped_one_group <- logical(nboot)

  for (ii in seq_len(nboot)) {
    if (ii %% 500 == 0) cat('iter', ii, 'of', nboot, '\n')

    picked_idx <- sample(seq_along(animal_list), size = n_animals, replace = TRUE)

    picked <- lapply(seq_along(picked_idx), function(j) {
      d <- animal_list[[picked_idx[j]]]
      d$Animal <- paste0(as.character(d$Animal[1]), '_boot', ii, '_', j)
      d$field  <- paste0(as.character(d$field), '_boot', ii, '_', j)
      d
    })

    dboot <- do.call(rbind, picked)
    dboot$Animal <- factor(dboot$Animal)
    dboot$field  <- factor(dboot$field)
    dboot$Group  <- factor(dboot$Group, levels = group_levels)

    if (length(unique(dboot$Group)) < 2) {
      est[ii] <- NA_real_
      dropped_one_group[ii] <- TRUE
      next
    }

    # Inference comes from the bootstrap CI/p-value and the Bayesian posterior,
    # both of which appropriately account for this. Singular-fit warnings during
    # bootstrap resampling are expected and can be ignored unless fits fail.

    fit <- tryCatch(
      suppressMessages(suppressWarnings(
        withCallingHandlers(
          rlmer(var ~ Group + (1 | Animal/field), data = dboot),
          warning = function(w) invokeRestart('muffleWarning'),
          message = function(m) invokeRestart('muffleMessage')
        )
      )),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      est[ii] <- NA_real_
      failed[ii] <- TRUE
    } else {
      est[ii] <- unname(fixef(fit)[coef_name])
    }
  }

  cat('N successful:', sum(is.finite(est)), 'out of', nboot, '\n')

  est <- est[is.finite(est)]
  if (!length(est)) stop('All bootstrap fits failed')
  ci <- quantile(est, probs = c(0.025, 0.975))
  names(ci) <- c('2.5%', '97.5%')
  p_lt0 <- mean(est < 0)
  pval  <- 2 * min(p_lt0, 1 - p_lt0)

  boot_summary <- data.frame(
    model       = 'rlmer bootstrap: var ~ Group + (1 | Animal/field)',
    data        = df_name,
    parameter   = coef_name,
    estimate    = unname(fixef(rmod)[coef_name]),
    p_value     = pval,
    ci_2.5      = unname(ci['2.5%']),
    ci_97.5     = unname(ci['97.5%']),
    n_boot      = nboot
  )
  rownames(boot_summary) <- NULL

  # Bayesian analysis
  options(mc.cores = max(1, parallel::detectCores() - 1))

  get_prior(bf(var ~ Group + (1 | Animal/field), center = TRUE), data = df)

  m1.prior <- c(
    prior(normal(0, 20), class = 'b')
  )

  mod <- var ~ Group + (1 | Animal/field)
  rstan:::rstan_options(disable_march_warning = TRUE)

  m1 <- brm(
    bf(var ~ Group + (1 | Animal/field)),
    data = df,
    family = student(),
    prior = c(
      prior(normal(0, 20), class = 'b'),
      prior(student_t(3, 0, 2.5), class = 'sigma')
    ),
    iter = 4000,
    chains = 4,
    seed = 42,
    control = list(adapt_delta = 0.99)
  )

  coef_name <- grep('^b_Group', names(as_draws_df(m1)), value = TRUE)[1]
  ps <- posterior_summary(m1, variable = coef_name)
  draws <- as_draws_df(m1)[[coef_name]]
  p_pos <- mean(draws < 0)

  post_diff <- -draws
  p1 <- mean(post_diff < 0)
  p2 <- 2 * min(mean(post_diff < 0), mean(post_diff > 0))

  new.data <- data.frame(
    field  = factor(c('new_control_field', 'new_6ohda_field')),
    Animal = factor(c('new_animal',        'new_animal')),
    Group  = factor(c('Control', '6OHDA'), levels = c('Control', '6OHDA'))
  )

  preds <- posterior_predict(
    m1,
    newdata    = new.data,
    re_formula = ~ (1 | Animal/field),
    allow_new_levels = TRUE
  )

  p3 <- mean(
    df$var[df$Group == 'Control'] > quantile(preds[, 1], 0.025) &
    df$var[df$Group == 'Control'] < quantile(preds[, 1], 0.975)
  )

  p4 <- mean(
    df$var[df$Group == '6OHDA'] > quantile(preds[, 2], 0.025) &
    df$var[df$Group == '6OHDA'] < quantile(preds[, 2], 0.975)
  )

  p5 <- mean(preds[, 2] < preds[, 1])

  bayes_summary <- data.frame(
    model                               = 'brm: var ~ Group + (1 | Animal/field), family = student()',
    data                                = df_name,
    parameter                           = coef_name,
    estimate                            = ps[1, 'Estimate'],
    ci_2.5                              = ps[1, 'Q2.5'],
    ci_97.5                             = ps[1, 'Q97.5'],
    p_6OHDA_less_than_0                 = p_pos,
    p_control_less_than_6OHDA           = p1,
    p_classic_2tail                     = p2,
    p_6OHDA_cell_less_than_control_cell = p5
  )
  rownames(bayes_summary) <- NULL

  setwd(analysis_path)
  save.image(file = RData_path)
}

# save all to single 'xlsx'
if (plotsave) {
  data_list <- setNames(
    list(df, fe, boot_summary, bayes_summary),
    c(
      paste0(cell_type, ' count data'),
      paste0(cell_type, ' lme fixed effect'),
      paste0(cell_type, ' robust lme bootstrap'),
      paste0(cell_type, ' bayes summary')
    )
  )
  list2excel(data_list, paste0(identifier, ' ', cell_type, ' RNAscope.xlsx'), wd = xlsx_path)
  list2csv(data_list,   paste0(identifier, ' ', cell_type, ' RNAscope.csv'),  wd = xlsx_path)
}

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
xrange <- c(0.75, 4.25)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 10
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5

BoxPlot(formula=var ~ Group + (1 | field), data=df[,c('field', 'Group', 'var')],
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length,
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width)

if (plotsave) save_graph(svg_path=svg_path, filename=paste0(identifier, '_RNAscope_boxplot_', cell_type, '.svg'), width=width, height=height, bg='transparent')

# Bayes

d.red <- aggregate(
  var ~ field + Group,
  data = df,
  FUN  = mean
)

names(d.red)[names(d.red) == 'var'] <- 'y.mean'

slice_ids   <- seq_along(sort(unique(df$field)))
ctrl_slices <- seq_along(sort(unique(df$field[df$Group=='Control'])))
OHDA_slices <- max(ctrl_slices) + seq_along(sort(unique(df$field[df$Group=='6OHDA'])))
divider     <- max(ctrl_slices) + 0.5
x_ctrl_txt  <- mean(range(ctrl_slices))
x_OHDA_txt   <- mean(range(OHDA_slices))
ylim        <- yrange
p.cex       <- 1.2
type        <- 16
wid         <- 0.3
lwd         <- 4/3

dev.new(width=10, height=3, noRStudioGD=TRUE)
par(las=1, mfrow=c(1, 3), mar=c(4.5, 4, 3.5, 1), cex=p.cex, tcl=-tick_length)
layout(matrix(c(1, 2, 3), nrow=1), width=c(0.5, 0.25, 0.25), height=1)

# panel A
stripchart(var~field, data=df, vertical=TRUE,
           pch=type, cex=p.cex, col='slateblue',
           ylim=ylim, xlab='SliceID', ylab='',
           frame.plot=FALSE, axes=FALSE)
axis(1, at=slice_ids, labels=slice_ids)
axis(2)
abline(v=divider, lty=2)
par(xpd=TRUE)
text(c(x_ctrl_txt, x_OHDA_txt), c(ylim[2], ylim[2]),
     labels=c('Control', '6OHDA'))
par(xpd=FALSE)
segments(slice_ids - wid, d.red$y.mean,
         slice_ids + wid, d.red$y.mean,
         lwd = 2 * lwd,
         col = 'indianred')
mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel B
par(bty='n')
lineplot.CI(Group, var, data=df, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=line_ylim, xlim=c(0.5, 2.5),
            xlab='', ylab=ylab,
            lwd=lwd, main='Pseudoreplicated')
axis(1, at=1:2, labels=levels(df$Group))
axis(2, at=seq(line_ylim[1], line_ylim[2]))
mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel C
par(bty='n')
lineplot.CI(Group, y.mean, data=d.red, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=line_ylim, xlim=c(0.5, 2.5),
            xlab='', ylab=ylab,
            lwd=lwd, main='Slice average')
axis(1, at=1:2, labels=levels(d.red$Group))
axis(2, at=seq(line_ylim[1], line_ylim[2]))
mtext('C', side=3, line=1.5, adj=0, font=2, cex=1.5)

# extract variances
vc <- as.data.frame(VarCorr(model))
σ2_slice    <- vc$vcov[vc$grp == 'field:Animal']
σ2_residual <- vc$vcov[vc$grp=='Residual']

ICC_slice <- σ2_slice / (σ2_slice + σ2_residual)
ICC_slice

dev.new(width=9, height=4.5, noRStudioGD=TRUE)

par(las=1,
    mfrow=c(1, 2),
    mar=c(4.5, 4, 3, 1),
    cex=1,
    lwd=1,
    xaxs='i',
    yaxs='i',
    tcl=-0.2)

# A. Group difference
dens <- density(post_diff, adjust=1.25)
plot(dens, main='Group difference',
     xlab='Difference (6OHDA - Control)',
     xlim=c(-5, 5), ylim=dens_ylim,
     lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
     bty='n', axes=FALSE)
axis(1, lwd=1, cex.axis=0.85)
axis(2, lwd=1, cex.axis=0.85)
box(bty='n')

polygon(c(-0.01, dens$x, 1.01),
        c(0, dens$y, 0),
        col=rgb(106/255, 90/255, 205/255, 0.6),
        border='slateblue', lwd=1)

abline(v=0, lty=2, lwd=1)
mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

legend('topleft',
       legend=c(paste0('P(diff < 0) = ', round(p1, 5)),
                paste0('two-sided p = ', round(p2, 5))),
       bty='n', box.lwd=0, cex=0.85)

# B. Posterior predictive
dens1 <- density(preds[, 1], adjust=1.25)
dens2 <- density(preds[, 2], adjust=1.25)

plot(dens1, xlab='Predicted signal size for future slice',
     main='Posterior predictive',
     ylim=pred_ylim, xlim=c(-10, 50),
     lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
     bty='n', axes=FALSE)
axis(1, lwd=1, cex.axis=0.85)
axis(2, lwd=1, cex.axis=0.85)
box(bty='n')

polygon(c(-0.01, dens1$x, 1.01),
        c(0, dens1$y, 0),
        col=rgb(106/255, 90/255, 205/255, 0.6),
        border='slateblue', lwd=1)

polygon(c(-0.01, dens2$x, 1.01),
        c(0, dens2$y, 0),
        col=rgb(205/255, 92/255, 92/255, 0.6),
        border='indianred', lwd=1)

mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

legend('topleft', col=c('slateblue', 'indianred'),
       fill=c(rgb(106/255, 90/255, 205/255, 0.6),
              rgb(205/255, 92/255, 92/255, 0.6)),
       legend=c('Control', '6-OHDA'),
       bty='n', box.lwd=0, cex=0.85)

legend('topright',
       legend=paste0('P(6-OHDA < Control) = ', round(p5, 5)),
       bty='n', box.lwd=0, cex=0.85)

if (plotsave) save_graph(svg_path=svg_path, filename=paste0(identifier, '_Bayesian_Analysis_', cell_type, '.svg'), width=width, height=height, bg='transparent')



# NDNF

# boot_summary
#                                               model    data parameter estimate     p_value    ci_2.5 ci_97.5 n_boot
# 1 rlmer bootstrap: var ~ Group + (1 | Animal/field) df_NDNF    Group1 1.454205 0.004002001 0.4468117  2.7196   1999

# bayes_summary
#                                                       model    data parameter estimate     ci_2.5  ci_97.5 p_6OHDA_less_than_0 p_control_less_than_6OHDA p_classic_2tail p_6OHDA_cell_less_than_control_cell
# 1 brm: var ~ Group + (1 | Animal/field), family = student() df_NDNF  b_Group1 1.604936 0.04065959 3.142714            0.022375                  0.977625         0.04475                             0.62025

# TH

# boot_summary
#                                               model  data parameter estimate p_value    ci_2.5  ci_97.5 n_boot
# 1 rlmer bootstrap: var ~ Group + (1 | Animal/field) df_TH    Group1 1.041789       0 0.6853093 2.281143   1999

# bayes_summary
#                                                       model  data parameter estimate    ci_2.5  ci_97.5 p_6OHDA_less_than_0 p_control_less_than_6OHDA p_classic_2tail p_6OHDA_cell_less_than_control_cell
# 1 brm: var ~ Group + (1 | Animal/field), family = student() df_TH  b_Group1 1.036545 0.2781221 1.802369             0.00425                   0.99575          0.0085                              0.6765

# fe
#                                    model  data   parameter estimate        se       df   t_value      p_value
# 1 lmer: var ~ Group + (1 | Animal/field) df_TH (Intercept) 7.070161 0.3551743 36.49453 19.906172 3.453134e-21
# 2 lmer: var ~ Group + (1 | Animal/field) df_TH      Group1 1.017490 0.3551743 36.49453  2.864761 6.881571e-03



