# ==============================================
# ANALYSE DATA FROM XLSX FILES
# Performs fits to averaged traces
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

source('/Users/euo9382/Documents/Repositories/analysis_Belal2026/R functions/setup.R')
load_required_packages(c('bayesplot', 'brms', 'lme4', 'parallel', 'robustlmm', 'sciplot'))

identifier <- 'Figure 11'
paths <- make_paths(identifier)
analysis_path <- paths$analysis_path
xlsx_path <- paths$xlsx_path
svg_path <- paths$svg_path

# Load data
reanalysis = FALSE

if (reanalysis){
  RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_reanalyzed_all')
}else{
  RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope')
}

RNAscope_data[1:10,]
#    condition cell_type slice_id              field hemisphere field_index replicate count session
# 1     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.01         UL           1         1    25  L1.ST8
# 2     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.01         UL           1         2    28  L1.ST8
# 3     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.01         UL           1         3    13  L1.ST8
# 4     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.01         UL           1         4     7  L1.ST8
# 5     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.01         UL           1         5    12  L1.ST8
# 6     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.02         UL           2         1    16  L1.ST8
# 7     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.02         UL           2         2    10  L1.ST8
# 8     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.03         UL           3         1    37  L1.ST8
# 9     Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.03         UL           3         2    25  L1.ST8
# 10    Intact     NDNF+   L1.ST8 L1.ST8_C.UL_60x.03         UL           3         3    21  L1.ST8


RNAscope_data  <- RNAscope_data[order(RNAscope_data$condition,
                                      RNAscope_data$cell_type,
                                      RNAscope_data$hemisphere,
                                      RNAscope_data$slice_id), ]

rownames(RNAscope_data) = 1:dim(RNAscope_data)[1]
# 
ndnf_data <- subset(RNAscope_data, cell_type == "NDNF+")
th_data <- subset(RNAscope_data, cell_type == "TH+")

wilcox.test(
  count ~ condition,
  data = ndnf_data,
  subset = condition %in% c("Intact", "Lesioned")
)

#   Wilcoxon rank sum test with continuity correction

# data:  count by condition
# W = 4280, p-value = 0.008855
# alternative hypothesis: true location shift is not equal to 0


wilcox.test(
  count ~ condition,
  data = th_data,
  subset = condition %in% c("Intact", "Lesioned")
)

#   Wilcoxon rank sum test with continuity correction

# data:  count by condition
# W = 1267.5, p-value = 0.005944
# alternative hypothesis: true location shift is not equal to 0


df_NDNF <- subset(RNAscope_data, cell_type == "NDNF+")
df_NDNF$Animal <- sub("^(L\\d+).*", "\\1", as.character(df_NDNF$slice_id))
df_NDNF$Group <- ifelse(df_NDNF$hemisphere == "UL", "Control", "6OHDA")

df_NDNF$Animal <- factor(df_NDNF$Animal)
df_NDNF$SliceID <- factor(df_NDNF$slice_id)
df_NDNF$Group <- factor(df_NDNF$Group, levels = c("Control", "6OHDA"))
df_NDNF$var <- df_NDNF$count


df_TH <- subset(RNAscope_data, cell_type == "TH+")
df_TH$Animal <- sub("^(L\\d+).*", "\\1", as.character(df_TH$slice_id))
df_TH$Group <- ifelse(df_TH$hemisphere == "UL", "Control", "6OHDA")

df_TH$Animal <- factor(df_TH$Animal)
df_TH$SliceID <- factor(df_TH$slice_id)
df_TH$Group <- factor(df_TH$Group, levels = c("Control", "6OHDA"))
df_TH$var <- df_TH$count


# Combine into one dataframe
df_NDNF <- rbind(control_df, lesioned_df)
rownames(df_NDNF) <- NULL


df_NDNF$SliceID <- factor(df_NDNF$SliceID)
df_NDNF$Group <- factor(df_NDNF$Group, levels = c("Control", "6OHDA"))

# TH+
control_df <- build_df(control_TH, "Control")
lesioned_df <- build_df(lesioned_TH, "6OHDA", 
                        start_animal = length(control_TH) + 1,
                        start_slice = max(control_df$SliceID) + 1)

# Combine into one dataframe
df_TH <- rbind(control_df, lesioned_df)
rownames(df_TH)  <- NULL


df1$SliceID <- factor(df1$SliceID)
df1$Group <- factor(df1$Group, levels = c("Control", "6OHDA"))




options(contrasts=c('contr.sum', 'contr.poly'))

# Fit mixed-effects model: ROI nested in slice
# raw ROI-level data: field is the more important clustering unit
# best if supported: both SliceID and field
# not ideal to use only SliceID when you have repeated ROIs within field

# model_NDNF <- lmerTest::lmer(var ~ Group + (1 | SliceID), data = df_NDNF)
model_NDNF <- lmerTest::lmer(var ~ Group + (1 | field), data = df_NDNF)
# model_NDNF <- lmerTest::lmer(var ~ Group + (1 | SliceID) + (1 | field), data = df_NDNF)

# Summary with p-values
summary(model_NDNF)

# Fit mixed-effects model: ROI nested in Animal
model_TH <- lmerTest::lmer(var ~ Group + (1 | field), data = df_TH)
# model_TH <- lmerTest::lmer(var ~ Group + (1 | SliceID) + (1 | field), data = df_TH)


# Summary with p-values
summary(model_TH)


# Bootstrap
name <- 'RNAscope_TH'
RData_path <- file.path(analysis_path, paste0('/', name, '.RData'))

# Check if the file exists
if (file.exists(RData_path)) {
  # if file exists then load it
  load(RData_path)

} else {

  # STATISTICS robust model
  df <- df_TH

  df$field <- factor(df$field)
  df$Group <- factor(df$Group, levels = c("Control", "6OHDA"))

  rmod <- rlmer(var ~ Group + (1 | field), data = df)

  set.seed(42)
  B <- 9999
  field_list <- split(df, df$field)
  n <- length(field_list)

  coef_name <- grep("^Group", names(fixef(rmod)), value = TRUE)[1]

  boot_one <- function(ii) {
    picked_idx <- sample(seq_along(field_list), size = n, replace = TRUE)

    picked <- lapply(seq_along(picked_idx), function(j) {
      d <- field_list[[picked_idx[j]]]
      d$field <- paste0(as.character(d$field[1]), "_boot", ii, "_", j)
      d
    })

    dboot <- do.call(rbind, picked)
    dboot$field <- factor(dboot$field)
    dboot$Group <- factor(dboot$Group, levels = levels(df$Group))

    fit <- try(rlmer(var ~ Group + (1 | field), data = dboot), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(NA_real_)
    }

    unname(fixef(fit)[coef_name])
  }

  n_cores <- max(1, detectCores() - 1)

  est <- unlist(
    mclapply(seq_len(B), boot_one, mc.cores = n_cores),
    use.names = FALSE
  )

  est <- est[is.finite(est)]

  ci <- quantile(est, probs = c(0.025, 0.975))
  names(ci) <- c("2.5%", "97.5%")
  ci

  p_lt0 <- mean(est < 0)
  pval <- 2 * min(p_lt0, 1 - p_lt0)
  pval

  # [1] 0.00480048
  
  # Preparation for Bayesian analysis
  # brms will automatically run chains in parallel across all available cores
  options(mc.cores = max(1, parallel::detectCores() - 1))

  # Proposed Bayesian analysis of all data
  # examine default priors
  # get_prior
  # Get information on all parameters (and parameter classes) for which priors may
  # be specified including default priors.

  get_prior(bf(var ~ Group + (1 | field), center = TRUE), data = df)

  # put prior on parameter for treatment effect
  # The brms() function does not calculate a prior for the difference in
  # group means directly, so we use a Normal(0, 20) prior on the Group coefficient.
  # This is wide relative to the expected effect size and therefore only weakly informative.
  m1.prior <- c(
    prior(normal(0, 20), class = "b")
  )

  # fit model
  mod <- var ~ Group + (1 | field)

  rstan:::rstan_options(disable_march_warning = TRUE)

  m1 <- brm(
    bf(var ~ Group + (1 | field)),
    data = df,
    family = student(),
    prior = c(
      prior(normal(0, 20), class = "b"),
      prior(student_t(3, 0, 2.5), class = "sigma")
    ),
    iter = 4000,
    chains = 4,
    seed = 42,
    control = list(adapt_delta = 0.99)
  )

  # extract the actual Group coefficient name safely
  coef_name <- grep("^b_Group", names(as_draws_df(m1)), value = TRUE)[1]

  posterior_summary(m1, variable = coef_name)
  # e.g. Estimate, Est.Error, 95% CI

  # Then P(Group effect < 0) or P(Group effect > 0)
  draws <- as_draws_df(m1)[[coef_name]]
  p_pos <- mean(draws < 0)
  p_pos

  # 0.0305

  ## posterior predictive checks
  # pp_check: Perform posterior predictive checks with the help of the bayesplot package

  # pp_check(m1, ndraws = 100)
  # pp_check(m1, type = "stat", binwidth = 0.002)
  # pp_check(m1, type = "stat", stat = "sd", binwidth = 0.001)
  # pp_check(m1, type = "stat", stat = "max", binwidth = 0.01)
  # pp_check(m1, type = "stat", stat = "min", binwidth = 0.01)
  # pp_check(m1, type = "intervals")

  ## extract parameter for the difference in group means
  # nb. the standard way to calculate diff is control - 6OHDA
  post_diff <- -draws

  ## probability that effect is less than zero
  # The interpretation is that there is a high probability that the control group
  # has a lower mean than the 6-OHDA group if post_diff < 0 with high probability.
  p1 <- mean(post_diff < 0)
  p1
  # 0.995125

  ## calculate a value similar to a classic p-value
  # calculate the probability in the other tail of the distribution (above zero),
  # and by multiplying this value by two we obtain a value often similar in magnitude
  # to a classic p-value
  p2 <- 2 * mean(post_diff > 0)
  p2
  # 0.00975

  ## new data used to make predictions
  new.data <- data.frame(
    field = factor(c("new_control_field", "new_6ohda_field")),
    Group = factor(c("Control", "6OHDA"), levels = c("Control", "6OHDA"))
  )

  # make predictions for new cells from new fields
  # If interested in making an inference about neurons, then use predictive perspective 
  # and ask: what is the probability that a randomly chosen neuron from a future field 
  # in the treated group will have a lower value than a randomly chosen neuron from a
  # future field in the control group?

  preds <- posterior_predict(
    m1,
    newdata = new.data,
    re_formula = ~(1 | field),
    allow_new_levels = TRUE
  )

  ## how many signals fall within prediction interval for new fields in each group
  p3 <- mean(
    df$var[df$Group == "Control"] > quantile(preds[, 1], 0.025) &
    df$var[df$Group == "Control"] < quantile(preds[, 1], 0.975)
  )

  p4 <- mean(
    df$var[df$Group == "6OHDA"] > quantile(preds[, 2], 0.025) &
    df$var[df$Group == "6OHDA"] < quantile(preds[, 2], 0.975)
  )

  (p3 + p4) / 2
  # 0.9650794

  p5 <- mean(preds[, 2] < preds[, 1])
  p5
  # 0.6785


  setwd(analysis_path)
  save.image(file = RData_path)

}

plotsave <- FALSE

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression('CHRNB2 Particles/TH Cell')
xrange <- c(0.75, 4.25)
yrange <- c(0, 20)
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
if (plotsave) save_graph(svg_path=svg_path, filename='RNAscope_boxplot_TH.svg', width=width, height=height, bg='transparent')

# Bayes

# compute mean((F1-F0)/F0) for each field × Condition
d.red <- aggregate(
  var ~ field + Group,
  data = df,
  FUN  = mean
)

# rename the aggregated column to y.mean
names(d.red)[names(d.red) == 'var'] <- 'y.mean'

slice_ids   <- as.numeric(sort(unique(df$field)))
ctrl_slices <- as.numeric(sort(unique(df$field[df$Group=='Control'])))
OHDA_slices  <- as.numeric(sort(unique(df$field[df$Group=='6OHDA'])))
divider     <- max(ctrl_slices) + 0.5
x_ctrl_txt  <- mean(range(ctrl_slices))
x_OHDA_txt   <- mean(range(OHDA_slices))
ylim        <- c(0, 20)
p.cex       <- 1.2
type        <- 16       # filled circle pch
wid         <- 0.3
lwd         <- 4/3      # original line width

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
         lwd = 2 * lwd,        # double thickness
         col = 'indianred')
mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel B
par(bty='n')
lineplot.CI(Group, var, data=df, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=c(4,12), xlim=c(0.5, 2.5),
            xlab='', ylab='CHRNB2 Particles/TH Cell',
            lwd=lwd, main='Pseudoreplicated')
axis(1, at=1:2, labels=levels(df$Group))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel C
par(bty='n')
lineplot.CI(Group, y.mean, data=d.red, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=c(4,12), xlim=c(0.5, 2.5),
            xlab='', ylab='CHRNB2 Particles/NDNF Cell',
            lwd=lwd, main='Slice average')
axis(1, at=1:2, labels=levels(d.red$Group))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('C', side=3, line=1.5, adj=0, font=2, cex=1.5)

# if (plotsave) save_graph(svg_path=svg_path, filename='pseudoreplication_TH.svg', width=width, height=height, bg='transparent')

# extract variances
model <- model_TH

vc <- as.data.frame(VarCorr(model))
σ2_slice    <- vc$vcov[vc$grp=='field']
σ2_residual <- vc$vcov[vc$grp=='Residual']
# ICC
ICC_slice <- σ2_slice / (σ2_slice + σ2_residual)
ICC_slice
# [1] 0.4059719

dev.new(width=9, height=4.5, noRStudioGD=TRUE)

par(las=1,
    mfrow=c(1, 2),
    mar=c(4.5, 4, 3, 1),
    cex=1,
    lwd=1,
    xaxs='i',
    yaxs='i',
    tcl=-0.2)  # Tick length

# A. Group difference
dens <- density(post_diff, adjust=1.25)
plot(dens, main='Group difference',
     xlab='Difference (6OHDA - Control)',
     xlim=c(-5, 5), ylim=c(0, 1.2),
     lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
     bty='n', axes=FALSE)
axis(1, lwd=1, cex.axis=0.85)
axis(2, lwd=1, cex.axis=0.85)
box(bty='n')

# slateblue fill with alpha 0.25
polygon(c(-0.01, dens$x, 1.01),
        c(0, dens$y, 0),
        col=rgb(106/255, 90/255, 205/255, 0.6),
        border='slateblue', lwd=1)

abline(v=0, lty=2, lwd=1)
mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

legend('topleft',
       legend=c(paste0('P(diff < 0) = ', round(p1, 5)),
                paste0('2 x P(diff > 0) = ', round(p2, 5))),
       bty='n', box.lwd=0, cex=0.85)

# B. Posterior predictive
dens1 <- density(preds[, 1], adjust=1.25)
dens2 <- density(preds[, 2], adjust=1.25)

plot(dens1, xlab='Predicted signal size for future slice',
     main='Posterior predictive',
     ylim=c(0, 0.16), xlim=c(-10, 50),
     lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
     bty='n', axes=FALSE)
axis(1, lwd=1, cex.axis=0.85)
axis(2, lwd=1, cex.axis=0.85)
box(bty='n')

# slateblue fill
polygon(c(-0.01, dens1$x, 1.01),
        c(0, dens1$y, 0),
        col=rgb(106/255, 90/255, 205/255, 0.6),
        border='slateblue', lwd=1)
# indianred fill
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

if (plotsave) save_graph(svg_path=svg_path, filename='Bayesian_Analysis_TH.svg', width=width, height=height, bg='transparent')








