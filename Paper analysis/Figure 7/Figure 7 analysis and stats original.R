# ==============================================
# CREATE GRAPHS AND PERFORM STATISTICAL TESTS
# Processed data in stored in '.RDATA' form
# '.RDATA' created by '~ analysis.R'
# ==============================================

rm(list = ls(all = TRUE))
graphics.off()

plotsave <- FALSE

# Load required packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('bayesplot', 'brms', 'dbscan', 'lme4', 'minpack.lm', 'openxlsx', 'Rcpp', 'robustbase', 'robustlmm', 'sciplot', 'signal')
load_required_packages(required.packages)

# Load user config
config_path <- file.path(Sys.getenv("HOME"), ".abf2nwb_config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found. Create ~/.abf2nwb_config.yaml with your settings.")
}
config <- yaml::read_yaml(config_path)

# Construct paths
username <- config$username
file_path1 <- paste0('/Users/', username, config$path_repository)
file_path2 <- paste0('/Users/', username, config$path_analysis)

source(paste0(file_path1, '/nNLS functions.R'))

# settings
identifier <- 'Figure 7'
analysis_path <- paste0(file_path2, '/', identifier)
xlsx_path <- paste0(analysis_path, '/xlsx')
# path where all graphs are stored
svg_path <- paste0(analysis_path, '/svg')
if (!dir.exists(svg_path)) {
  dir.create(svg_path, recursive = TRUE)
}

# import ACh_GRAB.xlsx from 'xlsx' folder

setwd(analysis_path)
name <- 'ACh_GRAB'

d <- load_data2(wd=paste0(analysis_path, '/xlsx'), name=name)[[1]]
names(d)[names(d) == "(F1-F0)/F0"] <- "dff"

single_examples <- load_data2(wd=paste0(analysis_path, '/xlsx'), name='single_examples')[[1]]

d$SliceID <- as.numeric(factor(with(d, paste(Animal, Slice, sep = '_'))))
d$Condition <- factor(d$Condition, levels = c('Control', 'MCI-Park'))
options(contrasts=c('contr.sum', 'contr.poly'))

# Fit mixed-effects model: ROI nested in Slice, nested in Animal
model <- lmerTest::lmer(dff ~ Condition + (1 | SliceID), data = d)

# Summary with p-values
summary(model)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: dff ~ Condition + (1 | SliceID)
#    Data: d

# REML criterion at convergence: -14.7

# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.75185 -0.50399 -0.04547  0.48874  2.12991 

# Random effects:
#  Groups   Name        Variance Std.Dev.
#  SliceID  (Intercept) 0.03710  0.1926  
#  Residual             0.02469  0.1571  
# Number of obs: 56, groups:  SliceID, 15

# Fixed effects:
#             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)  0.60908    0.05455 13.11035  11.166 4.54e-08 ***
# Condition1  -0.11900    0.05455 13.11035  -2.182   0.0479 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Correlation of Fixed Effects:
#            (Intr)
# Condition1 -0.072

# Bootstrap

RData_path <- file.path(analysis_path, paste0('/', name, '.RData'))

# Check if the file exists
if (file.exists(RData_path)) {
  # if file exists then load it
  load(RData_path)

} else {

  # STATISTICS robust model
  rmod <- rlmer(dff ~ Condition + (1|SliceID), data=d)
  library(robustlmm)

  # Bootstrap
  # Each bootstrap replicate samples entire slices with replacement

  set.seed(42)
  B <- 9999   # #bootstrap replicates
  slice_list <- split(d, d$SliceID)
  n_slices   <- length(slice_list)
  est        <- numeric(B)

  for (ii in seq_len(B)) {
    # resampling slices with replacement
    picked <- sample(slice_list, size = n_slices, replace = TRUE)
    dboot  <- do.call(rbind, picked)
    
    # refit the robust mixed‐effects model
    mboot  <- rlmer(dff ~ Condition + (1|SliceID), data = dboot)
    
    # save treatment‐effect estimate
    est[ii] <- fixef(mboot)['Condition1']
  }

  # percentile 95% CI
  ci <- quantile(est, probs = c(0.025, 0.975))
  names(ci) <- c('2.5%', '97.5%')
  ci

#        2.5%       97.5% 
# -0.20356472 -0.02276911 

  # two-tailed bootstrap p-value
  p_lt0 <- mean(est < 0)
  pval  <- 2 * min(p_lt0, 1 - p_lt0)

  pval
  # [1] 0.01420142
  
  
  # Preparatiom for Bayesian analysis
  #  brms will automatically run chains in parallel across all available cores
  options(mc.cores = max(1, parallel::detectCores() - 1))

  ## Proposed Bayesian analysis of all data
  ## examine default priors
  # get_prior
  # Get information on all parameters (and parameter classes) for which priors may
  # be specified including default priors.
  get_prior(bf(dff ~ Condition + (1 | SliceID), center = TRUE), data=d)


 #  get_prior(bf(dff ~ Condition + (1 | SliceID), center = TRUE), data=d)
 #                  prior     class       coef   group resp dpar nlpar lb ub tag       source
 #                 (flat)         b                                                   default
 #                 (flat)         b Condition1                                   (vectorized)
 # student_t(3, 0.6, 2.5) Intercept                                                   default
 #   student_t(3, 0, 2.5)        sd                                     0             default
 #   student_t(3, 0, 2.5)        sd            SliceID                  0        (vectorized)
 #   student_t(3, 0, 2.5)        sd  Intercept SliceID                  0        (vectorized)
 #   student_t(3, 0, 2.5)     sigma                                     0             default


  # put prior on parameter for treatment effect
  # The brms() function does not calculate a prior for the difference in
  # group means (); we therefore use a Normal(0, 20) prior, which is wide 
  # relative to the size of the expected effect and therefore has little 
  # influence on the results
  m1.prior <- c(
      prior(normal(0, 20), class='b', coef='Condition1')
  )

  # fit model
  mod <- dff ~ Condition + (1 | SliceID)

  rstan:::rstan_options(disable_march_warning = TRUE) 

  # m1 <- brm(mod, data=d,
  #           iter = 1e6, 
  #           chains = 3,
  #           seed=42,
  #           prior=m1.prior, 
  #           control=list(adapt_delta=0.99)
  #       )

  m1 <- brm(
    bf(dff ~ Condition + (1|SliceID)),
    data = d,
    family = student(), 
    prior = c(
      prior(normal(0, 20), class = 'b', coef = 'Condition1'),
      prior(student_t(3, 0, 2.5), class = 'sigma')
    ),
    iter = 1e6,   
    chains = 4,
    seed = 42,
    control = list(adapt_delta = 0.99)
  )

  posterior_summary(m1, variable = 'b_Condition1')
  #                Estimate Est.Error       Q2.5       Q97.5
  # b_Condition1 -0.1167961 0.0617521 -0.2408016 0.005042222
  # e.g. Estimate, Est.Error, 95% CI
  
  # Then P(Condition1 < 0) or P(Condition1 > 0)
  draws <- as_draws_df(m1)
  p_pos <- mean(draws$b_Condition1 > 0)
  p_pos
  # [1] 0.029468

  ## posterior predictive checks
  # pp_check: Perform posterior predictive checks with the help of the bayesplotpackage

  # pp_check(m1, ndraws = 100)
  # pp_check(m1, type='stat', binwidth=0.002)
  # pp_check(m1, type='stat', stat='sd', binwidth=0.001)
  # pp_check(m1, type='stat', stat='max', binwidth=0.01)
  # pp_check(m1, type='stat', stat='min', binwidth=0.01)
  # pp_check(m1, type='intervals')

  ## extract parameter for the difference in group means
  # post_diff <- posterior_samples(m1)$b_fa1

  # post_diff <- as_draws(m1)[[chain = 3]]$b_ConditionMCIMPark
  # nb. the standard way to calculate diff is control - MCIPark

  post_diff <- -as_draws_df(m1)$b_Condition1

  ## probability that effect is less than zero
  # The interpretation is that there is a 90% chance that the control group has a lower mean 
  # than the MCI-Park group.

  p1 <- mean(post_diff > 0)
  p1
  # [1] 0.970532

  ## calculate a value similar to a classic p-value
  # calculate the probability in the other tail of the distribution (above zero),
  # and by multiplying this value by two we obtain a value that is often a similar magnitude to a
  # classic p-value. This value is 0.20, and the p-value from the classical multilevel model is 0.18.

  p2 <- 2 * mean(post_diff < 0)
  p2
  # [1] 0.058936

  ## new data used to make predictions
  new.data <- data.frame(SliceID=c(1, 1), Condition=c('Control', 'MCI-Park'))

  ## make predictions for new cells from new animals
  # If we're interested in making an inference about neurons, we can turn to the predictive
  # perspective and ask: what is the probability that a randomly chosen neuron from a future
  # slice in the treated group will have a lower value that a randomly chosen neuron from a
  # future animal in the control group?   If the scientific interest is in the individual neurons, then
  # making a probabilistic prediction about as yet unseen neurons directly addresses this question.

  preds <- posterior_predict(m1, newdata = new.data,  re_formula =  ~ (1 | SliceID), allow_new_levels=TRUE)

  ## how many signals fall within prediction interval for new slices in each group
  p3 <- mean(d$'dff'[d$Condition == 'Control'] > quantile(preds[, 1], 0.025) &
             d$'dff'[d$Condition == 'Control'] < quantile(preds[, 1], 0.975))

  p4 <- mean(d$'dff'[d$Condition == 'MCI-Park'] > quantile(preds[, 1], 0.025) &
             d$'dff'[d$Condition == 'MCI-Park'] < quantile(preds[, 1], 0.975))

  (p3 + p4)/2
  # [1]0.8474359

  p5 <- mean(preds[, 2] > preds[, 1])
  p5
  # [1] 0.817608

  setwd(analysis_path)
  save.image(file = RData_path)

}

plotsave <- FALSE

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression(dff)
xrange <- c(0.75, 4.25)
yrange <- c(0, 1.5)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 0.25
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5


BoxPlot(formula=dff ~ Condition + (1 | SliceID), data=d[,c('SliceID', 'Condition', 'dff')],
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width)
if (plotsave) save_graph(svg_path=svg_path, filename='amplitude_SPN.svg', width=width, height=height, bg='transparent')



# Bayes

# compute mean(dff) for each SliceID × Condition
d.red <- aggregate(
  dff ~ SliceID + Condition,
  data = d,
  FUN  = mean
)

# rename the aggregated column to y.mean
names(d.red)[names(d.red) == 'dff'] <- 'y.mean'

slice_ids   <- sort(unique(d$SliceID))
ctrl_slices <- sort(unique(d$SliceID[d$Condition=='Control']))
mci_slices  <- sort(unique(d$SliceID[d$Condition=='MCI-Park']))
divider     <- max(ctrl_slices) + 0.5
x_ctrl_txt  <- mean(range(ctrl_slices))
x_mci_txt   <- mean(range(mci_slices))
ylim        <- c(0, 1.4)
p.cex       <- 1.2
type        <- 16       # filled circle pch
wid         <- 0.3
lwd         <- 4/3      # original line width

dev.new(width=10, height=3, noRStudioGD=TRUE)
par(las=1, mfrow=c(1, 3), mar=c(4.5, 4, 3.5, 1), cex=p.cex, tcl=-tick_length)
layout(matrix(c(1, 2, 3), nrow=1), width=c(0.5, 0.25, 0.25), height=1)

# panel A
stripchart(dff~SliceID, data=d, vertical=TRUE,
           pch=type, cex=p.cex, col='slateblue',
           ylim=ylim, xlab='SliceID', ylab='dff',
           frame.plot=FALSE, axes=FALSE)
axis(1, at=slice_ids, labels=slice_ids)
axis(2)
abline(v=divider, lty=2)
par(xpd=TRUE)
text(c(x_ctrl_txt, x_mci_txt), c(ylim[2], ylim[2]),
     labels=c('Control', 'MCI-Park'))
par(xpd=FALSE)
segments(slice_ids - wid, d.red$y.mean,
         slice_ids + wid, d.red$y.mean,
         lwd = 2 * lwd,        # double thickness
         col = 'indianred')
mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel B
par(bty='n')
lineplot.CI(Condition, dff, data=d, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=ylim, xlim=c(0.5, 2.5),
            xlab='', ylab='dff',
            lwd=lwd, main='Pseudoreplicated')
axis(1, at=1:2, labels=levels(d$Condition))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel C
par(bty='n')
lineplot.CI(Condition, y.mean, data=d.red, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=ylim, xlim=c(0.5, 2.5),
            xlab='', ylab='dff',
            lwd=lwd, main='Slice average')
axis(1, at=1:2, labels=levels(d.red$Condition))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('C', side=3, line=1.5, adj=0, font=2, cex=1.5)

if (plotsave) save_graph(svg_path=svg_path, filename='pseudoreplication.svg', width=width, height=height, bg='transparent')

# extract variances
vc <- as.data.frame(VarCorr(model))
σ2_slice    <- vc$vcov[vc$grp=='SliceID']
σ2_residual <- vc$vcov[vc$grp=='Residual']
# ICC
ICC_slice <- σ2_slice / (σ2_slice + σ2_residual)
ICC_slice
# [1] 0.6321741

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
     xlab='Difference (MCI-Park - Control)',
     xlim=c(-0.2, 0.4), ylim=c(0, 8),
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
       legend=c(paste0('P(diff > 0) = ', round(p1, 5)),
                paste0('2 x P(diff < 0) = ', round(p2, 5))),
       bty='n', box.lwd=0, cex=0.85)

# B. Posterior predictive
dens1 <- density(preds[, 1], adjust=1.25)
dens2 <- density(preds[, 2], adjust=1.25)

plot(dens1, xlab='Predicted signal size for future slice',
     main='Posterior predictive',
     ylim=c(0, 3.5), xlim=c(0, 2),
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
       legend=c('Control', 'MCI-Park'),
       bty='n', box.lwd=0, cex=0.85)

legend('topright',
       legend=paste0('P(MCI-Park > Control) = ', round(p5, 5)),
       bty='n', box.lwd=0, cex=0.85)
if (plotsave) save_graph(svg_path=svg_path, filename='Bayesian_Analysis.svg', width=width, height=height, bg='transparent')


# dev.off()

# There is strong posterior evidence (97.2%) that MCI-Park increases the signal relative to control.
# The predictive probability that a future MCI-Park value exceeds a control value is ~83%, consistent with the relevance of the effect.
# The model fits the data well for both groups, with 73–89% of points within their own 95% prediction intervals.



# Bootstrap p value vs the Bayesian approach
# Answers two fundamentally different questions:
#   1.  Cluster‐bootstrap p-value (≈0.014)
#   • This is a frequentist test of the null hypothesis 'no difference' (Condition1 = 0).
#   • Resampled entire slices, refit a robust mixed model N boot times, and found only about 1.4% of those bootstrap estimates were as extreme—or more extreme—than zero.
#   • A small p-value like 0.014 rejects the null at the 5% level.
#   2.  Bayesian posterior tail probability (≈0.055)
#   • This is the probability that the true effect is in the opposite direction, given your prior + data.
#   • Computed 2 * mean(post_diff < 0) which tells you there’s about a 5.5% chance the effect is ≤0, after seeing the data and your (weakly informative) Normal(0,20) prior.
#   • Directly interpretable as '94.5% sure the effect is positive'.

# ⸻

# appropriate?
#   • If goal is a classical hypothesis test ('Is there evidence beyond sampling variability to reject the NULL hypothesis?'), the bootstrap p-value is valid, especially since cluster structure was modelled and a robust model was employed.
#   • If goal is a probabilistic statement about the magnitude and direction of the effect e.g. 'what is the probability the treatment truly increases signal?'' — then the Bayesian posterior (0.055) is more appropriate.

# In summary
#   • Want a yes/no decision at α=0.05? Use p-value = 0.014.
#   • Want the probability the effect is positive, accounting for prior uncertainty? Use P(effect > 0) = 0.945 (i.e. 1 – 0.055).



