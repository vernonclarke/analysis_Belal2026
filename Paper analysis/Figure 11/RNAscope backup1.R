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
RNAscope_data <- load_data(wd=xlsx_path, name='RNAscope_reanalyzed_all')

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

RNAscope_data2 <- load_data(wd=xlsx_path, name='RNAscope')

RNAscope_data  <- RNAscope_data[order(RNAscope_data$condition,
                                      RNAscope_data$cell_type,
                                      RNAscope_data$session,
                                      RNAscope_data$hemisphere,
                                      RNAscope_data$field_index,
                                      RNAscope_data$replicate), ]

RNAscope_data2 <- RNAscope_data2[order(RNAscope_data2$condition,
                                       RNAscope_data2$cell_type,
                                       RNAscope_data2$slice_id,
                                       RNAscope_data2$hemisphere,
                                       RNAscope_data2$field_index,
                                       RNAscope_data2$replicate), ]

rownames(RNAscope_data) = 1:dim(RNAscope_data)[1]
rownames(RNAscope_data2) = 1:dim(RNAscope_data2)[1]



# 
ndnf_data <- subset(RNAscope_data, cell_type == "NDNF+")
th_data <- subset(RNAscope_data, cell_type == "TH+")

wilcox.test(
  count ~ condition,
  data = ndnf_data,
  subset = condition %in% c("Intact", "Lesioned")
)

wilcox.test(
  count ~ condition,
  data = th_data,
  subset = condition %in% c("Intact", "Lesioned")
)

slice_means1 <- aggregate(
  count ~ condition + cell_type + session + hemisphere + field,
  data = RNAscope_data,
  FUN = mean
)

names(slice_means1)[names(slice_means1) == "count"] <- "mean_obs1"

slice_means2 <- aggregate(
  count ~ condition + cell_type + slice_id + hemisphere + field,
  data = RNAscope_data2,
  FUN = mean
)

names(slice_means2)[names(slice_means2) == "count"] <- "mean_obs2"

slice_means1$key <- paste(slice_means1$condition,
                          slice_means1$cell_type,
                          slice_means1$field,
                          sep = "_")

slice_means2$key <- paste(slice_means2$condition,
                          slice_means2$cell_type,
                          slice_means2$field,
                          sep = "_")

merged <- merge(
  slice_means1[, c("condition", "cell_type", "session", "hemisphere", "field", "mean_obs1", "key")],
  slice_means2[, c("slice_id", "mean_obs2", "key")],
  by = "key"
)

merged$diff_mean <- merged$mean_obs1 - merged$mean_obs2
merged$abs_diff  <- abs(merged$diff_mean)

merged[order(merged$abs_diff, decreasing = TRUE),
       c("condition", "cell_type", "session", "slice_id", "hemisphere", "field", "mean_obs1", "mean_obs2", "diff_mean")][1:10, ]




make_nested_list <- function(df) {
  slice_ids <- unique(df$slice_id)

  out <- lapply(slice_ids, function(s) {
    sub <- df[df$slice_id == s, , drop = FALSE]
    sub <- sub[order(sub$field_index, sub$replicate), , drop = FALSE]

    field_ids <- seq_len(max(sub$field_index, na.rm = TRUE))

    fields <- lapply(field_ids, function(f) {
      sub$count[sub$field_index == f]
    })

    # drop empty fields
    fields[sapply(fields, length) > 0]
  })

  out
}

control_NDNF <- make_nested_list(
  RNAscope_data[RNAscope_data$hemisphere == "UL" &
                 RNAscope_data$cell_type == "NDNF+", ]
)

lesioned_NDNF <- make_nested_list(
  RNAscope_data[RNAscope_data$hemisphere == "L" &
                 RNAscope_data$cell_type == "NDNF+", ]
)

control_TH <- make_nested_list(
  RNAscope_data[RNAscope_data$hemisphere == "UL" &
                 RNAscope_data$cell_type == "TH+", ]
)

lesioned_TH <- make_nested_list(
  RNAscope_data[RNAscope_data$hemisphere == "L" &
                 RNAscope_data$cell_type == "TH+", ]
)


control_NDNF2 <- make_nested_list(
  RNAscope_data2[RNAscope_data2$hemisphere == "UL" &
                 RNAscope_data2$cell_type == "NDNF+", ]
)

lesioned_NDNF2 <- make_nested_list(
  RNAscope_data2[RNAscope_data2$hemisphere == "L" &
                 RNAscope_data2$cell_type == "NDNF+", ]
)



### Paste the lists first (from previous step) ###
# control <- list(...); lesioned <- list(...)

build_df <- function(group_list, group_name, start_animal = 1, start_slice = 1) {
  df <- data.frame()
  slice_id <- start_slice
  for (a in seq_along(group_list)) {
    animal_slices <- group_list[[a]]
    for (slice in animal_slices) {
      df <- rbind(df, data.frame(
        Animal = start_animal + a - 1,
        SliceID = slice_id,
        Group = group_name,
        var = slice
      ))
      slice_id <- slice_id + 1
    }
  }
  df$Animal <- factor(df$Animal)
  df$Group <- factor(df$Group, levels = c("Control", "6OHDA"))
  return(df)
}

# Convert both groups
control_df <- build_df(control_NDNF, "Control")
lesioned_df <- build_df(lesioned_NDNF, "6OHDA", 
                        start_animal = length(control_NDNF) + 1,
                        start_slice = max(control_df$SliceID) + 1)

# Combine into one dataframe
df1 <- rbind(control_df, lesioned_df)
rownames(df) <- NULL


df1$SliceID <- factor(df1$SliceID)
df1$Group <- factor(df1$Group, levels = c("Control", "6OHDA"))



# Convert both groups
control_df2 <- build_df(control_NDNF2, "Control")
lesioned_df2 <- build_df(lesioned_NDNF2, "6OHDA", 
                        start_animal = length(control_NDNF2) + 1,
                        start_slice = max(control_df2$SliceID) + 1)

# Combine into one dataframe
df2 <- rbind(control_df2, lesioned_df2)
rownames(df2) <- NULL


df2$SliceID <- factor(df2$SliceID)
df2$Group <- factor(df2$Group, levels = c("Control", "6OHDA"))


df1
df2


options(contrasts=c('contr.sum', 'contr.poly'))

# Fit mixed-effects model: ROI nested in Animal
model <- lmerTest::lmer(var ~ Group + (1 | SliceID), data = df1)

# Summary with p-values
summary(model)




# Fit mixed-effects model: ROI nested in Animal
model2 <- lmerTest::lmer(var ~ Group + (1 | SliceID), data = df2)

# Summary with p-values
summary(model2)


slice_means <- aggregate(var ~ Group + Animal + SliceID, data = df1, FUN = mean)
names(slice_means)[names(slice_means) == "var"] <- "mean_var"


slice_means2 <- aggregate(var ~ Group + Animal + SliceID, data = df2, FUN = mean)
names(slice_means2)[names(slice_means2) == "var"] <- "mean_var"


merged <- merge(slice_means, slice_means2,
                by = c("Group", "SliceID"),
                suffixes = c("_obs1", "_obs2"))
# correlation
cor(merged$mean_var_obs1, merged$mean_var_obs2)
# [1] 0.9335358

# difference
merged$diff_mean <- merged$mean_var_obs1 - merged$mean_var_obs2

summary(merged$diff_mean)



t.test(diff_mean ~ Group, data = merged)


# The difference between observers depends on Group

# Means:
#   • Control: -1.99
#   • 6OHDA: -0.14

# Observer 1 counts ~2 units lower than Observer 2 in Control
# but almost no difference in Lesioned

# merged$mean <- rowMeans(cbind(merged$mean_var_obs1, merged$mean_var_obs2))
# plot(merged$mean, merged$diff_mean)
# abline(h=mean(merged$diff_mean), col="red")


# model_diff <- lm(diff_mean ~ Group + mean, data = merged)
# summary(model_diff)

# Call:
# lm(formula = diff_mean ~ Group + mean, data = merged)

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -10.4529  -1.1481   0.1675   1.3053   4.8046 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  0.26842    1.09126   0.246   0.8067  
# Group1      -0.81133    0.37164  -2.183   0.0337 *
# mean        -0.06612    0.05101  -1.296   0.2007  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 2.646 on 51 degrees of freedom
# Multiple R-squared:  0.1393,  Adjusted R-squared:  0.1055 
# F-statistic: 4.126 on 2 and 51 DF,  p-value: 0.02183








# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: var ~ Group + (1 | SliceID)
#    Data: df

# REML criterion at convergence: 1148.4

# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.64970 -0.64996 -0.08153  0.56774  2.84522 

# Random effects:
#  Groups   Name        Variance Std.Dev.
#  SliceID  (Intercept) 28.10    5.301   
#  Residual             41.12    6.412   
# Number of obs: 167, groups:  SliceID, 54

# Fixed effects:
#             Estimate Std. Error     df t value Pr(>|t|)    
# (Intercept)   20.461      0.898 42.091  22.785   <2e-16 ***
# Group1         1.827      0.898 42.091   2.034   0.0483 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Correlation of Fixed Effects:
#        (Intr)
# Group1 0.078  

# 2 * pt(-abs(2.034), df = 42.091)
# [1] 0.04828563

# Bootstrap
name <- 'RNAscope'
RData_path <- file.path(analysis_path, paste0('/', name, '.RData'))

# Check if the file exists
if (file.exists(RData_path)) {
  # if file exists then load it
  load(RData_path)

} else {

  # STATISTICS robust model


  df = df1
  # df = df2

  rmod <- rlmer(var ~ Group + (1 | SliceID), data = df)

  set.seed(42)
  B <- 9999
  slice_list <- split(df, df$SliceID)
  n <- length(slice_list)

  coef_name <- grep("^Group", names(fixef(rmod)), value = TRUE)[1]

  boot_one <- function(ii) {
    picked_idx <- sample(seq_along(slice_list), size = n, replace = TRUE)

    picked <- lapply(seq_along(picked_idx), function(j) {
      d <- slice_list[[picked_idx[j]]]
      d$SliceID <- paste0(d$SliceID[1], "_boot", ii, "_", j)
      d
    })

    dboot <- do.call(rbind, picked)

    fit <- try(rlmer(var ~ Group + (1 | SliceID), data = dboot), silent = TRUE)
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

  # [1] 0.06740674

  # # compare
  # rmod <- rlmer(var ~ Group + (1 | SliceID), data = df2)

  # slice_list <- split(df2, df2$SliceID)
  # n <- length(slice_list)
  # B <- 999
  # n_cores <- max(1, detectCores() - 1)

  # coef_name <- grep("^Group", names(fixef(rmod)), value = TRUE)[1]

  # boot_old <- function(ii) {
  #   picked <- sample(slice_list, size = n, replace = TRUE)
  #   dboot <- do.call(rbind, picked)

  #   fit <- try(rlmer(var ~ Group + (1 | SliceID), data = dboot), silent = TRUE)
  #   if (inherits(fit, "try-error")) {
  #     return(NA_real_)
  #   }

  #   unname(fixef(fit)[coef_name])
  # }

  # boot_new <- function(ii) {
  #   picked_idx <- sample(seq_along(slice_list), size = n, replace = TRUE)

  #   picked <- lapply(seq_along(picked_idx), function(j) {
  #     d <- slice_list[[picked_idx[j]]]
  #     d$SliceID <- paste0(d$SliceID[1], "_boot", ii, "_", j)
  #     d
  #   })

  #   dboot <- do.call(rbind, picked)

  #   fit <- try(rlmer(var ~ Group + (1 | SliceID), data = dboot), silent = TRUE)
  #   if (inherits(fit, "try-error")) {
  #     return(NA_real_)
  #   }

  #   unname(fixef(fit)[coef_name])
  # }

  # set.seed(42)
  # est_old <- unlist(
  #   mclapply(seq_len(B), boot_old, mc.cores = n_cores),
  #   use.names = FALSE
  # )
  # est_old <- est_old[is.finite(est_old)]

  # set.seed(42)
  # est_new <- unlist(
  #   mclapply(seq_len(B), boot_new, mc.cores = n_cores),
  #   use.names = FALSE
  # )
  # est_new <- est_new[is.finite(est_new)]

  # ci_old <- quantile(est_old, probs = c(0.025, 0.975))
  # ci_new <- quantile(est_new, probs = c(0.025, 0.975))

  # p_old <- 2 * min(mean(est_old < 0), 1 - mean(est_old < 0))
  # p_new <- 2 * min(mean(est_new < 0), 1 - mean(est_new < 0))

  # list(
  #   n_old = length(est_old),
  #   n_new = length(est_new),
  #   ci_old = ci_old,
  #   ci_new = ci_new,
  #   p_old = p_old,
  #   p_new = p_new
  # )


  # [1] 0.008008008
  
  # Preparatiom for Bayesian analysis
  #  brms will automatically run chains in parallel across all available cores
  options(mc.cores = max(1, parallel::detectCores() - 1))

  ## Proposed Bayesian analysis of all data
  ## examine default priors
  # get_prior
  # Get information on all parameters (and parameter classes) for which priors may
  # be specified including default priors.

  get_prior(bf(var ~ Group + (1 | SliceID), center = TRUE), data=df)


   #                 prior     class      coef   group resp dpar nlpar lb ub       source
   #                (flat)         b                                              default
   #                (flat)         b    Group1                               (vectorized)
   # student_t(3, 19, 7.4) Intercept                                              default
   #  student_t(3, 0, 7.4)        sd                                    0         default
   #  student_t(3, 0, 7.4)        sd           SliceID                  0    (vectorized)
   #  student_t(3, 0, 7.4)        sd Intercept SliceID                  0    (vectorized)
   #  student_t(3, 0, 7.4)     sigma                                    0         default



  # put prior on parameter for treatment effect
  # The brms() function does not calculate a prior for the difference in
  # group means (); we therefore use a Normal(0, 20) prior, which is wide 
  # relative to the size of the expected effect and therefore has little 
  # influence on the results
  m1.prior <- c(
      prior(normal(0, 20), class='b', coef='Group1')
  )

  # fit model
  mod <- var ~ Group + (1 | SliceID)

  rstan:::rstan_options(disable_march_warning = TRUE) 

  # m1 <- brm(mod, data=d,
  #           iter = 1e6, 
  #           chains = 3,
  #           seed=42,
  #           prior=m1.prior, 
  #           control=list(adapt_delta=0.99)
  #       )

  m1 <- brm(
    bf(var ~ Group + (1|SliceID)),
    data = df,
    family = student(), 
    prior = c(
      prior(normal(0, 20), class = 'b', coef = 'Group1'),
      prior(student_t(3, 0, 2.5), class = 'sigma')
    ),
    iter = 1e6,   
    chains = 4,
    seed = 42,
    control = list(adapt_delta = 0.99)
  )

  posterior_summary(m1, variable = 'b_Group1')
  #                  Estimate Est.Error       Q2.5       Q97.5
  #   b_Condition1 -0.1224671  0.064238 -0.2508308 0.004945209
  # e.g. Estimate, Est.Error, 95% CI
  
  # Then P(Condition1 < 0) or P(Condition1 > 0)
  p_pos <- mean(as_draws_df(m1)$b_Group1 < 0)
  p_pos
  # [1] 0.031584

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

  post_diff <- -as_draws_df(m1)$b_Group1

  ## probability that effect is less than zero
  # The interpretation is that there is a 90% chance that the control group has a lower mean 
  # than the 6-OHDA group.

  p1 <- mean(post_diff < 0)
  p1
  # [1]  0.968416

  ## calculate a value similar to a classic p-value
  # calculate the probability in the other tail of the distribution (above zero),
  # and by multiplying this value by two we obtain a value that is often a similar magnitude to a
  # classic p-value. This value is 0.20, and the p-value from the classical multilevel model is 0.18.

  p2 <- 2 * mean(post_diff > 0)
  p2
  # [1] 0.063168

  ## new data used to make predictions
  new.data <- data.frame(SliceID=c(1, 1), Group=c('Control', '6OHDA'))

  ## make predictions for new cells from new animals
  # If we're interested in making an inference about neurons, we can turn to the predictive
  # perspective and ask: what is the probability that a randomly chosen neuron from a future
  # slice in the treated group will have a lower value that a randomly chosen neuron from a
  # future animal in the control group?   If the scientific interest is in the individual neurons, then
  # making a probabilistic prediction about as yet unseen neurons directly addresses this question.

  preds <- posterior_predict(m1, newdata = new.data,  re_formula =  ~ (1 | SliceID), allow_new_levels=TRUE)

  ## how many signals fall within prediction interval for new slices in each group
  p3 <- mean(df$var[df$Group == 'Control'] > quantile(preds[, 1], 0.025) &
             df$var[df$Group == 'Control'] < quantile(preds[, 1], 0.975))

  p4 <- mean(df$var[df$Group == '6OHDA'] > quantile(preds[, 1], 0.025) &
             df$var[df$Group == '6OHDA'] < quantile(preds[, 1], 0.975))

  (p3 + p4)/2
  # [1] 0.9026696

  p5 <- mean(preds[, 2] < preds[, 1])
  p5
  # 0.6421565

  setwd(analysis_path)
  save.image(file = RData_path)

}


plotsave <- FALSE

# graph settings
wid <- 0.3
cap <- 0.05
xlab <- ''
ylab <- expression('CHRNB2 Particles/NDNF Cell')
xrange <- c(0.75, 4.25)
yrange <- c(0, 60)
lwd <- 4/3
type <- 6
tick_length <- 0.2
y_tick_interval <- 10
amount <- 0.05
p.cex <- 0.6
width <- 3
height <- 3.5


BoxPlot(formula=var ~ Group + (1 | SliceID), data=df[,c('SliceID', 'Group', 'var')],
  wid=wid, cap=cap, xlab='', ylab = ylab, xrange=xrange, yrange=yrange,  xlabel_angle=45, tick_length=tick_length, 
  y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width)
if (plotsave) save_graph(svg_path=svg_path, filename='RNAscope_boxplot.svg', width=width, height=height, bg='transparent')

# Bayes

# compute mean((F1-F0)/F0) for each SliceID × Condition
d.red <- aggregate(
  var ~ SliceID + Group,
  data = df,
  FUN  = mean
)

# rename the aggregated column to y.mean
names(d.red)[names(d.red) == 'var'] <- 'y.mean'

slice_ids   <- as.numeric(sort(unique(df$SliceID)))
ctrl_slices <- as.numeric(sort(unique(df$SliceID[df$Group=='Control'])))
OHDA_slices  <- as.numeric(sort(unique(df$SliceID[df$Group=='6OHDA'])))
divider     <- max(ctrl_slices) + 0.5
x_ctrl_txt  <- mean(range(ctrl_slices))
x_OHDA_txt   <- mean(range(OHDA_slices))
ylim        <- c(0, 60)
p.cex       <- 1.2
type        <- 16       # filled circle pch
wid         <- 0.3
lwd         <- 4/3      # original line width

dev.new(width=10, height=3, noRStudioGD=TRUE)
par(las=1, mfrow=c(1, 3), mar=c(4.5, 4, 3.5, 1), cex=p.cex, tcl=-tick_length)
layout(matrix(c(1, 2, 3), nrow=1), width=c(0.5, 0.25, 0.25), height=1)

# panel A
stripchart(var~SliceID, data=df, vertical=TRUE,
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
            ylim=c(16,26), xlim=c(0.5, 2.5),
            xlab='', ylab='CHRNB2 Particles/NDNF Cell',
            lwd=lwd, main='Pseudoreplicated')
axis(1, at=1:2, labels=levels(df$Group))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

# panel C
par(bty='n')
lineplot.CI(Group, y.mean, data=d.red, type='p',
            pch=type, cex=p.cex, col='slateblue',
            ylim=c(16,26), xlim=c(0.5, 2.5),
            xlab='', ylab='CHRNB2 Particles/NDNF Cell',
            lwd=lwd, main='Slice average')
axis(1, at=1:2, labels=levels(d.red$Group))
axis(2, at=seq(ylim[1], ylim[2]))
mtext('C', side=3, line=1.5, adj=0, font=2, cex=1.5)

# if (plotsave) save_graph(svg_path=svg_path, filename='pseudoreplication.svg', width=width, height=height, bg='transparent')

# extract variances
vc <- as.data.frame(VarCorr(model))
σ2_slice    <- vc$vcov[vc$grp=='SliceID']
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
     xlim=c(-5, 5), ylim=c(0, 0.6),
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
     ylim=c(0, 0.08), xlim=c(-10, 50),
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
if (plotsave) save_graph(svg_path=svg_path, filename='Bayesian_Analysis.svg', width=width, height=height, bg='transparent')








