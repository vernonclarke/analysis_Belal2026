# nNLS functions
# Rewrite residFun in C++ and SS.fun in C++


# # Define the C++ code as a string
# cpp_code <- '
# #include <Rcpp.h>
# using namespace Rcpp;

# // Function to calculate residuals
# // [[Rcpp::export]]
# NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI) {
#   int n=y.size();
#   NumericVector residuals(n);
  
#   NumericVector fitted_values=as<NumericVector>(func(params, x, N, IEI));
  
#   for (int i=0; i < n; i++) {
#     residuals[i]=y[i] - fitted_values[i];
#   }
  
#   return residuals;
# }

# // Function to calculate sum of squares
# // [[Rcpp::export]]
# double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI) {
#   NumericVector residuals=residFunCpp(params, y, x, func, N, IEI);
#   return sum(pow(residuals, 2));
# }
# '

# # Define the C++ code as a string
# cpp_code <- '
# #include <Rcpp.h>
# using namespace Rcpp;

# // Function to calculate residuals
# // [[Rcpp::export]]
# NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI) {
#   NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));
#   return y - fitted_values;
# }
# // Function to calculate sum of squares
# // [[Rcpp::export]]
# double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI) {
#   NumericVector residuals=residFunCpp(params, y, x, func, N, IEI);
#   return sum(pow(residuals, 2));
# }
# '

# Define the C++ code as a string
cpp_code <- '
#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate residuals with weight_method
// [[Rcpp::export]]
NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));
  NumericVector residuals = y - fitted_values;
  
  // Apply weights based on the method
  if (weight_method == "~y_sqrt") {
    // Use the absolute value of y to avoid NaN issues
    NumericVector y_safe = abs(y);  // Ensure no negative values in y
    residuals = sqrt(y_safe) * residuals;  // Apply sqrt of y as weights
  } else if (weight_method == "~y") {
    residuals = y * residuals;  // 
  } else if (weights.isNotNull()) {
    NumericVector w = weights.get();
    residuals = sqrt(w) * residuals;  // Apply custom weights if provided
  }

  return residuals;
}

// Function to calculate weighted log-likelihood (posterior)
// [[Rcpp::export]]
double logLikPostCpp(NumericVector parameters, NumericVector x, NumericVector y, Function func, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  // Extract sigma (last parameter) and remaining params
  NumericVector params = parameters[seq(0, parameters.size() - 2)];
  double sigma = parameters[parameters.size() - 1];

  // Get the fitted values using the provided function
  NumericVector fitted_values = as<NumericVector>(func(params, x));

  // Initialize log-likelihood
  double loglik = 0.0;

  // Apply weights based on the method or the provided weights vector
  NumericVector w;
  if (weights.isNotNull()) {
    w = weights.get();
  } else {
    w = NumericVector(y.size(), 1.0);  // Default to 1 if no weights are provided
  }

  if (weight_method == "~y_sqrt") {
    w = sqrt(abs(y));  
  } else if (weight_method == "~y") {
    w = abs(y);  // Apply y as weights
  }

  // Calculate the log-likelihood assuming normal distribution with weights
  int n = y.size();
  for (int i = 0; i < n; i++) {
    double diff = y[i] - fitted_values[i];
    loglik += w[i] * (-0.5 * (log(2 * M_PI) + 2 * log(sigma) + pow(diff / sigma, 2)));
  }

  return loglik;
}

// Function to calculate weighted log-likelihood (posterior) with N and IEI
// [[Rcpp::export]]
double logLikPostCpp_N(NumericVector parameters, NumericVector x, NumericVector y, Function func, int N = 1, double IEI = 100, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  // Extract sigma (last parameter) and remaining params
  NumericVector params = parameters[seq(0, parameters.size() - 2)];
  double sigma = parameters[parameters.size() - 1];

  // Get the fitted values using the provided function with N and IEI
  NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));

  // Initialize log-likelihood
  double loglik = 0.0;

  // Apply weights based on the method or the provided weights vector
  NumericVector w;
  if (weights.isNotNull()) {
    w = weights.get();
  } else {
    w = NumericVector(y.size(), 1.0);  // Default to 1 if no weights are provided
  }

  if (weight_method == "~y_sqrt") {
    w = sqrt(abs(y));  // Apply sqrt of y as weights
  } else if (weight_method == "~y") {
    w = abs(y);  // Apply y^2 as weights
  }

  // Calculate the log-likelihood assuming normal distribution with weights
  int n = y.size();
  for (int i = 0; i < n; i++) {
    double diff = y[i] - fitted_values[i];
    loglik += w[i] * (-0.5 * (log(2 * M_PI) + 2 * log(sigma) + pow(diff / sigma, 2)));
  }

  return loglik;
}

// Function to calculate sum of squares
// [[Rcpp::export]]
double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI, Nullable<NumericVector> weights = R_NilValue) {
  NumericVector residuals = residFunCpp(params, y, x, func, N, IEI, weights);
  return sum(pow(residuals, 2));  // Sum of squared residuals
}
'

# Write the C++ code to a temporary file
cpp_file <- tempfile(fileext=".cpp")
writeLines(cpp_code, cpp_file)

# Source the C++ file
sourceCpp(cpp_file)

FITN <- function(response, dt=0.1, func=product2N, N=1, IEI=50, method=c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'), 
                 weight_method=c('none', '~y_sqrt', '~y'), stimulation_time=0, baseline=0, fast.decay.limit=NULL, latency.limit=NULL, 
                 lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
                 MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
                 MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'),
                 response_sign_method=c('smooth', 'regression', 'cumsum'), 
                 dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, 
                 return.output=FALSE, show.output=TRUE, show.plot=TRUE) {

  dx <- dt 
  method <- match.arg(method)
  weight_method <- match.arg(weight_method)
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  # sign <- sign_fun(y, direction_method=response_sign_method) 
  # y <- sign * y

  x <- seq(0, (length(y) - 1) * dx, by=dx)

  if (filter) {
    ind=20
    fs=1 / dx * 1000; bf <- butter(2, fc / (fs / 2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  x.orig <- x
  ind1 <- (stimulation_time - baseline) / dx
  ind2 <- baseline / dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by=dx)

  sign <- sign_fun(yfilter - mean(yfilter[1:ind2]), direction_method=response_sign_method) 
  y <- sign * y; yorig <- sign * yorig; yfilter <- sign * yfilter
  
  if (is_product_function(func)) {
    yorig <- yorig - mean(yorig[1:ind2])
    yfilter <- yfilter - mean(yfilter[1:ind2])

    y2fit <- yfilter[ind2:length(yfilter)]
    x2fit <- seq(0, dx * (length(y2fit) - 1), by=dx)

    upper <- adjust_product_bounds_N(bounds=upper, func=func, N=N, upper=TRUE)
    lower <- adjust_product_bounds_N(bounds=lower, func=func, N=N)

    if (is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit))) {
      if (identical(func, product1N)){
        upper <- c(rep(sign * Inf, N), rep(Inf, 3)) #  N amplitudes plus tau1, tau2, delay
      } else if (identical(func, product2N)) {
        upper <- rep(c(rep(sign * Inf, N), rep(Inf, 3)),2)
      } else if (identical(func, product3N)) {
        upper <- rep(c(rep(sign * Inf, N), rep(Inf, 3)),3)
      }
    
      if (!is.null(fast.decay.limit)){
        upper[N+2] <- fast.decay.limit[1]
        if (identical(func, product2N)){
          upper[2*(N+2)+1] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
        }
        if (identical(func, product3N)){
          upper[2*(N+2)+1] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
          upper[2*(N+2)+2] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[3]
        }
      }

      if (!is.null(latency.limit)){
          if (identical(func, product1N)){
            upper[N+3] <- latency.limit
          } else if (identical(func, product2N)){
            upper[c(N+3, 2*(N+3))] <- latency.limit
          } else if (identical(func, product3N)){
            upper[c(N+3, 2*(N+3), 3*(N+3))] <- latency.limit
          }
        }

    }

    # Define the formulas for each of the product functions
    if (identical(func, product1N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay')
      N.params <- length(param_names)
 
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) lower[1:N] <- sign * lower[1:N] 
      if (!is.null(upper)) upper[1:N] <- sign * upper[1:N] 
      
    } else if (identical(func, product2N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay', 
                       paste0('a', 1:N, '2'), 'tau1_2', 'tau2_2', 'delay_2')
      N.params <- length(param_names)
 
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1:N] <- sign * lower[1:N]; lower[(N+4):(2*N+3)] <- sign * lower[(N+4):(2*N+3)] }
      if (!is.null(upper)) { upper[1:N] <- sign * upper[1:N]; upper[(N+4):(2*N+3)] <- sign * upper[(N+4):(2*N+3)] }
     
    } else if (identical(func, product3N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay', 
                       paste0('a', 1:N, '2'), 'tau1_2', 'tau2_2', 'delay_2',
                       paste0('a', 1:N, '3'), 'tau1_3', 'tau2_3', 'delay_3')
      N.params <- length(param_names)

      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1:N] <- sign * lower[1:N]; lower[(N+4):(2*N+3)] <- sign * lower[(N+4):(2*N+3)]; lower[(2*N+7):(3*N+6)] <- sign * lower[(2*N+7):(3*N+6)] }
      if (!is.null(upper)) { upper[1:N] <- sign * upper[1:N]; upper[(N+4):(2*N+3)] <- sign * upper[(N+4):(2*N+3)]; upper[(2*N+7):(3*N+6)] <- sign * upper[(2*N+7):(3*N+6)] }

    }
  }else{
    fun.2.fit <- func
    y2fit <- yfilter[ind2:length(yfilter)]
    x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)
    param_names <- generate_param_names(fun.2.fit)
    N.params <- length(param_names)
    form.2.fit <- generate_formula(fun.2.fit, param_names, y_name = 'y', x_name = 'x')
  }


 if (method == 'MLE'){
    
    MLE.method <- match.arg(MLE.method)
    iter <- MLEsettings$iter
    metropolis.scale <- MLEsettings$metropolis.scale
    fit.attempts <- MLEsettings$fit.attempts
    RWm <- MLEsettings$RWm  
    sd.est <- sqrt(mean(yfilter[1:ind2]^2))

    logpost <- if (is_product_function(func)) logLikPostCpp_N else logLikPostCpp

    output <- FIT.MLE(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, sigma=sd.est, iter=iter, 
      metropolis.scale=metropolis.scale, logpost=logpost, params.start=NULL, 
      method=MLE.method, lower=lower, upper=upper, weight_method=weight_method, MLE.fun.attempts=100, fit.attempts=fit.attempts, RWm=RWm) 
  
  } else if (method == 'LM') {
    
    output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=NULL, upper=NULL, 
      weight_method=weight_method, fit.convergence.attempts=10, fit.attempts=10)
  
  } else if (method == 'BF.LM') {
    
    # default method
    bounds <- check_and_set_bounds(x=x2fit, y=y2fit, func=func, N.params=N.params, upper=upper, lower=lower)
    lower <- bounds$lower
    upper <- bounds$upper
    # print(upper)
    
    output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params,  N=N, IEI=IEI, lower=lower, upper=upper, 
      weight_method=weight_method, fit.convergence.attempts=10, fit.attempts=10)

  } else if (method == 'GN') {
    
    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=NULL, upper=NULL, weight_method=weight_method, algorithm='default')

  } else if (method == 'port') {
    
    bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
    upper <- bounds$upper
    lower <- if (method == 'port') bounds$lower else rep(-Inf, N.params)

    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, algorithm='port')
    
  } else if (method == 'robust') {
    
    bounds <- check_and_set_bounds(x=x2fit, y=y2fit, func=func, N.params=N.params, upper=upper, lower=lower)
    lower <- bounds$lower
    upper <- bounds$upper
    lat_win  <- if (is.null(latency.limit)) c(0.1, 10) else c(0.1, latency.limit)
    area_win <- c(0, max(x2fit))
    
    output <- FIT.NLSrobust(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, 
                        lower=lower, upper=upper, k=5.0, latency_window_ms=lat_win, area_window_ms=area_win)
  }
    
  if (is_product_function(func)) {
    
    if (identical(func, product1N)){
      
      df_output <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      output$fits[1:N] <- sign * output$fits[1:N]
    
    } else if (identical(func, product2N)){  
      df_output1 <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=sign)
            # Create a list of the outputs
      output_list <- list(df_output1, df_output2)
      # Extract the third elements from each output and determine the order
      order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
      # Reorder the output list based on the third element
      output_ordered <- output_list[order_indices]

      # Combine the outputs in the correct order
      df_output <- rbind('fast' = output_ordered[[1]], 'slow' = output_ordered[[2]])

      # Adjust the order of the fits and fits.se based on the sorted order
      output$fits <- c(output$fits[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])])

      output$fits.se <- c(output$fits.se[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits.se[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])])

      # adjust signs
      output$fits[1:N] <- sign * output$fits[1:N]
      output$fits[(N+4):(2*N+3)] <- sign * output$fits[(N+4):(2*N+3)]


    } else if (identical(func, product3N)){ 

      df_output1 <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=sign)
      df_output3 <- out.fun(params=output$fits[(2*N+7):(3*N+9)], interval=interval, dp=dp, sign=sign)

      # Create a list of the outputs
      output_list <- list(df_output1, df_output2, df_output3)
      # Extract the third elements from each output and determine the order
      order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
      # Reorder the output list based on the third element
      output_ordered <- output_list[order_indices]

      # Combine the outputs in the correct order
      df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])

      # Adjust the order of the fits and fits.se based on the sorted order
      output$fits <- c(output$fits[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])],
                       output$fits[((N+3)*order_indices[3]-(N+2)):((N+3)*order_indices[3])])

      output$fits.se <- c(output$fits.se[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits.se[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])],
                       output$fits.se[((N+3)*order_indices[3]-(N+2)):((N+3)*order_indices[3])])

      # adjust signs
      output$fits[1:N] <- sign * output$fits[1:N]
      output$fits[(N+4):(2*N+3)] <- sign * output$fits[(N+4):(2*N+3)]
      output$fits[(2*N+7):(3*N+6)] <- sign * output$fits[(2*N+7):(3*N+6)]
    }
  
    traces <- data.frame(x = xorig, y = sign * yorig, yfilter = sign * yfilter)
    fits <- output$fits

    # Define the list of functions
    func_list <- list(product1N, product2N, product3N)

    # Find the index of the matching function
    func_index <- sapply(func_list, function(f) identical(f, func))

    # Check if func matches one of the productN functions and apply the corresponding baseline adjustments
    if (any(func_index)) {
      index <- which(func_index)  # Get the index of the matching function
      for (i in 1:index) {
        fits[i * N + (3 * i)] <- fits[i * N + (3 * i)] + baseline
      }
    }

    traces$yfit <- func(params=fits, x=traces$x+dx, N=N, IEI=IEI) 

    if (identical(func, product2N)){
      traces$yfit1 <- product1N(params=fits[1:(N+3)], x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)], x=traces$x+dx, N=N, IEI=IEI) 
    } 
    if (identical(func, product3N)){
      traces$yfit1 <- product1N(params=fits[1:(N+3)],  x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)],  x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit3 <- product1N(params=fits[(2*N+7):(3*N+9)], x=traces$x+dx, N=N, IEI=IEI) 
    } 

  }else{
    traces=data.frame(x=xorig, y=sign * yorig, yfilter=sign * yfilter)
    fits <- output$fits
    traces$yfit <- func(fits,x2fit)  
    df_output <- data.frame(fits=fits)
    df_output <- t(df_output)

    df_output <- as.data.frame(df_output)
    colnames(df_output) <- param_names
  }

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  
  if (show.output){
    print(df_output)
  }

  if (return.output){
    return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces))
  }
}


eigen.values <- function(x, only.values = FALSE) {
  z <- .Internal(La_rg(x, only.values))
  ord <- order(Mod(z$values), decreasing = TRUE)
  z$values[ord]
}

residFun <- function(params, y, x, func) {
  y - func(params, x)
}

SS.fun <- function(params, x, y, func) {
  sum((y - func(params, x))^2)
}

# MLE function
MLE.fun <- function(logpost, par, x, y, N=1, IEI=100, gr=NULL, lower=NULL, upper=NULL, weight_method='none', method='Nelder-Mead', func=product1) {
  if (method %in% c('L-BFGS-B', 'Brent')) {
    suppressWarnings(
      MLE.fit <- optim(
        par=par, 
        fn=logpost, 
        gr=gr, 
        method=method, 
        hessian=TRUE, 
        lower=lower, 
        upper=upper, 
        control=list(fnscale=-1, maxit=2000, factr=1e7, pgtol=1e-8, parscale=rep(1, length(par))),
        x=x, 
        y=y, 
        func=func, 
        N=N,         
        IEI=IEI,
        weight_method=weight_method      
      )
    )
  } else {
    suppressWarnings(
      MLE.fit <- optim(
        par=par, 
        fn=logpost, 
        gr=gr, 
        method=method, 
        hessian=TRUE, 
        control=list(fnscale=-1, maxit=2000, reltol=1e-8, parscale=rep(1, length(par))),
        x=x, 
        y=y, 
        func=func, 
        N=N,         
        IEI=IEI,
        weight_method=weight_method      
      )
    )
  }

  fit <- MLE.fit$par
  h <- -solve(MLE.fit$hessian)
  p <- length(fit)
  int <- p / 2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(fit, x, y, func, N, IEI)
  list(fit=fit, fit.se=sqrt(diag(h)), var=h, int=int, converge=MLE.fit$convergence == 0)
}

# Optimized RWmetropolis function
RWmetropolis <- function(logpost, cov.mat, scale, start, m, x, y, func, N=1, IEI=50) {
  pb <- length(start)
  Mpar <- matrix(0, m, pb)
  b <- matrix(t(start))
  lb <- logpost(start, x, y, func)
  a <- chol(cov.mat)
  accept <- 0
  
  random_samples <- scale * t(a) %*% matrix(rnorm(pb * m), nrow = pb)
  
  for (ii in 1:m) {
    bc <- b + random_samples[, ii]
    lbc <- logpost(parameters=t(bc), x=x, y=y, func=func, N=N, IEI=IEI) # lbc <- logpost(t(bc), x, y, func)
    prob <- exp(lbc - lb)
    if (!is.na(prob) && runif(1) < prob) {
      lb <- lbc
      b <- bc
      accept <- accept + 1
    }
    Mpar[ii, ] <- b
  }
  list(par = Mpar, accept = accept / m)
}

# Log-likelihood function
log.lik.post <- function(parameters, x, y, func) {
  params <- parameters[-length(parameters)]
  sigma <- parameters[length(parameters)]
  sum(dnorm(y, mean=func(params, x), sd=sigma, log=TRUE))
}

log.lik.post_N <- function(parameters, x, y, func, N=1, IEI=100) {
  params <- parameters[-length(parameters)]
  sigma <- parameters[length(parameters)]
  sum(dnorm(y, mean=func(params=params, x=x, N=N, IEI=IEI), sd=sigma, log=TRUE))
}

# model.selection.criteria <- function(coeffs, x, y, func, N=1, IEI=100) {
#   res <- residFunCpp(coeffs, y, x, func, N, IEI)
#   k <- length(coeffs)
#   n <- length(res)
#   # definition: log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum(res^2)
#   # log(sigma2) == sum(res^2) / n and (1 / (2 * sigma2)) * sum(res^2) == 0.5 * n
#   loglik <- -0.5 * (n * log(2 * pi) + n * log(sum(res^2) / n) + n)
#   # loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sum(res^2))))
#   df <- k + 1
#   BIC <- df * log(n) - 2 * loglik
#   AIC <- df * 2 - 2 * loglik
#   c(AIC=AIC, BIC=BIC)
# }

# model.selection.criteria <- function(coeffs, x, y, func, N=1, IEI=100) {
#   res <- residFunCpp(coeffs, y, x, func, N, IEI)
#   k <- length(coeffs)
#   n <- length(res)
#   # definition: log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum(res^2)
#   # log(sigma2) == sum(res^2) / n and (1 / (2 * sigma2)) * sum(res^2) == 0.5 * n
#   loglik <- -0.5 * (n * log(2 * pi) + n * log(sum(res^2) / n) + n)
#   # loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sum(res^2))))
#   BIC <- k * log(n) - 2 * loglik
#   AIC <- k * 2 - 2 * loglik
#   c(AIC=AIC, BIC=BIC)
# }

model.selection.criteria <- function(coeffs, x, y, func, N=1, IEI=100, include.sigma=FALSE) {
  res <- residFunCpp(coeffs, y, x, func, N, IEI)
  k <- length(coeffs)
  n <- length(res)
  loglik <- -0.5 * (n * log(2 * pi) + n * log(sum(res^2) / n) + n)
  df <- if (include.sigma) k + 1 else k
  BIC <- df * log(n) - 2 * loglik
  AIC <- df * 2 - 2 * loglik
  c(AIC=AIC, BIC=BIC)
}

# Optimized fit.MLE function
fit.MLE <- function(x, y, func, N.params, N=1, IEI=100, sigma=5, iter=1e4, metropolis.scale=2, logpost=logLikPostCpp, 
  params.start=NULL, method='Nelder-Mead', lower=NULL, upper=NULL, weight_method='none', MLE.fun.attempts=100, RWm=TRUE) {
  
  if (is.null(params.start)) {
     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper)
     params.start <- c(start, sigma)
  }
  n <- length(y)
  k <- length(params.start) - 1
  attempts <- 0
  
  l.fit <- list(int=NaN, converge=FALSE)
  while (!(!is.nan(l.fit$int) && l.fit$converge) && (attempts < MLE.fun.attempts)) {
    st1 <- if (is_product_function(func)) params.start * runif(k + 1) else params.start
    l.fit <- suppressWarnings(
      MLE.fun(logpost=logpost, par=st1, x=x, y=y, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, method=method, func=func)
    )
    attempts <- attempts + 1
  }
  
  if (!l.fit$converge) stop('MLE.fun does not converge')
  
  e.values <- Re(eigen.values(l.fit$var))
  if (any(e.values <= 0)) stop('MLE.fun covariance matrix not positive definite')

  if (RWm) {
    rw.fit <- RWmetropolis(logpost=logpost, cov.mat=l.fit$var, scale=metropolis.scale, start=l.fit$fit, m=iter, x=x, y=y, func=func, N=N, IEI=IEI)
    
    acceptance.rate <- rw.fit$accept
    parameters <- apply(rw.fit$par, 2, median)
    
    fits <- parameters[1:k]
    fits.se <- apply(rw.fit$par, 2, sd)[1:k]
  } else {
    fits <- l.fit$fit[1:k]
    fits.se <- sqrt(diag(l.fit$var))[1:k]
  }
  
  paramsfit <- fits
  res <- residFunCpp(paramsfit, y, x, func, N, IEI)
  gof.se <- (sum(res^2) / (n - k))^0.5

  msc <- model.selection.criteria(coeffs=paramsfit, x=x, y=y, func, N=N, IEI=IEI)

  if (RWm) {
    list(fits=fits, fits.se=fits.se, acceptance.rate=acceptance.rate, MLE.fun.fit=l.fit$fit, MLE.fun.fit.se=l.fit$fit.se, MLE.fun.var=l.fit$var, MLE.fun.int=l.fit$int, MLE.fun.converge=l.fit$converge, MLE.fun.convergence.attempts=attempts, gof=gof.se, AIC=msc[1], BIC=msc[2])
  } else {
    list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=l.fit$converge, model.message=l.fit$converge)
  }
}

# Wrapper for fit.MLE
FIT.MLE <- function(x, y, func, N.params, N=1, IEI=100, sigma=5, iter=1e4, metropolis.scale=2, logpost=logLikPostCpp_N, params.start=NULL, 
  method='Nelder-Mead', lower=NULL, upper=NULL,  weight_method='none', MLE.fun.attempts=100, fit.attempts=10, RWm=TRUE) {

  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper # some solutions were sitting on this limit so increased 2-fold
  if (method %in% c('L-BFGS-B', 'Brent')) {
    tol_low <- 1e-08
    lower <- if (any(lower <= 0)) lower + tol_low
  }

  output.MLE <- NULL
  attempts <- 0
  while (is.null(output.MLE) && attempts < fit.attempts) {
    tryCatch({
      output.MLE <- fit.MLE(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, sigma=sigma, iter=iter, metropolis.scale=metropolis.scale, logpost=logpost, params.start=params.start, 
        method=method, lower=lower, upper=upper, weight_method=weight_method, MLE.fun.attempts=MLE.fun.attempts, RWm=RWm)
    }, error=function(e) {})
    attempts <- attempts + 1
  }
  c(output.MLE, fit.attempts=attempts)
}

# Function to try optim with N retries
optim.N <- function(N.params, SS.fun, lower, upper, y, x, func, max_attempts=10) {
  attempt <- 1
  success <- FALSE
  est <- NULL

  while (attempt <= max_attempts && !success) {
    tryCatch({
      est <- optim(runif(N.params), SS.fun, method='L-BFGS-B', lower=lower, upper=upper, hessian=TRUE, y=y, x=x, func=func)
      success <- TRUE
    }, error = function(e) {
      if (grepl("L-BFGS-B needs finite values of 'fn'", e$message)) {
        attempt <- attempt + 1
      } else {
        stop(e)
      }
    })
  }

  if (!success) {
    stop("Failed to find a valid solution with 'L-BFGS-B' method after ", max_attempts, " attempts.")
  }

  return(est)
}

# Function to try optim with N retries
start_optimiser <- function(N.params, SS.fun, lower, upper, y, x, func, max_attempts=10, start.method=c('uniform', 'lognormal'), cv=0.4){
  start.method <- match.arg(start.method)
  if (identical(func, product3) || identical(func, product3N)) {
    N <- (N.params - 9)/3
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars <- c(rep(ests[1], N), ests[2],  ests[3]) 
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    pars2 <- c(pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[(N+4):(2*(N+3))]) 
    pars3 <- c(pars, 15)
    if (!is.null(upper)) pars3 <- pmin(pars3, upper[(2*(N+3)+1):(3*(N+3))]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2, pars3) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=c(pars1, pars2, pars3), cv=cv, n=1)
    }

  } else if (identical(func, product2) || identical(func, product2N)) {
    N <- (N.params - 6)/2
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars <- c(rep(ests[1], N), ests[2],  ests[3]) 
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    pars2 <- c(pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[(N+4):(2*(N+3))]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=c(pars1, pars2), cv=cv, n=1)
    }

  } else if (identical(func, product1) || identical(func, product1N)) {
    N <- (N.params - 3)
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars1 <- c(rep(ests[1], N), ests[2],  ests[3], 15) 
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    if (start.method=='uniform'){
      par <- pars1 * runif(N.params)
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=pars1, cv=cv, n=1)
    }

  }
  est <- list(par=par)
  return(est)
}

# function to check if a given function is product1 or product2
is_product_function <- function(func) {
  valid_funcs <- list(product1, product2, product3, product1N, product2N, product3N)
  any(vapply(valid_funcs, identical, logical(1), func))
}

check_and_set_bounds <- function(x, y, func, N.params, upper=NULL, lower=NULL) {
  if (is_product_function(func)) {
    if (is.null(upper)) {
      if (identical(func, product2) || identical(func, product2N)) {
        N <- (N.params - 6)/2
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- rep(c(rep(ests[1], N), Inf,  ests[3], Inf) ,2) 
        upper <- sapply(6 * ests, round_up)
      } else if (identical(func, product1) || identical(func, product1N)) {
        N <- N.params - 3
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- c(rep(ests[1], N), Inf,  ests[3], Inf)
        upper <- sapply(6 * ests, round_up)
      } else if (identical(func, product3) || identical(func, product3N)){
        N <- (N.params - 9)/3
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- rep(c(rep(ests[1], N), Inf,  ests[3], Inf) ,3) 
        upper <- sapply(6 * ests, round_up)
      }  
    }
    if (is.null(lower)) lower <- rep(0, N.params)
  } else {
    upper <- rep(Inf, N.params)
    lower <- rep(-Inf, N.params)
  }
  return(list(upper=upper, lower=lower))
}

start_optimization <- function(x, y, func, N.params, SS.fun, lower, upper) {
  if (is_product_function(func)) {
    est <- start_optimiser(N.params = N.params, SS.fun = SS.fun, lower = lower, upper = upper, y = y, x = x, func = func, max_attempts = 10)
    start <- est$par
  } else {
    est <- optim(par = runif(N.params), SS.fun, method = "L-BFGS-B", lower= lower, upper = upper, hessian = TRUE, y = y, x = x, func = func) 
    start <- est$par
  }
  return(start)
}

fit.LM <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none') {
  # est <- optim(runif(N.params), SS.fun, method='L-BFGS-B', lower=lower, upper=upper, hessian=TRUE, y=y, x=x, func=func)
 
  start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
  model <- minpack.lm::nls.lm(
    par=start, 
    fn=residFunCpp, 
    y=y, 
    x=x, 
    lower=lower, 
    upper=upper, 
    control=list(maxiter=1024, tol=1e-05, warnOnly=TRUE),
    func=func,  
    N=N,             
    IEI=IEI,
    weight_method=weight_method     
  )
  e.values <- Re(eigen.values(model$hessian))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)
  list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$info, model.message=model$message)
}

round_up <- function(x, factor=10) factor * ceiling( x / factor)
  
# Wrapper for fit.LM
FIT.LM <- function(x, y, func, N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', fit.convergence.attempts=10, fit.attempts=10) {
  convergence.attempts <- 0
  modelinfo <- NULL

  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper
  
  while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.LM <- NULL
    attempts <- 0
    while (is.null(output.LM) && attempts < fit.attempts) {
      tryCatch({
        output.LM <- fit.LM(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method)
      }, error=function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.LM$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.LM, convergence.attempts=convergence.attempts)
}

# FIT.bf.LM <- function(x, y, func, N.params, N=1, IEI=100, lower=NULL, upper=NULL, fit.convergence.attempts = 10, fit.attempts = 10) {
#   convergence.attempts <- 0
#   modelinfo <- NULL

#   # bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
#   # lower <- bounds$lower
#   # upper <- bounds$upper

#   while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
#     output.LM <- NULL
#     attempts <- 0
#     while (is.null(output.LM) && attempts < fit.attempts) {
#       tryCatch({
#         output.LM <- fit.bf.LM(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper)
#       }, error = function(e) {})
#       attempts <- attempts + 1
#     }
#     modelinfo <- output.LM$model.info
#     convergence.attempts <- convergence.attempts + 1
#   }
#   c(output.LM, convergence.attempts = convergence.attempts)
# }

# # Wrapper for FIT.bf.LM
# FIT.BF.LM <- function(x, y, func, N.params, lower = NULL, upper = NULL, fit.convergence.attempts = 10, fit.attempts = 10, N.repeat = 25) {
#   fits <- matrix(NA, N.repeat, N.params)
#   fits.se <- matrix(NA, N.repeat, N.params)
#   gof.se <- rep(NA, N.repeat)
#   BIC <- rep(NA, N.repeat)
#   AIC <- rep(NA, N.repeat)
#   fit.convergence.attempts <- rep(NA, N.repeat)
#   model.info <- rep(NA, N.repeat)
#   model.message <- rep(NA, N.repeat)
#   convergence.attempts <- rep(NA, N.repeat)
  
#   for (iii in 1:N.repeat) {
#     output1 <- FIT.bf.LM(x=x, y=y, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=fit.convergence.attempts, fit.attempts=fit.attempts)
#     if (!is.null(output1)) {
#       fits[iii, ] <- output1$fits
#       fits.se[iii, ] <- output1$fits.se
#       gof.se[iii] <- output1$gof
#       AIC[iii] <- output1$AIC
#       BIC[iii] <- output1$BIC
#       fit.convergence.attempts[iii] <- output1$convergence.attempts
#       model.info[iii] <- output1$model.info
#       model.message[iii] <- output1$model.message
#       convergence.attempts[iii] <- output1$convergence.attempts
#     }
#   }
  
#   ind <- which.min(gof.se)
  
#   list(fits = fits[ind, ], fits.se = fits.se[ind, ], gof = gof.se[ind], AIC = AIC[ind], BIC = BIC[ind], model.info = model.info[ind], model.message = model.message[ind], convergence.attempts = convergence.attempts[ind])
# }

# # Algorithm = c('default', 'port')
# fit.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, algorithm='default') {

#   if (is_product_function(func)){
#     if (algorithm == 'port'){ 
#       start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
#     }else{          # default is Gauss-Newton
#       est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
#       start <- est$par
#     }
#   }else{
#     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   
#   }
  
#   # names(start) <- param.names[!param.names == 'x']
  
#   # must declare form.2.fit formula in local environment
  
#   form.2.fit <- as.formula(y ~ func(params, x, N, IEI))
  
#   model <- suppressWarnings(
#           nls(
#             form.2.fit, 
#             start=list(params = start),
#             lower=lower,
#             upper=upper,
#             algorithm=algorithm
#             ) 
#           )

#   e.values <- Re(eigen.values(summary(model)$cov.unscaled))
#   if (any(e.values <= 0)) stop('Hessian not positive definite')
  
#   fits <- summary(model)$coefficients[, 'Estimate']
#   fits.se <- summary(model)$coefficients[, 'Std. Error']
#   gof.se <- summary(model)$sigma
#   msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func, N=N, IEI=IEI)
#   list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$convInfo$stopCode, model.message=model$convInfo$stopMessage)
# }

# Algorithm = c('default', 'port')
fit.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', algorithm='default') {

  if (is_product_function(func)){
    if (algorithm == 'port'){ 
      start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
    }else{          # default is Gauss-Newton
      est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
      start <- est$par
    }
  }else{
    start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   
  }
  
  # weights
  if (weight_method == 'none') {
    weights <- rep(1, length(y)) 
  } else if (weight_method == '~y_sqrt') {
    weights <- abs(y)  
  } else if (weight_method == '~y') {
    weights <- y^2
  }
  
  # must declare form.2.fit formula in local environment
  
  form.2.fit <- as.formula(y ~ func(params, x, N, IEI))
  
  model <- suppressWarnings(
          nls(
            form.2.fit, 
            start=list(params = start),
            lower=lower,
            upper=upper,
            weights=weights,
            algorithm=algorithm
            ) 
          )

  e.values <- Re(eigen.values(summary(model)$cov.unscaled))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func, N=N, IEI=IEI)
  list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$convInfo$stopCode, model.message=model$convInfo$stopMessage)
}

# Wrapper for fit.NLS
FIT.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', algorithm='default', fit.convergence.attempts=10, fit.attempts=10) {
  convergence.attempts <- 0
  modelinfo <- NULL
  
  while ((is.null(modelinfo) || !any(modelinfo == c(0, 1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.NLS <- NULL
    attempts <- 0
    while (is.null(output.NLS) && attempts < fit.attempts) {
      tryCatch({
        output.NLS <- fit.NLS(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, algorithm=algorithm)
      }, error = function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.NLS$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.NLS, convergence.attempts = convergence.attempts)
}


# fit.NLSrobust <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, algorithm='default', method = 'M') {

#   # Define the formula, but only in terms of params and x
#   func_wrapper <- function(params, x) {
#     func(params, x, N, IEI)  # Pass N and IEI directly here
#   }

#   if (method == 'M') {
#     if (is_product_function(func)){
#       start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
    
#       form.2.fit <- as.formula(y ~ func_wrapper(params, x))

#     }else{   
#       est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
#       start <- est$par
#       form.2.fit <- as.formula(y ~ func_wrapper(params, x))
#     }


#     control=robustbase:::nlrob.control(method)

#     model <- suppressWarnings(
#           nlrob(
#             formula=form.2.fit, 
#             start = list(params = start),
#             data=data.frame(x=x, y=y),
#             algorithm=algorithm, 
#             method=method, 
#             control=control, 
#             lower=lower, 
#             upper=upper
#             )
#           )  

#   }else{
#     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   

#     model <- suppressWarnings(
#           robustbase::nlrob(
#             formula=form.2.fit, 
#             data=data.frame(x=x, y=y),
#             algorithm=algorithm, 
#             control=control, 
#             lower=lower, 
#             upper=upper
#             )
#           )  

#     }

#   e.values <- Re(eigen.values(vcov(model)))
#   if (any(e.values <= 0)) {
#     warning('Hessian not positive definite - using anyway')
#     # Don't stop, just continue
#   }
  
#   fits <- summary(model)$coefficients[, 'Estimate']
#   fits.se <- summary(model)$coefficients[, 'Std. Error']
#   gof.se <- summary(model)$Scale
#   msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)
#   list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.message=summary(model)$status)
# }

# fit.NLSrobust <- function(x, y, func, N.params, N=1, IEI=100,
#                           lower=NULL, upper=NULL,
#                           algorithm='default', method='M',
#                           params.start=NULL) {

#   func_wrapper <- function(params, x) {
#     func(params, x, N, IEI)
#   }

#   form.2.fit <- as.formula(y ~ func_wrapper(params, x))

#   control <- robustbase:::nlrob.control(method)

#   if (!is.null(params.start)) {
#     start <- params.start
#   } else if (is_product_function(func)) {
#     start <- start_optimization(x=x, y=y, func=func, N.params=N.params,
#                                 SS.fun=SS.fun, lower=lower, upper=upper)
#   } else {
#     est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params),
#                            upper=upper, y=y, x=x, func=func, max_attempts=10)
#     start <- est$par
#   }

#   model <- suppressWarnings(
#     robustbase::nlrob(
#       formula   = form.2.fit,
#       start     = list(params = start),
#       data      = data.frame(x=x, y=y),
#       algorithm = algorithm,
#       method    = method,
#       control   = control,
#       lower     = lower,
#       upper     = upper
#     )
#   )

#   e.values <- Re(eigen.values(vcov(model)))
#   if (any(e.values <= 0)) stop('Hessian not positive definite')

#   fits    <- summary(model)$coefficients[, 'Estimate']
#   fits.se <- summary(model)$coefficients[, 'Std. Error']
#   gof.se  <- summary(model)$Scale
#   msc     <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)

#   list(
#     fits          = unname(fits),
#     fits.se       = unname(fits.se),
#     gof           = gof.se,
#     AIC           = msc[1],
#     BIC           = msc[2],
#     model.message = summary(model)$status
#   )
# }

# # Wrapper for fit.NLSrobust
# FIT.NLSrobust <- function(x, y, func, N.params, N=1, IEI=100,
#                           lower=NULL, upper=NULL, algorithm='default',
#                           method='M', params.start=NULL,
#                           fit.convergence.attempts=10, fit.attempts=10) {

#   convergence.attempts <- 0
#   modelinfo <- NULL

#   while ((is.null(modelinfo) || modelinfo != 'converged') &&
#          convergence.attempts < fit.convergence.attempts) {
#     output.NLS <- NULL
#     attempts <- 0
#     while (is.null(output.NLS) && attempts < fit.attempts) {
#       output.NLS <- tryCatch({
#         fit.NLSrobust(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI,
#                       lower=lower, upper=upper, algorithm=algorithm,
#                       method=method, params.start=params.start)
#       }, error = function(e) {
#         NULL
#       })
#       # only use params.start on the first attempt; after that randomise
#       params.start <- NULL
#       attempts <- attempts + 1
#     }
#     modelinfo <- output.NLS$model.message
#     convergence.attempts <- convergence.attempts + 1
#   }
#   c(output.NLS, convergence.attempts = convergence.attempts)
# }

psi.huber <- function(u, k = 1.345, deriv = 0) {
  if (!deriv) return(pmin(1, k / abs(u)))
  abs(u) <= k
}

fit.NLSrobust.lm <- function(x, y, func, N.params, N=1, IEI=100,
                             lower=NULL, upper=NULL, params.start=NULL,
                             psi=psi.huber, k=2.0, max.iter=10, tol=1e-5,
                             latency_window_ms=c(1,10), area_window_ms=c(0,40),
                             protect_floor=0.4) {

  if (!is.null(params.start)) {
    start <- params.start
  } else {
    lm_fit <- fit.LM(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI,
                     lower=lower, upper=upper, weight_method='none')
    start <- lm_fit$fits
  }

  is_signal <- (x >= latency_window_ms[1]) & (x <= area_window_ms[2])
  if (!any(is_signal)) is_signal <- rep(TRUE, length(x))

  params <- start
  w <- rep(1, length(y))

  for (iter in 1:max.iter) {

    resid <- y - func(params, x, N, IEI)
    sigma <- mad(resid)
    if (sigma < .Machine$double.eps) sigma <- 1e-6
    r <- resid / sigma

    w.new <- psi(r, k = k) / pmax(abs(r), .Machine$double.eps)
    w.new[!is.finite(w.new)] <- 1
    w.new[is_signal] <- pmax(w.new[is_signal], protect_floor)

    model <- minpack.lm::nls.lm(
      par = params,
      fn = function(params, y, x, func, N, IEI, w) {
        sqrt(w) * (y - func(params, x, N, IEI))
      },
      y = y, x = x, func = func, N = N, IEI = IEI, w = w.new,
      lower = lower, upper = upper,
      control = list(maxiter = 200, tol = tol)
    )

    params.new <- model$par
    if (max(abs(params.new - params), na.rm = TRUE) < tol) {
      params <- params.new
      w <- w.new
      break
    }

    params <- params.new
    w <- w.new
  }

  fits <- params
  resid <- y - func(fits, x, N, IEI)
  sigma <- mad(resid)
  if (sigma < .Machine$double.eps) sigma <- 1e-6

  fits.se <- rep(NA_real_, length(fits))
  if (!is.null(model$hessian)) {
    e.values <- Re(eigen.values(model$hessian))
    if (all(is.finite(e.values)) && all(e.values > 0)) {
      fits.se <- sqrt(diag(solve(model$hessian))) * sigma
    }
  }

  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)

  list(fits = unname(fits), fits.se = unname(fits.se), gof = sigma,
       AIC = msc[1], BIC = msc[2], model.message = "converged",
       weights = w)
}

FIT.NLSrobust <- function(x, y, func, N.params, N=1, IEI=100,
                          lower=NULL, upper=NULL, params.start=NULL,
                          psi=psi.huber, k=1.345, max.iter=10, tol=1e-5,
                          latency_window_ms=c(1,10), area_window_ms=c(0,40),
                          protect_floor=0.4,
                          fit.convergence.attempts=10, fit.attempts=10) {

  convergence.attempts <- 0
  modelinfo <- NULL

  while ((is.null(modelinfo) || modelinfo != 'converged') &&
         convergence.attempts < fit.convergence.attempts) {

    output.NLS <- NULL
    attempts <- 0

    while (is.null(output.NLS) && attempts < fit.attempts) {

      output.NLS <- tryCatch({
        fit.NLSrobust.lm(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI,
                         lower=lower, upper=upper, params.start=params.start,
                         psi=psi, k=k, max.iter=max.iter, tol=tol,
                         latency_window_ms=latency_window_ms, area_window_ms=area_window_ms,
                         protect_floor=protect_floor)
      }, error = function(e) {
        NULL
      })

      params.start <- NULL
      attempts <- attempts + 1
    }

    modelinfo <- output.NLS$model.message
    convergence.attempts <- convergence.attempts + 1
  }

  c(output.NLS, convergence.attempts = convergence.attempts)
}

sort2.fun <- function(x){
  index = seq(1:length(x))
  if (any(is.na(x))) {
    return(x)
  }
  if (x[4] > x[2]) {
    x <- x[c(3,4,1,2,index[5:length(x)])]
  }
  x
}

sort2 <- function(mat, col.names=NULL){
  mat <- t(apply(mat, 1, function(row) {
    if (any(is.na(row))) {
      return(row)
    } else {
      return(sort2.fun(row))
    }
  }))
  colnames(mat) <- col.names
  mat
}

sort3.fun <- function(x){
  index=seq(1:length(x))
  if ( (x[2] > x[4]) && (x[4] > x[6]) ) {
    x <- x[c(1,2,3,4,5,6,7)]
  } else if ( (x[2] > x[6]) && (x[6] > x[4]) ) {
    x <- x[c(1,2,5,6,3,4,7)]
  } else if ( (x[4] > x[2]) && (x[2] > x[6]) ) {
    x <- x[c(3,4,1,2,5,6,7)]
  } else if ( (x[4] > x[6]) && (x[6] > x[2]) ) {
    x <- x[c(3,4,5,6,1,2,7)]
  } else if ( (x[6] > x[2]) && (x[2] > x[4]) ) {
    x <- x[c(5,6,1,2,3,4,7)]
  } else if ( (x[6] > x[4]) && (x[4] > x[2]) ) {
    x <- x[c(5,6,3,4,1,2,7)]
  }
  x
}

sort3 <- function(mat, col.names=NULL){
  mat <- t(apply(mat,1,sort3.fun))
  colnames(mat) <- col.names
  mat
  }

generate_formula <- function(fun, param_names, y_name = 'yfilter', x_name = 'x') {
  # Extract the body of the function as a character string
  fun_body <- deparse(body(fun))
  
  # Replace 'params' with the corresponding parameter names
  for (i in seq_along(param_names)) {
    fun_body <- gsub(paste0('params\\[', i, '\\]'), param_names[i], fun_body)
  }
  
  # Remove the function declaration and braces
  fun_body <- gsub('\\{', '', fun_body)
  fun_body <- gsub('\\}', '', fun_body)
  fun_body <- trimws(fun_body)
  fun_body <- paste(fun_body, collapse = ' ')
  
  # Create the formula
  formula_str <- paste(y_name, '~', fun_body)
  formula <- as.formula(formula_str)
  
  return(formula)
}

# generate parameter names
generate_param_names <- function(fun) {
  # Get the function body
  fun_body <- body(fun)
  
  # Convert the function body to a character vector
  fun_body_char <- deparse(fun_body)
  
  # Find all the params indices used in the function
  param_indices <- regmatches(fun_body_char, gregexpr("params\\[[0-9]+\\]", fun_body_char))
  param_indices <- unlist(param_indices)
  param_indices <- as.numeric(gsub("params\\[([0-9]+)\\]", "\\1", param_indices))
  
  # Handle any NAs or duplicates
  param_indices <- unique(param_indices[!is.na(param_indices)])
  
  # Generate names for the parameters
  param_names <- letters[seq_along(param_indices)]
  
  # Return a named vector of parameter names
  setNames(param_names, paste0("params[", param_indices, "]"))
}

tau_rise <- function(tau1, tau2) tau1 * tau2 / ( tau1 + tau2 )

tau1_fun <- function(tau_rise, tau_decay) tau_rise * tau_decay / ( tau_decay - tau_rise )

# Define the function to solve
find_tau1 <- function(tpeak, tau2) {
  # Define the equation to solve
  equation <- function(tau1) {
    tau1 * (exp(tpeak / tau1) - 1) - tau2
  }
  
  # Use uniroot to solve for tau1, give an interval of reasonable guesses
  result <- uniroot(equation, interval = c(1e-6, 100), tol = 1e-9)
  
  # Return the root found for tau1
  return(result$root)
}

product_area <- function(Apeak, tau1, tau2){
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2) 
  abs(Apeak / f * tau2^2 / (tau1 + tau2))
}

product_tpeak <- function(tau1, tau2) tau1 * log( ( tau1 + tau2 ) / tau1 )

# calculate interval 10-90% rise and 90-10% decay based on default interval <- c(0.1, 0.9) from tau1 and tau2
product_rise_and_decay_percent <- function(tau1, tau2, interval=c(0.1, 0.9), showplot=FALSE, verbose=FALSE) {
  # Calculate tpeak
  tpeak <- product_tpeak(tau1, tau2)

  # Calculate f
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)

  # Define the function to find the root of
  target_function <- function(x, p) {
    exp(-x / tau2) - exp(-x * (tau1 + tau2) / (tau1 * tau2)) - f * p
  }

  # Function to find roots
  find_roots <- function(p, end_time) {
    # Define the intervals based on the plot
    interval1 <- c(0, tpeak)  # Interval before the peak
    interval2 <- c(tpeak, end_time)  # Interval after the peak
    
    # Find the roots
    root1 <- tryCatch(uniroot(function(x) target_function(x, p), interval1)$root, 
                      error = function(e) NA)
    root2 <- tryCatch(uniroot(function(x) target_function(x, p), interval2)$root, 
                      error = function(e) NA)
    
    if (is.na(root1) | is.na(root2)) {
      if (verbose) {
        cat("Failed to find roots for p =", p, "\n")
        cat("Interval1:", interval1, "\n")
        cat("Interval2:", interval2, "\n")
      }
    }
    
    return(c(root1, root2))
  }

  # determine the end time of the signal
  end_time <- 10 * max(tau1, tau2) #  adjust based on your signal's characteristics

  # roots for both p values 
  if (showplot) {
    x <- seq(0, end_time, length.out = 1000)
    yplot <- exp(-x / tau2) - exp(-x * (tau1 + tau2) / (tau1 * tau2))

    plot(x, yplot, type = 'l', xlab = 'x', ylab = 'f(x)', bty='l', las=1, main = '')
    abline(h = interval[1] * f, col = 'Indianred', lty = 3)
    abline(h = interval[2] * f, col = 'Indianred', lty = 3)

    # Add text labels
    text(x = max(x) * 0.95, y = interval[1] * f, labels = 'p=0.1', pos = 3, col = 'Indianred')
    text(x = max(x) * 0.95, y = interval[2] * f, labels = 'p=0.9', pos = 3, col = 'Indianred')
  }

  # Find roots for both p values
  roots_lower <- find_roots(interval[1], end_time)
  roots_upper <- find_roots(interval[2], end_time)

  # Compute the absolute differences
  diffs <- abs(roots_lower - roots_upper)

  return(diffs)
}

# product1 <- function(params, x) {
#   a1_max <- params[1]
#   tau_rise <- params[2]
#   tau2 <- params[3]
#   delay <- params[4]
  
#   tau1 <- tau1_fun(tau_rise, tau2) # tau2 = tau_decay

#   x_adjusted <- pmax(0, x - delay)
#   f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)
#   a1 <- a1_max / f
#   a1 * (1 - exp(-x_adjusted / tau1)) * exp(-x_adjusted / tau2)
# }

product1 <- function(params, x) {
  a1_max <- params[1]
  tau1 <- params[2]
  tau2 <- params[3]
  delay <- params[4]
  
  x_adjusted <- pmax(0, x - delay)
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)
  a1 <- a1_max / f
  a1 * (1 - exp(-x_adjusted / tau1)) * exp(-x_adjusted / tau2)
}

product2 <- function(params, x) {
  product1(params[1:4], x) + product1(params[5:8], x)
}

# product1N <- function(params, x, N=1, IEI=50) {
#   # Extract parameters
#   amplitudes <- params[1:N]  
#   shared_params <- params[(N + 1):length(params)]  
#   out <- sapply(1:N, function(ii){
#     shifted_x <- x - (ii - 1) * IEI 
#     current_params <- c(amplitudes[ii], shared_params)  
#     product1(current_params, shifted_x)  
#     })
#   response <- apply(out, 1, sum) 
#   return(response)
# }

# product2N <- function(params, x, N=1, IEI=50) {
  
#   response1 <- product1N(params[1:(N+3)], x=x, N=N, IEI=IEI)
#   response2 <- product1N(params[(N+4):(2*N+6)], x=x, N=N, IEI=IEI)
#   response <- response1 + response2
  
#   return(response) 

# }

product1N <- function(params, x, N=1, IEI=50) {
  # Extract parameters
  amplitudes <- params[1:N]  
  shared_params <- params[(N + 1):length(params)]  
  
  response <- numeric(length(x))  
  
  # Loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    current_params <- c(amplitudes[i], shared_params)  
    current_response <- product1(current_params, shifted_x)  
    response <- response + current_response  
  }
  
  return(response)
}

product2N <- function(params, x, N=1, IEI=50) {
  # parameters for the first product1
  amplitudes1 <- params[1:N]                
  tau1_1 <- params[N + 1]                   
  tau2_1 <- params[N + 2]                   
  delay1 <- params[N + 3]                   
  
  # parameters for the second product1
  amplitudes2 <- params[(N + 4):(2 * N + 3)]  
  tau1_2 <- params[2 * N + 4]                
  tau2_2 <- params[2 * N + 5]                
  delay2 <- params[2 * N + 6]                 
  
  response <- numeric(length(x)) 
  
  # loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    
    # first product1
    current_params1 <- c(amplitudes1[i], tau1_1, tau2_1, delay1)
    current_response1 <- product1(current_params1, shifted_x)
    
    # second product1
    current_params2 <- c(amplitudes2[i], tau1_2, tau2_2, delay2)
    current_response2 <- product1(current_params2, shifted_x)
    
    response <- response + current_response1 + current_response2
  }
  
  return(response) 
}

product3 <- function(params, x) {
  product1(params[1:4], x) + product1(params[5:8], x) + product1(params[9:12], x)
}

product3N <- function(params, x, N=1, IEI=50) {
  # parameters for the first product1
  amplitudes1 <- params[1:N]   
  tau1_1 <- params[N + 1]
  tau2_1 <- params[N + 2]
  delay1 <- params[N + 3]
  
  # parameters for the second product1
  amplitudes2 <- params[(N + 4):(2 * N + 3)] 
  tau1_2 <- params[2 * N + 4]                
  tau2_2 <- params[2 * N + 5]                
  delay2 <- params[2 * N + 6]                
  
  # parameters for the third product1
  amplitudes3 <- params[(2 * N + 7):(3 * N + 6)] 
  tau1_3 <- params[3 * N + 7]                    
  tau2_3 <- params[3 * N + 8]                    
  delay3 <- params[3 * N + 9]                    
  
  response <- numeric(length(x))  
  
  # loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    
    # first product1
    current_params1 <- c(amplitudes1[i], tau1_1, tau2_1, delay1)
    current_response1 <- product1(current_params1, shifted_x)
    
    # second product1
    current_params2 <- c(amplitudes2[i], tau1_2, tau2_2, delay2)
    current_response2 <- product1(current_params2, shifted_x)
    
    # third product1
    current_params3 <- c(amplitudes3[i], tau1_3, tau2_3, delay3)
    current_response3 <- product1(current_params3, shifted_x)
    
    response <- response + current_response1 + current_response2 + current_response3
  }
  
  return(response)  # Return the total response
}

# sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k=5) {
#   method <- match.arg(direction_method)
  
#   # Check if the length of y is sufficient for processing
#   if (length(y) < 10) {
#     stop("The length of y must be at least 10.")
#   }
  
#   n <- length(y)
#   peak_value <- NA
  
#   if (method == 'smooth') {
#     # Calculate the smoothed signal using a simple moving average
#     smoothed_signal <- rep(NA, n)
#     for (i in 1:(n - k + 1)) {
#       smoothed_signal[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
#     }
#     # Find the peak value
#     peak_value <- max(abs(smoothed_signal), na.rm = TRUE)
#     peak_value <- smoothed_signal[which.max(abs(smoothed_signal))]
    
#   } else if (method == 'diff') {
#     # Calculate the differences of the signal
#     diff_signal <- diff(y)
#     # Identify both the maximum and minimum differences
#     max_diff <- max(diff_signal, na.rm = TRUE)
#     min_diff <- min(diff_signal, na.rm = TRUE)
    
#     # Determine the peak value considering the direction
#     peak_value <- ifelse(abs(max_diff) > abs(min_diff), max_diff, min_diff)
     
#   } else if (method == 'regression') {
#     # Fit a quadratic regression to the signal
#     x <- 1:n
#     model <- lm(y ~ I(x^2) + x)
    
#     # Use the fitted values to find the peak
#     fitted_values <- predict(model)
#     peak_value <- max(abs(fitted_values), na.rm = TRUE)
#     peak_value <- fitted_values[which.max(abs(fitted_values))]
    
#   } else if (method == 'cumsum') {
#     # Calculate the cumulative sum of the signal
#     cumsum_signal <- cumsum(y)
#     # Find the peak value of the cumulative sum
#     peak_value <- max(abs(cumsum_signal), na.rm = TRUE)
#     peak_value <- cumsum_signal[which.max(abs(cumsum_signal))]
#   }
  
#   # Determine if the peak is positive or negative
#   peak_direction <- ifelse(peak_value > 0, 1, -1)
  
#   # Return the sign of the peak direction
#   return(peak_direction[[1]])
# }


sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k = 5) {
  
  method <- match.arg(direction_method)
  
  # Check if the length of y is sufficient for processing
  if (length(y) < 10) {
    stop("The length of y must be at least 10.")
  }
  
  n <- length(y)
  peak_value <- NA
  
  if (method == 'smooth') {
    # Calculate the smoothed signal using a simple moving average
    smoothed_signal <- rep(NA, n)
    for (i in 1:(n - k + 1)) {
      smoothed_signal[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
    }
    # Remove NA values before finding the indices of the maximum and minimum values
    valid_indices <- which(!is.na(smoothed_signal))
    max_idx <- valid_indices[which.max(smoothed_signal[valid_indices])]
    min_idx <- valid_indices[which.min(smoothed_signal[valid_indices])]
    
    # Determine the peak value based on the magnitude comparison
    if (abs(smoothed_signal[max_idx]) > abs(smoothed_signal[min_idx])) {
      peak_value <- smoothed_signal[max_idx]
    } else {
      peak_value <- smoothed_signal[min_idx]
    }
    
  } else if (method == 'diff') {
    # Calculate the differences of the signal
    diff_signal <- diff(y)
    # Identify both the maximum and minimum differences
    max_diff <- max(diff_signal, na.rm = TRUE)
    min_diff <- min(diff_signal, na.rm = TRUE)
    
    # Determine the peak value considering the direction
    peak_value <- ifelse(abs(max_diff) > abs(min_diff), max_diff, min_diff)
     
  } else if (method == 'regression') {
    # Fit a quadratic regression to the signal
    x <- 1:n
    model <- lm(y ~ I(x^2) + x)
    
    # Use the fitted values to find the peak
    fitted_values <- predict(model)
    max_idx <- which.max(fitted_values)
    min_idx <- which.min(fitted_values)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(fitted_values[max_idx]) > abs(fitted_values[min_idx])) {
      peak_value <- fitted_values[max_idx]
    } else {
      peak_value <- fitted_values[min_idx]
    }
    
  } else if (method == 'cumsum') {
    # Calculate the cumulative sum of the signal
    cumsum_signal <- cumsum(y)
    # Find the indices of the maximum and minimum values
    max_idx <- which.max(cumsum_signal)
    min_idx <- which.min(cumsum_signal)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(cumsum_signal[max_idx]) > abs(cumsum_signal[min_idx])) {
      peak_value <- cumsum_signal[max_idx]
    } else {
      peak_value <- cumsum_signal[min_idx]
    }
  }
  
  # Determine if the peak is positive or negative
  peak_direction <- ifelse(peak_value > 0, 1, -1)
  
  # Return the sign of the peak direction
  return(peak_direction)
}

out.fun <- function(params, interval = c(0.1, 0.9), dp = 3, sign=1) {
  
  # Calculate derived metrics
  N <- length(params)-3

  A <- sign*params[1:N]
  tau.rise <- tau_rise(params[N+1], params[N+2])
  tau.decay <- params[N+2]
  tpeak <- product_tpeak(params[N+1], params[N+2])
  
  diffs <- product_rise_and_decay_percent(tau1 = params[N+1], tau2 = params[N+2], interval = interval, showplot = FALSE, verbose=FALSE)
  trise.percent <- diffs[1]
  tdecay.percent <- diffs[2]
  
  area <- sapply(1:N, function(ii) product_area(params[ii], params[N+1], params[N+2]))
  delay <- params[[N+3]]
  
  # Calculate rise and decay interval labels
  r_percent_start <- interval[1] * 100
  r_percent_end <- interval[2] * 100
  r_label <- paste0('r', r_percent_start, '_', r_percent_end)
  
  d_label <- paste0('d', interval[2] * 100, '_', interval[1] * 100)
  
  # Create data frame with calculated values
  output <- c(round(A, dp), round(tau.rise, dp), round(tau.decay, dp), round(tpeak, dp), round(trise.percent, dp), round(tdecay.percent, dp), round(delay, dp), round(area, dp))
  col.names <- c(paste0('A', 1:N), 'τrise', 'τdecay', 'tpeak', r_label, d_label, 'delay', paste0('area', 1:N))
  names(output) <- col.names
  df_output <- data.frame(t(output))
  
  return(df_output)
}
  
adjust_product_bounds <- function(bounds, func, upper=FALSE) {
  tau1_fun <- function(tau_rise, tau_decay) {
    return(tau_rise * tau_decay / (tau_decay - tau_rise))
  }
  adjust_tau <- function(value1, value2, upper=FALSE) {
    tau <- tau1_fun(value1, value2)
    if (is.nan(tau)) {
      tau <- if (upper) Inf else 0
    } else if (!upper && tau == Inf) {
      tau <- 0
    }
    return(tau)
  }
  
  if (!is.null(bounds)) {
    if (identical(func, product1)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
    } else if (identical(func, product2)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6]  <- adjust_tau(bounds[6], bounds[7], upper=upper)
    } else if (identical(func, product3)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6]  <- adjust_tau(bounds[6], bounds[7], upper=upper)
      bounds[10] <- adjust_tau(bounds[10], bounds[11], upper=upper)
    }
  }
  return(bounds)
}

# adjust_product_bounds_N <- function(bounds, func, N, upper=FALSE) {
#   # Function to calculate tau1 based on tau_rise and tau_decay
#   tau1_fun <- function(tau_rise, tau_decay) {
#     return(tau_rise * tau_decay / (tau_decay - tau_rise))
#   }
  
#   # Function to adjust tau1 and handle upper or lower bounds
#   adjust_tau <- function(value1, value2, upper=FALSE) {
#     tau <- tau1_fun(value1, value2)
#     if (is.nan(tau)) {
#       tau <- if (upper) Inf else 0
#     } else if (!upper && tau == Inf) {
#       tau <- 0
#     }
#     return(tau)
#   }
  
#   # Adjust bounds based on the product function type
#   if (!is.null(bounds)) {
#     if (identical(func, product1N)) {
#       for (i in 1:length(bounds[1:N])) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#       }
#     } else if (identical(func, product2N)) {
#       for (i in 1:N) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#         bounds[(N + 5) + (i - 1) * 3] <- adjust_tau(bounds[(N + 5) + (i - 1) * 3], bounds[(N + 6) + (i - 1) * 3], upper=upper)
#       }
#     } else if (identical(func, product3N)) {
#       for (i in 1:N) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#         bounds[(N + 5) + (i - 1) * 3] <- adjust_tau(bounds[(N + 5) + (i - 1) * 3], bounds[(N + 6) + (i - 1) * 3], upper=upper)
#         bounds[(2 * N + 7) + (i - 1) * 3] <- adjust_tau(bounds[(2 * N + 7) + (i - 1) * 3], bounds[(2 * N + 8) + (i - 1) * 3], upper=upper)
#       }
#     }
#   }
  
#   return(bounds)
# }


adjust_product_bounds_N <- function(bounds, func, N, upper = FALSE) {
  # Function to calculate tau1 based on tau_rise and tau_decay
  tau1_fun <- function(tau_rise, tau_decay) {
    tau_rise * tau_decay / (tau_decay - tau_rise)
  }

  # Function to adjust tau1 and handle upper or lower bounds
  adjust_tau <- function(tau_rise, tau_decay, upper = FALSE) {
    tau <- tau1_fun(tau_rise, tau_decay)
    if (is.nan(tau) || tau < 0) {
      tau <- if (upper) Inf else 0
    }
    tau
  }

  # Adjust bounds based on the product function type
  if (!is.null(bounds)) {
    if (identical(func, product1N)) {
      # Adjust tau1
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
    } else if (identical(func, product2N)) {
      # Adjust tau1 for first component
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
      # Adjust tau1 for second component
      bounds[2 * N + 4] <- adjust_tau(bounds[2 * N + 4], bounds[2 * N + 5], upper = upper)
    } else if (identical(func, product3N)) {
      # Adjust tau1 for first component
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
      # Adjust tau1 for second component
      bounds[2 * N + 4] <- adjust_tau(bounds[2 * N + 4], bounds[2 * N + 5], upper = upper)
      # Adjust tau1 for third component
      bounds[3 * N + 7] <- adjust_tau(bounds[3 * N + 7], bounds[3 * N + 8], upper = upper)
    }
  }

  return(bounds)
}

# FIT <- function(response, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
#   stimulation_time=0, baseline=0, fast.decay.limit=NULL, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
#   MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
#   MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'),
#   response_sign_method = c('smooth', 'regression', 'cumsum'), 
#   dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, 
#   return.output=FALSE, show.output=TRUE, show.plot=TRUE){

#   dx <- dt 
#   method <- match.arg(method)
#   y <- response
#   if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
#     y <- y[!is.na(y)]
#   }

#   sign <- sign_fun(y, direction_method=response_sign_method) 
#   y <- sign * y

#   x <- seq(0, (length(y) - 1) * dx, by = dx)

#   if (filter){
#     ind = 20
#     fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
#     yfilter <- signal::filter(bf, y)
#   } else {
#     ind=1
#     yfilter=y
#   }
#   x.orig <- x

#   ind1 <- (stimulation_time - baseline)/dx
#   ind2 <- baseline/dx
  
#   yorig <- y[ind1:length(y)]
#   yfilter <- yfilter[ind1:length(yfilter)]
#   xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

#   if (is_product_function(func)) {
#     yorig <- yorig - mean(yorig[1:ind2])
#     yfilter <- yfilter - mean(yfilter[1:ind2])

#     y2fit <- yfilter[ind2:length(yfilter)]
#     x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)

#     # check if both upper and fast.decay.limit are specified
#     if ( !is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit)) ) {
#       warning("both 'upper' boundaries and at either 'fast decay limit' and/or 'latency limit' are specified:\ninputs 'fast.decay.limit' and/or 'latency.limit' will be ignored")
#       fast.decay.limit <- NULL
#       latency.limit <- NULL
#     }

#     # adjusts bounds for tau_rise to correct form involving tau1
#     upper <- adjust_product_bounds(bounds=upper, func=func, upper=TRUE)
#     lower <- adjust_product_bounds(bounds=lower, func=func)
      
#     if (is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit))){

#       if (identical(func, product1)){
#         upper <- c(sign * Inf, Inf, Inf, Inf) 
#       } else if (identical(func, product2)) {
#         upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
#       } else if (identical(func, product3)) {
#         upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
#       }

#       if (!is.null(fast.decay.limit)){
#         upper[3] <- fast.decay.limit[1]
#         if (identical(func, product3)){
#           upper[7] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
#         }
#       }

#       if (!is.null(latency.limit)){
#         if (identical(func, product1)){
#           upper[4] <- latency.limit
#         } else if (identical(func, product2)){
#           upper[c(4,8)] <- latency.limit
#         } else if (identical(func, product3)){
#           upper[c(4,8,12)] <- latency.limit
#         }
#       }

#     }
  
#     if (identical(func, product1)){
#       param_names <- c('a', 'b', 'c', 'd')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d)

#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) lower[1] <- sign * lower[1] 
#       if (!is.null(upper)) upper[1] <- sign * upper[1] 

#     } else if (identical(func, product2)){
#       param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
#                         (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h)
#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5] }
#       if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5] }
    
#     } else if (identical(func, product3)){
#       param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
#                         (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h) +
#                         (i / (((j / (j + k)) ^ (j / k)) * g / (j + k))) * (1 - exp(-(x - l) / j)) * exp(-(x - l) / k) * (x >= l)
#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5]; lower[9] <- sign * lower[9] }
#       if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5]; upper[9] <- sign * upper[9] }
#     }
#   }else{
#     fun.2.fit <- func
#     y2fit <- yfilter[ind2:length(yfilter)]
#     x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)
#     param_names <- generate_param_names(fun.2.fit)
#     N.params <- length(param_names)
#     form.2.fit <- generate_formula(fun.2.fit, param_names, y_name = 'y', x_name = 'x')
#   }

#  if (method == 'MLE'){
    
#     MLE.method <- match.arg(MLE.method)
#     iter <- MLEsettings$iter
#     metropolis.scale <- MLEsettings$metropolis.scale
#     fit.attempts <- MLEsettings$fit.attempts
#     RWm <- MLEsettings$RWm  
#     sd.est <- sqrt(mean(yfilter[1:ind2]^2))

#     output <- FIT.MLE(x=x2fit, y=y2fit, func=func, N.params=N.params, sigma=sd.est, iter=iter, metropolis.scale=metropolis.scale, 
#       logpost=log.lik.post, params.start=NULL, method=MLE.method, lower=lower, upper=upper, MLE.fun.attempts=100, fit.attempts=fit.attempts, RWm=RWm) 

#   } else if (method == 'LM'){
#      output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

#   } else if (method == 'BF.LM'){
#     output <- FIT.bf.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

#   } else if (method == 'GN'){
#     output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper)

#   } else if (method == 'port'){
#     output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm='port')

#   } else if (method == 'robust'){
#     algorithm <- if (!is.null(lower) || !is.null(upper)) 'port' else 'default'
#     # currently only 'M' working
#     algorithm <-  'port'
#     output <- FIT.NLSrobust(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm=algorithm, method='M')
#   }
    

#   if (is_product_function(func)) {
#     if (identical(func, product1)){
#       df_output <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
#       output$fits[1] <- sign * output$fits[1]
#     } else if (identical(func, product2)){  
#       df_output1 <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
#       df_output2 <- out.fun(params=output$fits[5:8], interval=interval, dp=dp, sign=sign)
#       df_output <- if (df_output1[3] < df_output2[3]) rbind('fast' = df_output1, 'slow' = df_output2) else rbind('fast' = df_output2, 'slow' = df_output1)
#       output$fits[1] <- sign * output$fits[1]; output$fits[5] <- sign * output$fits[5]
#       if (output$fits[3] > output$fits[7]){
#         output$fits <- c(output$fits[5:8], output$fits[1:4])
#         output$fits.se <- c(output$fits.se[5:8], output$fits.se[1:4])
#       }
#     } else if (identical(func, product3)){  
#       # Compute the output data frames
#       df_output1 <- out.fun(params=output$fits[1:4],  interval=interval, dp=dp, sign=sign)
#       df_output2 <- out.fun(params=output$fits[5:8],  interval=interval, dp=dp, sign=sign)
#       df_output3 <- out.fun(params=output$fits[9:12], interval=interval, dp=dp, sign=sign)

#       # Create a list of the outputs
#       output_list <- list(df_output1, df_output2, df_output3)

#       # Extract the third elements from each output and determine the order
#       order_indices <- order(sapply(output_list, function(x) x[[3]]))

#       # Reorder the output list based on the third element
#       output_ordered <- output_list[order_indices]

#       # Combine the outputs in the correct order
#       df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])

#       # Adjust the order of the fits and fits.se based on the sorted order
#       output$fits <- c(output$fits[(4*order_indices[1]-3):(4*order_indices[1])],
#                        output$fits[(4*order_indices[2]-3):(4*order_indices[2])],
#                        output$fits[(4*order_indices[3]-3):(4*order_indices[3])])

#       output$fits.se <- c(output$fits.se[(4*order_indices[1]-3):(4*order_indices[1])],
#                           output$fits.se[(4*order_indices[2]-3):(4*order_indices[2])],
#                           output$fits.se[(4*order_indices[3]-3):(4*order_indices[3])])

#       # adjust signs
#       output$fits[1] <- sign * output$fits[1]
#       output$fits[5] <- sign * output$fits[5]
#       output$fits[9] <- sign * output$fits[9]
#     }

#     traces=data.frame(x=xorig, y= sign * yorig, yfilter= sign * yfilter)
#     fits <- output$fits
#     if (identical(func, product1)){
#       fits[4] <- fits[4] + baseline
#     } else if (identical(func, product2)){
#       fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline
#     } else if (identical(func, product3)){
#       fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline; fits[12] <- fits[12] + baseline
#     }    
#     traces$yfit <- func(fits,traces$x+dx)  
#     if (identical(func, product2)){
#       traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
#       traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
#     } 
#     if (identical(func, product3)){
#       traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
#       traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
#       traces$yfit3 <- product1(fits[9:12],traces$x+dx) 
#     } 

#   }else{
#     traces=data.frame(x=xorig, y=sign * yorig, yfilter=sign * yfilter)
#     fits <- output$fits
#     traces$yfit <- func(fits,x2fit)  
#     df_output <- data.frame(fits=fits)
#     df_output <- t(df_output)

#     df_output <- as.data.frame(df_output)
#     colnames(df_output) <- param_names
#   }

#   if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)

#   if (show.output){
#     print(df_output)
#   }

#   if (return.output){
#     return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces))
#   }
# }


fit_plot <- function(traces, func=product2, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, main='', bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$x, traces$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main=main)
  
  if (filter) {
    lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$x, traces$yfit, col='indianred', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2) || identical(func, product2N)) {
    lines(traces$x, traces$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3) || identical(func, product3N)) {
    lines(traces$x, traces$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}

nFIT <- function(response, n=30, N=1, IEI=50, dt=0.1, func=product2N, method= c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'), weight_method=c('none', '~y_sqrt', '~y'),
  stimulation_time=0, baseline=0, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
  latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE),  
  MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), response_sign_method = c('smooth', 'regression', 'cumsum'), half_width_fit_limit=500, dp=3, lwd=1.2, 
  xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42) {
  
  set.seed(seed)
  
  output <- NULL
  gof <- Inf
  fit_results <- list()
  df_output <- NULL
  sign <- NULL
  traces <- NULL

  fast.constraint.method <- match.arg(fast.constraint.method)

  for (i in 1:n) {
    fit_result <- tryCatch({
      
      FITN(response=response, dt=dt, func=func, N=N, IEI=IEI, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=fast.decay.limit, 
                      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, 
                      MLE.method=MLE.method, response_sign_method=response_sign_method, dp=10, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, 
                      return.output=TRUE, show.output=FALSE, show.plot=FALSE)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(fit_result)) {
      fit_results[[i]] <- fit_result[!(names(fit_result) %in% c('sign', 'traces'))]
      
      if (first.delay.constraint){
        d1 <- fit_result$output[1, 'delay'] 
        d2 <- fit_result$output[2, 'delay'] 
      }
      
      if (identical(func, product2N)){

        if (fast.constraint){

          r_percent_start <- interval[1] * 100
          r_percent_end   <- interval[2] * 100
          r_label <- paste0('r', r_percent_start, '_', r_percent_end)

          if (fast.constraint.method == 'rise'){
            fast_ <- fit_result$output[1, r_label]
            slow_ <- fit_result$output[2, r_label]
          }else if (fast.constraint.method == 'peak'){
            fast_ <- fit_result$output[1, 'tpeak'] 
            slow_ <- fit_result$output[2, 'tpeak'] 
          }
          
          if (first.delay.constraint){
            if (!(fast_ < slow_ && d1 < d2)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }else{
            if (fast_ > slow_) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }

        }

      }else if (identical(func, product3N)){
      
        if (fast.constraint){
          r_percent_start <- interval[1] * 100
          r_percent_end   <- interval[2] * 100
          r_label <- paste0('r', r_percent_start, '_', r_percent_end)

          if (fast.constraint.method == 'rise'){
            fast_  <- fit_result$output[1, r_label]
            medium_ <- fit_result$output[2, r_label]
            slow_   <- fit_result$output[3, r_label]
          }else if (fast.constraint.method == 'peak'){
            fast_   <- fit_result$output[1, 'tpeak'] 
            medium_ <- fit_result$output[2, 'tpeak'] 
            slow_   <- fit_result$output[3, 'tpeak']
          }
          
          if (first.delay.constraint){
            d3 <- fit_result$output[3, 'delay'] 
            if (!(fast_ < medium_ && medium_ < slow_ && d1 < d2 && d1 < d3)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }else{
            if (!(fast_ < medium_ && medium_ < slow_)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }
         
        }   
      }
  
      if (fit_result$gof < gof) {
        df_output <- fit_result$output
        output <- fit_result
        gof <- fit_result$gof
        sign <- fit_result$sign
        traces <- fit_result$traces      
      }
    }
  }
  
  if (is.null(output)) {
    stop("All fit attempts failed")
  }
  
  # apply half_width function and extract its value
  df_output$half_width <- mapply(
    function(A, tau1, tau2) half_width(A, tau1, tau2, limit = half_width_fit_limit)[['half_width']],
    A = sign * df_output$A1, 
    tau1 = df_output$τrise,
    tau2 = df_output$τdecay
  )

  # reorder
  # df_output <- df_output[, c(names(df_output)[1:6], 'half_width', 'delay', 'area1')]
  # df_output <- round(df_output, digits=dp)

  cols <- names(df_output)
  area_idx <- which(grepl("^area", cols))[1]
  new_order <- append(cols[-which(cols == "half_width")], "half_width", after = area_idx - 1)
  df_output <- df_output[, new_order]
  df_output <- round(df_output, digits=dp)

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  
  if (show.output) print(df_output)
  
  if (return.output) {
    return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces, fit_results=fit_results))
  }
}

first.peak <- function(y, threshold=0.1) {
  
  sign <- sign_fun(y)
  y <- sign * y
  
  # Define a simple threshold to ignore tiny fluctuations (for example, 10% of max)
  threshold_value <- threshold * max(abs(y), na.rm=TRUE)

  # Initialize peak_index to NA, which will be returned if no peak is found
  peak_index <- NA
  
  # Loop through the data and find the first peak (positive value followed by a negative one)
  for (i in 2:(length(y) - 1)) {
    if (!is.na(y[i]) && y[i] > y[i - 1] && y[i] > y[i + 1] && y[i] > threshold_value) {
      peak_index <- i
      break  # Exit the loop once the first peak is found
    }
  }

  # Return the index of the first peak or NA if none is found
  return(peak_index)
}



last.peak <- function(y, threshold=0.9) {
  
  sign <- sign_fun(y)
  y <- sign * y

  # Define a simple threshold to ignore tiny fluctuations (for example, 10% of max)
  threshold_value <- threshold * max(abs(y), na.rm=TRUE)

  # Reverse the data to find the "first peak" in the reversed vector, which corresponds to the last peak
  y_rev <- rev(y)

  peak_index_rev <- NA
  for (i in 2:(length(y_rev) - 1)) {
    if (!is.na(y_rev[i]) && y_rev[i] > y_rev[i - 1] && y_rev[i] > y_rev[i + 1] && y_rev[i] > threshold_value) {
      peak_index_rev <- i
      break
    }
  }

  # Convert the index from the reversed vector back to the original vector's index
  if (!is.na(peak_index_rev)) {
    peak_index <- length(y) - peak_index_rev + 1
  } else {
    peak_index <- NA
  }

  return(peak_index)
}

ests_fun <- function(x, y, N=1, showplot=FALSE) {
  n <- length(y)
  
  dx=x[2]-x[1]
  k <- round(ceiling(1/dx))

  ysmooth <- rep(NA, n)
  for (i in 1:(n - k + 1)) {
    ysmooth[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
  }

  idx.peak <- which.max(abs(ysmooth))
  t.peak <- x[idx.peak]
  A <- ysmooth[idx.peak]

  # Find the index for peak value; if N=1 then indices are the same else indices for first and last peak
  idx.1st <- ifelse(N == 1 || is.na(first.peak(ysmooth)), which.max(abs(ysmooth)), first.peak(ysmooth))
  idx.Nth <- ifelse(N == 1 || is.na(last.peak(ysmooth)), idx.1st, last.peak(ysmooth))

  t.peak1 <- x[idx.1st]
  A1 <- ysmooth[idx.1st]
  A.Nth <- ysmooth[idx.Nth]
  
  # Calculate 20% and 80% of the 1st and Nth peak amplitude
  A1.20 <- 0.20 * A1; A.Nth.20 <- 0.20 * A.Nth 
  A1.80 <- 0.80 * A1; A.Nth.80 <- 0.80 * A.Nth 
  
  # Find all indices where y crosses 20% and 80% amplitude levels during rising phase
  rise_idx.20 <- which(diff(sign(ysmooth[1:idx.1st] - A1.20)) != 0)
  rise_idx.80 <-  which(diff(sign(ysmooth[1:idx.1st] - A1.80)) != 0)
  
  # Find all indices where y crosses 20% and 80% amplitude levels during decaying phase
  decay_idx.80 <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A.Nth.80)) != 0) + idx.Nth
  decay_idx.20 <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A.Nth.20)) != 0) + idx.Nth
  
  # Handle the case where decay_idx.20 does not exist
  if (length(decay_idx.20) == 0) {
    # Perform linear extrapolation
    if (length(decay_idx.80) > 1) {
      # Get the last two points near the 80% level
      x1 <- x[decay_idx.80[length(decay_idx.80) - 1]]
      x2 <- x[decay_idx.80[length(decay_idx.80)]]
      y1 <- y[decay_idx.80[length(decay_idx.80) - 1]]
      y2 <- y[decay_idx.80[length(decay_idx.80)]]
      
      # Calculate the slope of the line
      slope <- (y2 - y1) / (x2 - x1)
      
      # Extrapolate to the 20% level
      t.decay.20 <- x2 + (A1.20 - y2) / slope
    } else {
      # If we don't have enough points to extrapolate, use the last time point
      t.decay.20 <- x[length(x)]
    }
  } else {
    t.decay.20 <- x[decay_idx.20]
  }
  
  # Calculate times for these points
  t.rise.20 <- x[rise_idx.20]
  t.rise.80 <- x[rise_idx.80]
  t.decay.80 <- x[decay_idx.80]
  
  # Average the times if there are multiple points
  avg_t.rise.20 <- median(t.rise.20)
  avg_t.rise.80 <- median(t.rise.80)
  avg_t.decay.80 <- median(t.decay.80)
  avg_t.decay.20 <- median(t.decay.20)
  
  # Calculate 20-80% rise time and 80-20% decay time
  rise.time <- avg_t.rise.80 - avg_t.rise.20

  if (rise.time==0) rise.time <- 0.5 * t.peak
  decay.time <- avg_t.decay.20 - avg_t.decay.80

  # Plot the original and smoothed signals
  if (showplot) {
    plot(x, y, col='indianred', xlab='x', type='l', bty='l', las=1, main='')
    lines(x, ysmooth, col='lightgray')
    # abline(h=A1.20, col='black', lty=3)
    # abline(h=A1.80, col='black', lty=3)
    abline(h= A.Nth.20, col='black', lty=3)
    abline(h= A.Nth.80, col='black', lty=3)
  }

  return(c(A, rise.time, decay.time))
}


# Use generate_lognormal_samples to produce random starting points balsed on estimated starting values
generate_lognormal_samples <- function(means, cv=0.4, n=1) {
  # Calculate parameters for the lognormal distribution
  calculate_lognormal_params <- function(mean, cv) {
    mu <- log(mean) - 0.5 * log(1 + cv^2)
    sigma <- sqrt(log(1 + cv^2))
    list(mu = mu, sigma = sigma)
  }
  
  # Generate samples for each mean
  samples <- sapply(means, function(mean) {
    params <- calculate_lognormal_params(mean, cv)
    rlnorm(n, meanlog = params$mu, sdlog = params$sigma)
  })
  return(samples)
}

adjust_upper1 <- function(upper, ests, factor) {
  changed <- FALSE
  for (i in 1:length(upper)) {
    if (upper[i] < ests[i]) {
      upper[i] <- factor * ceiling(5 * ests[i] / factor)
      changed <- TRUE
    }
  }
  return(list(upper = upper, changed = changed))
}

bounds_check <- function(ests, model='product2', lower=NULL, upper=NULL, fast.decay.limit=NULL, latency.limit=NULL){
  if (model == 'product') {
    result <- adjust_upper1(upper, ests, 10)
    N.params <- 4
    upper <- result$upper
  } else if (model == 'product2') {
    result <- adjust_upper1(upper, ests, 10)
    N.params <- 8
    upper <- result$upper
    if (upper[1] < upper[5]) upper[1] <- upper[5] else if (upper[5] < upper[1]) upper[5] <- upper[1]
  }
  if (result$changed) {
    message("upper bounds have been changed to: [", paste(upper, collapse = ", "), "]")
  }

  if (is.null(lower)) lower <- rep(0, N.params) 
  return(list(lower=lower, upper=upper))
  }

adjust_upper <- function(upper, ests, indices, factor) {
  changed <- FALSE
  for (i in indices) {
    est_index <- ifelse(i %in% c(5, 6, 7), (i - 4), (i - 1) %% length(ests) + 1)
    if (upper[i] < ests[est_index]) {
      upper[i] <- factor * round(5 * ests[est_index] / factor)
      changed <- TRUE
    }
  }
  return(list(upper = upper, changed = changed))
}

start_fun <- function(x, y, cv=0.4, showplot=FALSE, model='product2', lower=NULL, upper=NULL, fast.decay.limit=NULL, latency.limit=NULL){

  ests <- ests_fun(x, y)

  if (!is.null(upper)) {
    if (model == 'product') {
      result <- adjust_upper(upper, ests, c(1, 2, 3), 10)
    } else if (model == 'product2') {
      result <- adjust_upper(upper, ests, c(1, 2, 3, 5, 6, 7), 10)
    }
    
    upper <- result$upper
    if (result$changed) {
      message("upper bounds have been changed to: [", paste(upper, collapse = ", "), "]")
    }

    if (model == 'product'){
      success<- FALSE
      for (i in 1:1e4) {
        st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
        if (all(st[1:3] < upper[1:3])){
          success <- TRUE
          break
        }
      }
      if (!success) {
        stop("failed to generate initial starting values (1e4 attempts)")
      }

      st <- c(st, runif(1)*upper[4])

    }  else if (model == 'product2'){
        success <- FALSE
        ests1 <- ests; ests1 <- ests1/2
        ests2 <- ests; ests2[1] <- ests2[1]/2

        for (i in 1:1e4) {
        st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
        if ((all(st[1:3] < upper[1:3])) && all(st[4:6] < upper[5:7])) {
          success <- TRUE
          break
        }
      }
      if (!success) {
        stop("failed to generate initial starting values (1e4 attempts)")
      }
      st <- c(st[1:3], runif(1)*upper[4], st[4:6], runif(1)*upper[4])
    } 

  } else {
  
    if (model == 'product'){
        success <- FALSE
        if (!is.null(fast.decay.limit)) {
          for (i in 1:1e4) {
            st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
            if (st[3] < fast.decay.limit) {
              success <- TRUE
              break
            }
          }
          if (!success) {
            stop("failed to generate initial starting values where τdecay is less than the fast decay limit (1e4 attempts)")
          }
        }else{
          st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
        }

        st <- if (!is.null(latency.limit)) c(st, runif(1)*latency.limit) else c(st, runif(1)*5)

        if (is.null(lower)) lower= c(0, 0, 0, 0) # rep(0, N.params)
        if (is.null(upper)) upper = c(ests * 10, Inf)
        if (!is.null(fast.decay.limit)) upper[3] <- fast.decay.limit
        if (!is.null(latency.limit)) upper[4] <- latency.limit

      } else if (model == 'product2'){
        success <- FALSE
        # 
        ests1 <- ests; ests1 <- ests1/2
        ests2 <- ests; ests2[1] <- ests2[1]/2
        
        if (!is.null(fast.decay.limit)) {
          for (i in 1:1e4) {
            st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
            if (st[3] < fast.decay.limit) {
              success <- TRUE
              break
            }
          }
          if (!success) {
            stop("failed to generate initial starting values where τdecay is less than the fast decay limit (1e4 attempts)")
          }
        } else {
          st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
        }

        st <- if (!is.null(latency.limit)) c(st[1:3], runif(1)*latency.limit, st[4:6], runif(1)*latency.limit) else c(st[1:3], runif(1)*5, st[4:6], runif(1)*5)

        if (is.null(lower)) lower <- c(0, 0, 0, 0, 0, 0, 0, 0) # rep(0, N.params)
        if (is.null(upper)) upper <- rep(c(ests * 10, Inf),2)
        if (!is.null(fast.decay.limit)) upper[3] <- fast.decay.limit
        if (!is.null(latency.limit)) upper[c(4,8)] <- latency.limit

      }
    }
    return(list(start=st, lower=lower, upper=upper))
  }

# analysis functions
load_data <- function(wd, name) {
    
  # Create the file path
  file_path <- file.path(wd, paste0(name, '.', 'csv'))
  
  # Load the csv data
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  return(data)
}

moving_avg <- function(y, n = 5) {
  sign <- sign_fun(y)
  y <- y * sign
  y_length <- length(y)
  result <- rep(NA, y_length)
  
  for (i in 1:y_length) {
    # Determine the start and end indices for the window
    start_idx <- max(1, i - floor(n / 2))
    end_idx <- min(y_length, i + floor(n / 2))
    
    # Calculate the mean for the current window
    result[i] <- mean(y[start_idx:end_idx], na.rm = TRUE)
  }
  
  return(max(result) * sign)
}

# load_data2 <- function(wd, name) {
#   # Create the file path
#   file_path <- file.path(wd, paste0(name, '.', 'xlsx'))
  
#   # Load the Excel file
#   workbook <- openxlsx::loadWorkbook(file_path)
  
#   # Get the sheet names
#   sheet_names <- openxlsx::getSheetNames(file_path)
  
#   # Initialize an empty list to store each sheet's data
#   data_list <- list()
  
#   # Loop through each sheet and read the data into the list
#   for (sheet in sheet_names) {
#     data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet)
#   }
  
#   return(data_list)
# }


load_data2 <- function(wd, name, header = TRUE) {
  # Create the file path
  file_path <- file.path(wd, paste0(name, '.', 'xlsx'))
  
  # Load the Excel file
  workbook <- openxlsx::loadWorkbook(file_path)
  
  # Get the sheet names
  sheet_names <- openxlsx::getSheetNames(file_path)
  
  # Initialize an empty list to store each sheet's data
  data_list <- list()
  
  # Loop through each sheet and read the data into the list
  for (sheet in sheet_names) {
    data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet, colNames = header)
  }
  
  return(data_list)
}

# save list to excel spreadsheet
list2excel <- function(data_list, file_name, wd = getwd(), center_align = TRUE) {
  # Load the openxlsx library
  library(openxlsx)
  
  # Create a new workbook
  workbook <- createWorkbook()
  
  # Define styles: bold and centered for headers, and centered for data
  header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
  center_style <- createStyle(halign = "center", valign = "center")
  
  # Loop over each element in the list
  for (i in seq_along(data_list)) {
    # Use the name of the list element as the sheet name
    sheet_name <- names(data_list)[i]
    
    # Default to "Sheet1", "Sheet2", etc., if the name is missing
    if (is.null(sheet_name) || sheet_name == "") {
      sheet_name <- paste0("Sheet", i)
    }
    
    # Add a new sheet with the specified name to the workbook
    addWorksheet(workbook, sheet_name)
    
    # Write the data to the sheet
    writeData(workbook, sheet_name, data_list[[i]])
    
    # Apply header style to make headers bold and centered
    addStyle(workbook, sheet_name, style = header_style, rows = 1, cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    
    # Apply center alignment to all cells if center_align is TRUE
    if (center_align) {
      addStyle(workbook, sheet_name, style = center_style, rows = 1:(nrow(data_list[[i]]) + 1), 
               cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    }
  }
  
  # Create the full file path
  file_path <- file.path(wd, file_name)
  
  # Save the workbook
  saveWorkbook(workbook, file_path, overwrite = TRUE)
}

# Save each data frame as individual CSV files
list2csv <- function(data_list, filename, wd) {
  if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
  
  # Extract prefix from filename (remove .xlsx or .csv extension)
  prefix <- gsub("\\.(xlsx|csv)$", "", filename)
  
  lapply(names(data_list), function(name) {
    filepath <- file.path(wd, paste0(prefix, " ", name, ".csv"))
    write.csv(data_list[[name]], filepath, row.names = FALSE)
  })
  
  invisible(NULL)
}


peak.fun <- function(y, dt, stimulation_time, baseline, smooth=5){
  
  idx1 <- (stimulation_time - baseline) / dt
  idx2 <- baseline / dt

  y1 <- y[idx1:length(y)]

  window_idx <- smooth/dt
  y1 <- y1 - mean(y1[0:idx2])
  return(moving_avg(y1, n = smooth)) 
}



raw_plot <- function(response, dt=0.1, stimulation_time=0, baseline=0, smooth=5, y_abline=0.1, height=5, width=5){

    y <- response

    if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
      y <- y[!is.na(y)]
    }
  
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

    ind1 <- (stimulation_time - baseline)/dt
    ind2 <- stimulation_time/dt
    y2plot <- y - mean(y[ind1:ind2])

    dev.new(width=width, height=height, noRStudioGD=TRUE)
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by = dt)
    plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty=1)
    abline(h = peak * y_abline, col = 'black', lwd = 1, lty=3)

    # Add a label to the abline
    text(x=max(X[ind1:length(X)]) * 0.95, y=peak * (y_abline - 0.025), labels=y_abline, pos=4)

}

# determine_tmax <- function(y, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, height=5, width=5){ 
#   if (is.null(tmax)){
#     # Plot the data
#     # plot(x, y, col='indianred', xlab='x', type='l', bty='l', las=1, main=paste('trace', ii))
#     peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

#     # Get one fit and assess
#     ind1 <- (stimulation_time - baseline)/dt
#     ind2 <- stimulation_time/dt
#     y2plot <- y - mean(y[ind1:ind2])

#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#     Y <- y2plot[ind1:length(y2plot)]
#     X <- seq(0, dt * (length(Y) - 1), by = dt)
#     plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
#     abline(h = 0, col = 'black', lwd = 1, lty=1)
#     abline(h = peak * y_abline, col = 'black', lwd = 1, lty=3)

#     # Add a label to the abline
#     text(x=max(X[ind1:length(X)]) * 0.95, y=peak * (y_abline - 0.025), labels=y_abline, pos=4)

#     # Prompt user for the range of x to use for nFIT
#     x_limit <- NA
#     while (is.na(x_limit)) {
#       cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
#       x_limit <- as.numeric(readLines(n = 1))
#       if (is.na(x_limit)) {
#         cat('\nInvalid input. Please enter a numeric value.\n')
#       }
#     }
#     dev.off()
#   }else{
#     x_limit <- tmax
#   }

#   x_limit <- x_limit + stimulation_time - baseline
#   return(x_limit)
# }

# determine_tmax <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, height=5, width=5, prompt=TRUE) { 
#   if (is.null(tmax)){
#     peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

#     ind1 <- as.integer((stimulation_time - baseline)/dt)
#     ind2 <- as.integer(stimulation_time/dt)
#     y2plot <- y - mean(y[ind1:ind2])

#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#     Y <- y2plot[ind1:length(y2plot)]
#     X <- seq(0, dt * (length(Y) - 1), by = dt)

#     out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
#     A_abline <- out[1]
#     avg_t.abline <- out[2]
#     avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline

#     plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
#     abline(h = 0, col = 'black', lwd = 1, lty=1)
    
#     # Get the left and bottom of the plot
#     left_axis <- par("usr")[1]
#     bottom_axis <- par("usr")[3]
    
#     # Horizontal dotted line from X=0 to X=avg_t.abline at height A_abline
#     lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
#     # Vertical dotted line down to the bottom of the plot
#     lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)

#     # Add a label to the abline
#     ind3 <- as.integer(avg_t.abline/dt)
#     text(x=max(X[ind1:ind3])*1.05, y=A_abline * 1.2, labels=paste0(y_abline*100, ' %'), pos=4, cex=0.6)

#     text(x=max(X[ind1:ind3])*1.05, y=bottom_axis*0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=0.6)

#     if (prompt) {
#       x_limit <- NA
#       while (is.na(x_limit)) {
#         cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
#         x_limit <- as.numeric(readLines(n = 1))
#         if (is.na(x_limit)) {
#           cat('\nInvalid input. Please enter a numeric value.\n')
#         }
#       }
#       dev.off()
#     } else {
#       x_limit <- avg_t.abline
#     }
#   } else {
#     x_limit <- tmax
#   }

#   x_limit <- x_limit + stimulation_time - baseline
#   return(x_limit)
# }


determine_tmax <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, ylab=NULL, height=5, width=5, prompt=TRUE) { 
  if (is.null(tmax)){
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
    
    ind1 <- as.integer((stimulation_time - baseline)/dt)
    ind2 <- as.integer(stimulation_time/dt)
    y2plot <- y - mean(y[ind1:ind2])
    
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by = dt)
    
    out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
    A_abline <- out[1]
    avg_t.abline <- out[2]
    avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline
    
    plot(X, Y, col='indianred', xlab='time (ms)', ylab=ylab, type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty=1)
    
    # Get the left and bottom of the plot
    left_axis <- par("usr")[1]
    bottom_axis <- par("usr")[3]
    
    # Horizontal dotted line from X=0 to X=avg_t.abline at height A_abline
    lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
    # Vertical dotted line down to the bottom of the plot for the abline
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
    # Add a label to the abline
    ind3 <- as.integer(avg_t.abline/dt)
    text(x=max(X[ind1:ind3])*1.05, y=A_abline * 1.2, labels=paste0(y_abline*100, ' %'), pos=4, cex=0.6)
    text(x=max(X[ind1:ind3])*1.05, y=bottom_axis*0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=0.6)
    
    # asterisk and label stim
    stim_index <- round(baseline/dt) + 1
    if (stim_index > length(X)) stim_index <- length(X)
    points(X[stim_index], Y[stim_index], pch=8, col='darkgray', cex=1)  # pch=8 is a star symbol.

    # Place the label "stim" on the same y level as the star, slightly to the right.
    x_offset <- 0.02 * diff(range(X))  # horizontal offset based on the range of X
    text(x = X[stim_index] + x_offset, y = Y[stim_index], labels = "stim", pos = 4, col = 'darkgray', cex = 0.6)    

    if (prompt) {
      x_limit <- NA
      while (is.na(x_limit)) {
        cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
        x_limit <- as.numeric(readLines(n = 1))
        if (is.na(x_limit)) {
          cat('\nInvalid input. Please enter a numeric value.\n')
        }
      }
      dev.off()
    } else {
      x_limit <- avg_t.abline
    }
  } else {
    x_limit <- tmax
  }
  
  x_limit <- x_limit + stimulation_time - baseline
  return(x_limit)
}


abline_fun <- function(x, y, N = 1, y_abline = 0.1) {
  sign <- sign_fun(y)
  y <- y * sign
  n <- length(y)
  
  if (n < 2) return(NA)  # Ensure input length is valid
  
  dx <- x[2] - x[1]
  k <- round(ceiling(1 / dx))

  # Ensure we have enough points to perform the moving average
  if (n < k) return(NA)
  
  # Vectorized smoothing using a moving average filter
  ysmooth <- stats::filter(y, rep(1/k, k), sides = 2) # Use two-sided filtering for better alignment
  
  # Ensure the smoothed signal is of the correct length
  ysmooth <- as.numeric(ysmooth) # Convert from 'ts' object to numeric
  ysmooth[is.na(ysmooth)] <- 0 # Replace NAs from edges with 0 (or you can opt for another value)
  
  # Identify the first and Nth peaks
  idx.1st <- ifelse(N == 1 || is.na(first.peak(ysmooth)), which.max(abs(ysmooth)), first.peak(ysmooth))
  idx.Nth <- ifelse(N == 1 || is.na(last.peak(ysmooth)), idx.1st, last.peak(ysmooth))

  # Check if peaks were found
  if (is.na(idx.1st) || is.na(idx.Nth)) return(NA)

  A.Nth <- ysmooth[idx.Nth]
  
  # Calculate threshold amplitude level based on y_abline
  A_abline <- y_abline * A.Nth 
    
  # Find indices where the smoothed signal crosses the threshold during the decaying phase
  abline_idx <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A_abline)) != 0) + idx.Nth
  
  # If no crossing points are found, return NA
  if (length(abline_idx) == 0) return(NA)

  # Find corresponding times and average if there are multiple crossings
  t.abline <- x[abline_idx]
  avg_t.abline <- median(t.abline)
  A_abline <- sign * A_abline

  return(c(A_abline, avg_t.abline))
}

sequential_fit <- function(response, n=30, dt=0.1, func=product2N, N=1, IEI=50, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), weight_method=c('none', '~y_sqrt', '~y'),
    stimulation_time=0, baseline=0, tmax=NULL, y_abline=0.1, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
    latency.limit=NULL, lower=NULL, upper=NULL,  filter=FALSE, fc=1000, interval=c(0.1, 0.9), MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
    MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), smooth=5, response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, 
    lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42){
  
  y <- response
  x_limit <- determine_tmax(y=y, N=N, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=y_abline, width=width, height=height)

  x <- seq(0, (length(y) - 1) * dt, by = dt)

  # Adjust response based on user input
  y2fit <- y[x < x_limit]
  x2fit <- seq(0, dt * (length(y2fit)-1), by = dt)

  ind1 <- (stimulation_time - baseline)/dt
  ind2 <- stimulation_time/dt
  ind3 <- length(y2fit)
  # ind4 <- length(y)
   
  # bl <- mean(y[ind1:ind2])

  # xfit <- seq(0, dt * (ind3 - ind2), by = dt)
  out <- nFIT(response=y2fit, n=n, dt=dt, N=N, IEI=IEI, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, 
    fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, 
    first.delay.constraint=first.delay.constraint, latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, 
    fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, response_sign_method=response_sign_method, 
    dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=return.output, show.output=FALSE, show.plot=FALSE, seed=seed)

  if (show.plot) fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  if (show.output) print(out$output)

  if (!fast.constraint && (!identical(func, product1))){
      
      # Define the specific message based on the fast.constraint.method
      message_part <- switch(fast.constraint.method,
                       rise = "the fastest rise",
                       peak = "the fastest time to peak")

      # Prompt user if they want to repeat with fast.rise.constraint=TRUE
      cat('\nDo you want to repeat with "fast constraint" turned on?',
          '\nThis constraint ensures the response with the fastest decay also has', message_part, "(y/n): ")
      
      repeat_with_constraint <- tolower(readLines(n = 1))
      
      if (repeat_with_constraint == 'y') {
      dev.off()
      out <- nFIT(response=y2fit, n=n, dt=dt, N=N, IEI=IEI, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, 
        fast.decay.limit=fast.decay.limit, fast.constraint=TRUE, fast.constraint.method=fast.constraint.method, 
        first.delay.constraint=first.delay.constraint,latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, 
        fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, response_sign_method=response_sign_method, 
        dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=return.output, show.output=FALSE, show.plot=FALSE, seed=seed)

      if (show.plot) fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
      if (show.output) print(out$output)
    }
  }
 
  N <- length(out$fits)
  params <- out$fits[(N-3):N]

  xfit <- seq(0, dt * (ind3 - ind2), by = dt)
  yfit <- product1(params, xfit)

  ynew <- c(y2fit[1:(ind2-1)], (y2fit[ind2:ind3] - yfit))

  traces <- data.frame(x=x2fit, y=y2fit, yreduced=ynew)
    
  # plot(x2fit,y2fit,type='l')
  # lines(x2fit, ynew, type='l')

  return(list(params=params, traces=traces, tmax=x_limit))

}


nFIT_sequential <- function(response, n=30, N=1, IEI=50, dt=0.1, func=product2N, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
    weight_method = c('none', '~y_sqrt', '~y'), stimulation_time=0, baseline=0, fit.limits=NULL, fast.decay.limit=NULL, fast.constraint=FALSE, 
    fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
    MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), 
    response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, 
    show.output=TRUE, show.plot=TRUE, seed=42){

  fast.constraint.method <- match.arg(fast.constraint.method)

  y <- response
  if (identical(func, product1N)){
    functions <- list(product1N)
  }else if (identical(func, product2N)){
    functions <- list(product2N, product1N)
  }else if (identical(func, product3N)){
    functions <- list(product3N, product2N, product1N)
  }

  if (!is.null(fit.limits) && length(fit.limits) != length(functions)){
    warning("fit.limits must contain same number of elements as fitted product equations; fit.limits will be ignored")
    fit.limits=NULL
  }
  outputs <- vector("list", length(functions))

  for (ii in 1:length(functions)) {
    if (ii == 1) {
        y2fit <- y
    } else {
        y2fit <- outputs[[ii - 1]]$traces$yreduced
        if (!is.null(upper)) upper <- upper[1:(length(upper) - 4)]
        if (!is.null(lower)) lower <- lower[1:(length(lower) - 4)]
    }
    
    # fit.limits is used to reproduce the fits rapidly; do not show plots - stimulation_time + baseline
    tmax <- if (is.null(fit.limits)) NULL else fit.limits[ii] #  + stimulation_time - baseline
    output_logic <- if (is.null(tmax)) TRUE else FALSE 
    outputs[[ii]] <- sequential_fit(response=y2fit, n=n, N=N, IEI=IEI, dt=dt, func=functions[[ii]], method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
                    fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
                    latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, 
                    response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=output_logic, show.plot=output_logic, seed=seed)

    if (is.null(tmax)) { # Prompt user if they want to repeat
      cat('Do you want to repeat fit with new time base? (y/n): ')
      repeat_fit <- tolower(readLines(n = 1))
      while (repeat_fit == 'y') {
          dev.off()
          outputs[[ii]] <- sequential_fit(response=y2fit, n=n, N=N, IEI=IEI, dt=dt, func=functions[[ii]], method=method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
                        fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
                        latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, 
                        response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=output_logic, show.plot=output_logic, seed=seed)

          cat('Do you want to repeat fit with new time base? (y/n): ')
        repeat_fit <- tolower(readLines(n = 1))
        }
        dev.off()
    }
  }
  
  fits_list <- lapply(outputs, function(out) out$params)
  # Sort the list by the increasing decay (3rd element)
  fits_list_sorted <- fits_list[order(sapply(fits_list, function(x) x[3]))]
  fits <- unlist(fits_list_sorted)

  t_limits <- sapply(outputs, function(out) out$tmax)

  # if (identical(func, product1N)){
  #   df_output <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
  # } else if (identical(func, product2)){  
  #   df_output1 <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
  #   df_output2 <- out.fun(params=fits[5:8], interval=interval, dp=dp, sign=1)
  #   df_output <- rbind('fast' = df_output1, 'slow' = df_output2)
  # } else if (identical(func, product3N)){  
  #   # Compute the output data frames
  #   df_output1 <- out.fun(params=fits[1:4],  interval=interval, dp=dp, sign=1)
  #   df_output2 <- out.fun(params=fits[5:8],  interval=interval, dp=dp, sign=1)
  #   df_output3 <- out.fun(params=fits[9:12], interval=interval, dp=dp, sign=1)
  #   df_output <- rbind('fast' = df_output1, 'medium' =  df_output2, 'slow' =  df_output3)
  # }
  
  if (identical(func, product1N)){
    
    df_output <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
  
  } else if (identical(func, product2N)){  
    df_output1 <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=1)
      # Create a list of the outputs
    output_list <- list(df_output1, df_output2)
    # Extract the third elements from each output and determine the order
    order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
    # Reorder the output list based on the third element
    output_ordered <- output_list[order_indices]

    # Combine the outputs in the correct order
    df_output <- rbind('fast' = output_ordered[[1]], 'slow' = output_ordered[[2]])



  } else if (identical(func, product3N)){ 

    df_output1 <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=1)
    df_output3 <- out.fun(params=fits[(2*N+7):(3*N+9)], interval=interval, dp=dp, sign=1)

    # Create a list of the outputs
    output_list <- list(df_output1, df_output2, df_output3)
    # Extract the third elements from each output and determine the order
    order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
    # Reorder the output list based on the third element
    output_ordered <- output_list[order_indices]

    # Combine the outputs in the correct order
    df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])
  }

  traces <- traces_fun2(y=y, fits=fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  if (show.output) print(df_output)
    
  idx1 <- baseline/dt
  idx2 <- max(t_limits)/dt

  k = length(fits)
  gof.se <- sqrt(sum((traces$y[idx1:idx2] - traces$yfit[idx1:idx2])^2) / (length(traces$y[idx1:idx2])-k))
  msc <- model.selection.criteria(coeffs=fits, x=traces$x[idx1:idx2]-traces$x[idx1], y=traces$y[idx1:idx2], func=func, N=N, IEI=IEI)

  if (return.output) {
    out <- list(output=df_output, fits=fits, gof=gof.se, AIC=msc[1], BIC=msc[2], traces=traces, fit.limits=t_limits - stimulation_time + baseline)
    return(out)
  }
}

analyse_PSC <- function(response, dt=0.1, n=30, N=1, IEI=50, stimulation_time=150, baseline=50, smooth=5, func=product2N,  method=c("BF.LM", "LM", "GN", "port", "robust", "MLE"), 
  weight_method=c('none', '~y_sqrt', '~y'), sequential.fit=FALSE, fit.limits=NULL, MLEsettings=list(iter=1e3, metropolis.scale=1.5, fit.attempts=10, RWm=FALSE), 
  filter=FALSE, fc=1000, interval=c(0.1, 0.9), lower=NULL, upper=NULL,  fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), 
  first.delay.constraint=FALSE, latency.limit=NULL, rel.decay.fit.limit=0.1, half_width_fit_limit=500, dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', 
  downsample=1, return.output=TRUE, height=5, width=5, seed=42) {
  
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  y <- y[seq(1, length(y), by = downsample)]
  dt <- dt * downsample

  x <- seq(0, (length(y) - 1) * dt, by = dt)

  if (!sequential.fit){
    
    tmax <- fit.limits
    x_limit <- determine_tmax(y=y, N=N, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=rel.decay.fit.limit, ylab=ylab, width=width, height=height)   
   
    adjusted_response <- y[x < x_limit]
    
    # Execute nFIT
    out <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func, method=method, weight_method=weight_method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
      filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method,
      first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, half_width_fit_limit=half_width_fit_limit, dp=dp, height=height, width=width, seed=seed)

    out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)
    fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
    
    if (!fast.constraint && (!identical(func, product1)) && is.null(fit.limits)){
      
      fast.constraint.method <- match.arg(fast.constraint.method)
      # Define the specific message based on the fast.constraint.method
      message_part <- switch(fast.constraint.method,
                       rise = "the fastest rise",
                       peak = "the fastest time to peak")

      # Prompt user if they want to repeat with fast.rise.constraint=TRUE
      cat('\nDo you want to repeat with "fast constraint" turned on?',
          '\nThis constraint ensures the response with the fastest decay also has', message_part, "(y/n): ")

      repeat_with_constraint <- tolower(readLines(n = 1))
      
      if (repeat_with_constraint == 'y') {
        dev.off()
        out <- nFIT(response=adjusted_response, n=n,  N=N, IEI=IEI,dt=dt, func=func, method=method, weight_method=weight_method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
          filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=TRUE, fast.constraint.method=fast.constraint.method, half_width_fit_limit=half_width_fit_limit,
          first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, dp=dp, height=height, width=width, seed=seed)

        out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)
        fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
    
      }
    }
    
  }else{
    
    out <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, fit.limits=fit.limits, 
      fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method,  half_width_fit_limit=half_width_fit_limit, 
      dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=TRUE, show.plot=TRUE, seed=seed)
  }

  if (return.output) return(out)

}

traces_fun <- function(y, fits, dt=0.1,  stimulation_time=150, baseline=50, func=product2, filter=FALSE, fc=1000, height=5, width=5){
  
  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)


  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces=data.frame(x=xorig, y=yorig, yfilter=yfilter)
  if (identical(func, product1)){
    fits[4] <- fits[4] + baseline
  } else if (identical(func, product2)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline
  } else if (identical(func, product3)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline; fits[12] <- fits[12] + baseline
  }    
  traces$yfit <- func(fits,traces$x+dx)  
  if (identical(func, product2)){
    traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
    traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
  } 
  if (identical(func, product3)){
    traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
    traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
    traces$yfit3 <- product1(fits[9:12],traces$x+dx) 
  } 

  return(traces)
}

traces_fun2 <- function(y, fits, dt=0.1,  N=1, IEI=50, stimulation_time=150, baseline=50, func=product2N, filter=FALSE, fc=1000){
  
  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces <- data.frame(x = xorig, y = yorig, yfilter = yfilter)
  
  # Define the list of functions
  func_list <- list(product1N, product2N, product3N)

  # Find the index of the matching function
  func_index <- sapply(func_list, function(f) identical(f, func))

  # Check if func matches one of the productN functions and apply the corresponding baseline adjustments
  if (any(func_index)) {
    index <- which(func_index)  # Get the index of the matching function
    for (i in 1:index) {
      fits[i * N + (3 * i)] <- fits[i * N + (3 * i)] + baseline
    }
  }

  traces$yfit <- func(params=fits, x=traces$x+dx, N=N, IEI=IEI) 

  if (identical(func, product2N)){
    traces$yfit1 <- product1N(params=fits[1:(N+3)], x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)], x=traces$x+dx, N=N, IEI=IEI) 
  } 
  if (identical(func, product3N)){
    traces$yfit1 <- product1N(params=fits[1:(N+3)],  x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)],  x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit3 <- product1N(params=fits[(2*N+7):(3*N+9)], x=traces$x+dx, N=N, IEI=IEI) 
  } 

  return(traces)
}

# Loop through each dataset
datasets2list <- function(datasets, id='_fits'){
  out <- list()
  for (dataset in datasets) {
    # Construct the name of the fits object
    name <- paste0(dataset, id)
    
    # Check if the fits object exists
    if (exists(name)) {
      # Add the fits object to the list
      out[[dataset]] <- get(name)
    } else {
      cat(paste("Object", name, "does not exist.\n"))
    }
  }
  return(out)
}


# function to extract from list of summaries
extract_variable <- function(data_list, id = 'area', rename_datasets = TRUE) {
  # Initialize empty vectors to store the combined data
  dataset_names <- c()
  values1 <- c()
  values2 <- c()
  subjects <- c()
  
  # Create a mapping for dataset renaming based on the order in the data_list
  dataset_mapping <- setNames(seq_along(names(data_list)), names(data_list))
  
  # Loop through each item in data_list
  for (name in names(data_list)) {
    # Extract the dataframe
    df <- data_list[[name]]
    
    # Extract the columns with the specified id
    variable_columns <- df[, grep(paste0("^", id, "$"), colnames(df))]
    
    # Ensure we have exactly two columns with the specified id
    if (ncol(variable_columns) == 2) {
      # Determine the dataset name
      dataset_name <- if (rename_datasets) {
        dataset_mapping[[name]]
      } else {
        name
      }
      
      # Append data to vectors
      dataset_names <- c(dataset_names, rep(dataset_name, nrow(df)))
      values1 <- c(values1, variable_columns[, 1])
      values2 <- c(values2, variable_columns[, 2])
    } else {
      cat(paste("Dataframe", name, "does not contain exactly two '", id, "' columns.\n"))
    }
  }

  # Create the combined dataframe
  combined_df <- data.frame(
    dataset = dataset_names,
    s = 1:length(dataset_names),
    value1 = values1,
    value2 = values2
  )
  
  # Rename the value1 and value2 columns to id1 and id2
  colnames(combined_df)[3:4] <- paste0(id, 1:2)

  return(combined_df)
}


boxplot_calculator_ver0 <- function(data, type = 6, na.rm=FALSE) {
  unique_x <- unique(data$x)
  result <- data.frame(x = numeric(), Q1 = numeric(), Q3 = numeric(), Median = numeric(), Min = numeric(), Max = numeric(), MAD = numeric())
  
  for (i in 1:length(unique_x)) {
    current_x <- unique_x[i]
    d <- data$y[data$x == current_x]
    
    q1 <- quantile(d, probs = 0.25, type = type, na.rm = na.rm)
    q3 <- quantile(d, probs = 0.75, type = type, na.rm = na.rm)
    iqr <- q3 - q1  # Calculate IQR
    
    lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
    upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
    
    # Exclude outliers
    d_filtered <- d[d >= lower_bound & d <= upper_bound]
    
    median_val <- median(d, na.rm = na.rm)
    min_val <- min(d_filtered, na.rm = na.rm)
    max_val <- max(d_filtered, na.rm = na.rm)
    
    # Calculate MAD
    mad <- median(abs(d - median_val), na.rm = na.rm)
    
    result <- rbind(result, data.frame(x = current_x, Q1 = q1, Q3 = q3, Median = median_val, Min = min_val, Max = max_val, MAD = mad))
  }
  
  rownames(result) <- NULL  # Remove row names
  return(result)
}

boxplot_calculator <- function(data, type = 6, na.rm=FALSE) {
  unique_x <- unique(data$x)
  result <- data.frame(x = numeric(), Q1 = numeric(), Q3 = numeric(), Median = numeric(), Min = numeric(), Max = numeric(), MAD = numeric())
  
  for (i in 1:length(unique_x)) {
    current_x <- unique_x[i]
    d <- data$y[data$x == current_x]
    
    q1 <- quantile(d, probs = 0.25, type = type, na.rm = na.rm)
    q3 <- quantile(d, probs = 0.75, type = type, na.rm = na.rm)
    iqr <- q3 - q1  # Calculate IQR
    
    lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
    upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
    
    # Exclude outliers
    d_filtered <- d[d >= lower_bound & d <= upper_bound]
    
    median_val <- median(d, na.rm = na.rm)
    min_val <- min(d_filtered, na.rm = na.rm)
    max_val <- max(d_filtered, na.rm = na.rm)

    # added to prevewnt 'internal' whiskers
    min_val <- min(min_val, q1)
    max_val <- max(max_val, q3)
    
    # Calculate MAD
    mad <- median(abs(d - median_val), na.rm = na.rm)
    
    result <- rbind(result, data.frame(x = current_x, Q1 = q1, Q3 = q3, Median = median_val, Min = min_val, Max = max_val, MAD = mad))
  }
  
  rownames(result) <- NULL  # Remove row names
  return(result)
}
WBplot_ver0 <- function(data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSP amplitude (mV)', 
                   xrange = c(0.75, 2.25), yrange = c(0, 400), main = '', tick_length = 0.02, 
                   x_tick_interval = NULL, y_tick_interval = 100, lwd = 0.8, type = 6, na.rm=FALSE) {
  
  boxplot_values <- boxplot_calculator_ver0(data=data, type=type, na.rm=na.rm)
  
  if (is.null(x_tick_interval)) {
    x_ticks <- unique(data$x)
  } else {
    x_ticks <- seq(xrange[1], xrange[2], by = x_tick_interval)
  }
  xrange <- xrange + c(-wid, wid)
  
  # Ensure background is off and plot area is clear
  par(bg = NA)
  plot(1, type = 'n', ylim = yrange, xlim = xrange, xlab = xlab, ylab = ylab, 
       main = main, xaxt = 'n', yaxt = 'n', bty = 'n', lwd = lwd)
  
  for (i in 1:nrow(boxplot_values)) {
    # Convert current_x to numeric for arithmetic operations
    current_x <- as.numeric(boxplot_values$x[i])
    
    rect(current_x - wid, boxplot_values$Q1[i], current_x + wid, boxplot_values$Q3[i], col = 'white', lwd = lwd)
    segments(current_x, boxplot_values$Q1[i], current_x, boxplot_values$Min[i], lwd = lwd)
    segments(current_x, boxplot_values$Q3[i], current_x, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Min[i], current_x + cap, boxplot_values$Min[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Max[i], current_x + cap, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - wid * 1.1, boxplot_values$Median[i], current_x + wid * 1.1, boxplot_values$Median[i], col = 'black', lwd = 3 * lwd)
  }
  
  # Set the x-axis ticks
  axis(1, at = x_ticks, labels = x_ticks, tcl = -tick_length, lwd = lwd)
  
  # Set the y-axis ticks
  y_ticks <- seq(yrange[1], yrange[2], by = y_tick_interval)
  axis(2, at = y_ticks, tcl = -tick_length, las = 1, lwd = lwd)
}

WBplot <- function(data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSP amplitude (mV)', 
                   xrange = c(0.75, 2.25), yrange = c(0, 400), main = '', tick_length = 0.02, 
                   x_tick_interval = NULL, y_tick_interval = 100, lwd = 0.8, type = 6, log_y=FALSE, na.rm=FALSE) {
  
  boxplot_values <- boxplot_calculator(data=data, type=type, na.rm=na.rm)
  if (log_y){
    num_cols <- sapply(boxplot_values, is.numeric)
    boxplot_values[num_cols] <- lapply(boxplot_values[num_cols], log10)
  }

  if (is.null(x_tick_interval)) {
    x_ticks <- unique(data$x)
  } else {
    x_ticks <- seq(xrange[1], xrange[2], by = x_tick_interval)
  }
  xrange <- xrange + c(-wid, wid)
  
  # Ensure background is off and plot area is clear
  par(bg = NA)
  plot(1, type = 'n', ylim = yrange, xlim = xrange, xlab = xlab, ylab = ylab, 
       main = main, xaxt = 'n', yaxt = 'n', bty = 'n', lwd = lwd)
  
  for (i in 1:nrow(boxplot_values)) {
    # Convert current_x to numeric for arithmetic operations
    current_x <- as.numeric(boxplot_values$x[i])
    
    rect(current_x - wid, boxplot_values$Q1[i], current_x + wid, boxplot_values$Q3[i], col = 'white', lwd = lwd)
    segments(current_x, boxplot_values$Q1[i], current_x, boxplot_values$Min[i], lwd = lwd)
    segments(current_x, boxplot_values$Q3[i], current_x, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Min[i], current_x + cap, boxplot_values$Min[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Max[i], current_x + cap, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - wid * 1.1, boxplot_values$Median[i], current_x + wid * 1.1, boxplot_values$Median[i], col = 'black', lwd = 3 * lwd)
  }
  
  # Set the x-axis ticks
  axis(1, at = x_ticks, labels = x_ticks, tcl = -tick_length, lwd = lwd)
  
  # Set the y-axis ticks
  y_ticks <- seq(yrange[1], yrange[2], by = y_tick_interval)
  axis(2, at = y_ticks, tcl = -tick_length, las = 1, lwd = lwd)
}


# BoxPlot4 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
#                      xrange = NULL, yrange = c(-400, 0), tick_length = 0.2, x_tick_interval = NULL, y_tick_interval = 100,
#                      xlabel_angle = NULL, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5, filename = 'boxplot.svg', 
#                      height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, na.rm = FALSE, save = FALSE) {

#   response   <- as.character(formula[[2]])
#   predictors <- all.vars(formula[[3]])
  
#   if (!all(c(response, predictors) %in% names(data))) stop('Response or predictor not found in data.')
  
#   if (length(predictors) == 1) {
#     data$x <- factor(data[[predictors]])
#   } else {
#     data$x <- interaction(data[[predictors[1]]], data[[predictors[2]]], sep=' : ')
#   }
#   data$y <- data[[response]]
#   data_orig <- data

#   if (log_y) {
#     data$y[data$y <= 0] <- NA
#     if (any(yrange <= 0, na.rm=TRUE)) stop('yrange must be > 0 when log_y=TRUE')
#     data$y <- log10(data$y)
#     yrange <- log10(yrange)
#   }

#   if (is.null(xrange)) {
#     xrange <- range(as.numeric(data$x), na.rm=TRUE) + c(-wid, wid)
#   }

#   if (save) {
#     svg(filename, width=width, height=height, bg=bg)
#     on.exit(dev.off(), add=TRUE)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   orig_axis <- graphics::axis
#   if (log_y) {
#     assign('axis', function(...) {}, envir=.GlobalEnv)
#   } else {
#     assign('axis', function(side, at, labels=TRUE, tcl=NA, ...){
#       if (side==1) orig_axis(side, at=at, labels=FALSE, tcl=-tick_length, ...)
#       else         orig_axis(side, at=at, labels=labels, tcl=-tick_length, ...)
#     }, envir=.GlobalEnv)
#   }
#   on.exit(assign('axis', orig_axis, envir=.GlobalEnv), add=TRUE)

#   WBplot(data=data_orig, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main,
#          xrange=xrange, yrange=yrange, tick_length=tick_length,
#          x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval,
#          lwd=lwd, type=type, log_y=log_y, na.rm=na.rm)

#   assign('axis', orig_axis, envir=.GlobalEnv)

#   labs <- sub('\\s*:\\s*.*$', '', levels(data$x))
#   at <- seq_along(labs)
#   if (is.null(xlabel_angle)) {
#     axis(1, at=at, labels=labs, tcl=-tick_length, lwd=lwd)
#   } else {
#     axis(1, at=at, labels=FALSE, tcl=-tick_length, lwd=lwd)
#     usr <- par('usr'); y0 <- usr[3] - 0.1*diff(usr[3:4])
#     text(x=at, y=y0, labels=labs, srt=xlabel_angle, adj=1, xpd=TRUE)
#   }

#   set.seed(42)
#   data$x_jitter <- jitter(as.numeric(data$x), amount=amount)
#   if ('col' %in% names(data)) {
#     points(data$x_jitter, data$y, pch=19, col=data$col, cex=p.cex)
#   } else {
#     if ('s' %in% names(data)) {
#       counts      <- ave(!is.na(data$y), data$s, FUN=sum)
#       data$paired <- counts >= 2
#       for (subj in unique(data$s[data$paired])) {
#         sd <- subset(data, s==subj & !is.na(y))
#         sd <- sd[order(as.numeric(sd$x)), ]
#         lines(sd$x_jitter, sd$y, col='darkgray', lwd=lwd, lty=3)
#       }
#     } else {
#       data$paired <- FALSE
#     }
#     gray_col <- adjustcolor('#A9A9A9', alpha.f=alpha)
#     red_col  <- adjustcolor('#CD5C5C', alpha.f=alpha)
#     pts <- data$paired | !('s' %in% names(data))
#     points(data$x_jitter[pts], data$y[pts], pch=19, col=gray_col, cex=p.cex, lwd=lwd/3)
#     if ('s' %in% names(data)) {
#       up <- !data$paired & !is.na(data$y)
#       points(data$x_jitter[up], data$y[up], pch=19, col=red_col, cex=p.cex, lwd=lwd/3)
#     }
#   }

#   if (log_y) {
#     major_tick_len <- -0.2
#     minor_tick_len <- -0.1
#     d0 <- floor(yrange[1]):ceiling(yrange[2])
#     axis(2, at=d0, labels=10^d0, tcl=major_tick_len, las=1)
#     minor_vals <- unlist(lapply(d0, function(d) (2:9)*10^d))
#     mp <- log10(minor_vals)
#     mp <- mp[mp >= yrange[1] & mp <= yrange[2]]
#     axis(2, at=mp, labels=NA, tcl=minor_tick_len)
#   }

#   if (save) dev.off()
# }

BoxPlot4 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), tick_length = 0.2, x_tick_interval = NULL, y_tick_interval = 100,
                     xlabel_angle = NULL, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5, filename = 'boxplot.svg', 
                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, na.rm = FALSE, save = FALSE) {

  response   <- as.character(formula[[2]])
  predictors <- all.vars(formula[[3]])
  
  if (!all(c(response, predictors) %in% names(data))) stop('Response or predictor not found in data.')
  
  if (length(predictors) == 1) {
    data$x <- factor(data[[predictors]])
  } else {
    data$x <- interaction(data[[predictors[1]]], data[[predictors[2]]], sep=' : ')
  }
  data$y <- data[[response]]
  data_orig <- data

  if (log_y) {
    data$y[data$y <= 0] <- NA
    if (any(yrange <= 0, na.rm=TRUE)) stop('yrange must be > 0 when log_y=TRUE')
    data$y <- log10(data$y)
    yrange <- log10(yrange)
  }

  if (is.null(xrange)) {
    xrange <- range(as.numeric(data$x), na.rm=TRUE) + c(-wid, wid)
  }

  if (save) {
    svg(filename, width=width, height=height, bg=bg)
    on.exit(dev.off(), add=TRUE)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  orig_axis <- graphics::axis
  if (log_y) {
    assign('axis', function(...) {}, envir=.GlobalEnv)
  } else {
    assign('axis', function(side, at, labels=TRUE, tcl=NA, ...){
      if (side==1) orig_axis(side, at=at, labels=FALSE, tcl=-tick_length, ...)
      else         orig_axis(side, at=at, labels=labels, tcl=-tick_length, ...)
    }, envir=.GlobalEnv)
  }
  on.exit(assign('axis', orig_axis, envir=.GlobalEnv), add=TRUE)

  WBplot(data=data_orig, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main,
         xrange=xrange, yrange=yrange, tick_length=tick_length,
         x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval,
         lwd=lwd, type=type, log_y=log_y, na.rm=na.rm)

  assign('axis', orig_axis, envir=.GlobalEnv)

  labs <- sub('\\s*:\\s*.*$', '', levels(data$x))
  at <- seq_along(labs)
  if (is.null(xlabel_angle)) {
    axis(1, at=at, labels=labs, tcl=-tick_length, lwd=lwd)
  } else {
    axis(1, at=at, labels=FALSE, tcl=-tick_length, lwd=lwd)
    usr <- par('usr'); y0 <- usr[3] - 0.1*diff(usr[3:4])
    text(x=at, y=y0, labels=labs, srt=xlabel_angle, adj=1, xpd=TRUE)
  }

  set.seed(42)
  data$x_jitter <- jitter(as.numeric(data$x), amount = amount)

  if ("s" %in% names(data)) {
    counts <- ave(!is.na(data$y), data$s, FUN = sum)
    data$paired <- counts >= 2

    for (subj in unique(data$s[data$paired])) {
      sd <- subset(data, s == subj & !is.na(y))
      sd <- sd[order(as.numeric(sd$x)), ]

      # line color: if col exists, use the first color for that subject (e.g., slice color)
      line_col <- adjustcolor("darkgray", alpha.f = alpha)

      lines(sd$x_jitter, sd$y, col = line_col, lwd = lwd, lty = 3)
    }
  } else {
    data$paired <- FALSE
  }

  if ("col" %in% names(data)) {
    # colored points for all observations
    points(data$x_jitter, data$y, pch = 19, col = data$col, cex = p.cex)
  } else {
    # original gray/red behavior when no col exists
    gray_col <- adjustcolor("#A9A9A9", alpha.f = alpha)
    red_col  <- adjustcolor("#CD5C5C", alpha.f = alpha)

    pts <- data$paired | !("s" %in% names(data))
    points(data$x_jitter[pts], data$y[pts], pch = 19, col = gray_col, cex = p.cex, lwd = lwd/3)

    if ("s" %in% names(data)) {
      up <- !data$paired & !is.na(data$y)
      points(data$x_jitter[up], data$y[up], pch = 19, col = red_col, cex = p.cex, lwd = lwd/3)
    }
  }

  if (log_y) {
    major_tick_len <- -0.2
    minor_tick_len <- -0.1
    d0 <- floor(yrange[1]):ceiling(yrange[2])
    axis(2, at=d0, labels=10^d0, tcl=major_tick_len, las=1)
    minor_vals <- unlist(lapply(d0, function(d) (2:9)*10^d))
    mp <- log10(minor_vals)
    mp <- mp[mp >= yrange[1] & mp <= yrange[2]]
    axis(2, at=mp, labels=NA, tcl=minor_tick_len)
  }

  if (save) dev.off()
}

BoxPlot <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
  xrange = NULL, yrange = c(-400, 0), xlabel_angle = NULL, tick_length = 0.2, x_tick_interval = NULL, y_tick_interval = 100,
  lwd = 1, type = 6, amount = 0.05, p.cex = 0.5, height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, 
  na_rm_subjects = FALSE, test_result, color_palette = c('viridis', 'jet', 'darkgray'), alpha_level = 0.05, sig_offset = NULL, reverse = FALSE) {

  color_palette <- match.arg(color_palette)

  response_name <- deparse(formula[[2]])
  predictors    <- all.vars(formula[[3]])
  f_str         <- paste(deparse(formula), collapse = '')
  has_re        <- grepl('\\(1\\s*\\|', f_str)
  has_err       <- grepl('Error\\(', f_str)

  if (has_re) {
    group_var <- gsub('\\s+', '', sub('.*\\(1\\s*\\|\\s*([^\\)]+)\\).*', '\\1', f_str))
  }
  if (has_err) {
    subject_var <- gsub('\\s+', '', sub('.*Error\\(([^\\)]+)\\).*', '\\1', f_str))
  }

  fixed_str <- gsub('\\+?\\s*\\(1\\s*\\|[^\\)]+\\)', '', f_str)
  fixed_str <- gsub('\\+?\\s*Error\\([^\\)]+\\)', '', fixed_str)
  main_fml  <- as.formula(fixed_str)

  if (has_err && na_rm_subjects) {
    data <- data[ !ave(is.na(data[[response_name]]), data[[subject_var]], FUN = any), ]
  }

  keep <- c(response_name, predictors)
  if (has_err) keep <- c(keep, subject_var)
  if (has_re)  keep <- c(keep, group_var)
  data_sub <- data[, unique(keep), drop = FALSE]

  if (has_err) {
    data_sub$s <- data[[subject_var]]
  } else {
    data_sub$s <- NULL
  }

  if (has_re) {
    levs <- unique(data_sub[[group_var]])
    n   <- length(levs)
    pal <- switch(color_palette,
      viridis = {
        stops <- c('#440154', '#3B528B', '#21918C', '#5DC963', '#FDE725')
        p <- colorRampPalette(stops)(n)
        adjustcolor(p, alpha.f = alpha)
      },
      jet = {
        stops <- c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow','#FF7F00','red','#7F0000')
        p <- colorRampPalette(stops)(n)
        adjustcolor(p, alpha.f = alpha)
      },
      darkgray = {
        rep(adjustcolor('darkgray', alpha.f = alpha), n)
      }
    )
    if (reverse) pal <- rev(pal)
    names(pal)   <- levs
    data_sub$col <- pal[as.character(data_sub[[group_var]])]
    
  } else {
    data_sub$col <- NULL
  }

  BoxPlot4(formula=main_fml, data=data_sub, wid=wid, cap=cap, xlab=xlab, ylab=ylab, xlabel_angle=xlabel_angle, 
    main=main, xrange=xrange, yrange=yrange, tick_length=tick_length, y_tick_interval=y_tick_interval, 
    lwd=lwd, type=type, amount=amount, p.cex=p.cex, height=height, width=width, log_y=log_y, bg=bg, na.rm=TRUE)

  if (missing(test_result) || is.null(test_result) || nrow(test_result) == 0) return()

  data_sub$y <- data_sub[[response_name]]
  if (length(predictors) == 1) {
    data_sub$x <- factor(data_sub[[predictors[1]]])
    single_factor <- TRUE
  } else {
    data_sub$x <- interaction(data_sub[[predictors[1]]], data_sub[[predictors[2]]], sep = ' : ')
    single_factor <- FALSE
  }
  x_labels <- levels(data_sub$x)
  x_positions <- setNames(seq_along(x_labels), x_labels)
  y_span <- diff(range(data_sub$y, na.rm = TRUE))
  offset <- if (is.null(sig_offset)) 0.05 * y_span else sig_offset
  tick <- 0.25 * offset

  for (i in seq_len(nrow(test_result))) {
    p_adj <- test_result$`p adjusted`[i]
    if (is.na(p_adj) || p_adj >= alpha_level) next
    parts <- strsplit(as.character(test_result$contrast[i]), ' vs ')[[1]]
    if (length(parts) != 2) next
    lev1 <- trimws(parts[1]); lev2 <- trimws(parts[2])
    if (single_factor) {
      label1 <- lev1; label2 <- lev2; paired <- FALSE
    } else {
      comp      <- as.character(test_result$comparison[i])
      outer_lev <- sub('.*? ([^ ]+) \\(.*', '\\1', comp)
      if (grepl('\\(paired\\)', comp)) {
        label1 <- paste0(lev1, ' : ', outer_lev); label2 <- paste0(lev2, ' : ', outer_lev); paired <- TRUE
      } else {
        label1 <- paste0(outer_lev, ' : ', lev1); label2 <- paste0(outer_lev, ' : ', lev2); paired <- FALSE
      }
    }
    if (!(label1 %in% x_labels) || !(label2 %in% x_labels)) next
    x1 <- x_positions[label1]; x2 <- x_positions[label2]
    if (log_y) {
      yvals <- c(data_sub$y[data_sub$x == label1], data_sub$y[data_sub$x == label2])
      yvals <- yvals[!is.na(yvals) & yvals > 0]; if (!length(yvals)) next
      yvals <- log10(yvals); y_max <- max(yvals); y_min <- min(yvals)
      y_span_log <- diff(range(log10(data_sub[[response_name]][data_sub[[response_name]] > 0]), na.rm = TRUE))
      offset <- if (is.null(sig_offset)) 0.05 * y_span_log else sig_offset
      tick <- 0.25 * offset; shift_amt <- if (!single_factor && !paired) 6 * tick else 0
      if (y_max > 0) {
        y_line <- y_max + offset + shift_amt
        segments(x1, y_line, x1, y_line - tick, lwd = lwd)
        segments(x2, y_line, x2, y_line - tick, lwd = lwd)
        text_y <- y_line + tick
      } else {
        y_line <- y_min - offset - shift_amt
        segments(x1, y_line, x1, y_line + tick, lwd = lwd)
        segments(x2, y_line, x2, y_line + tick, lwd = lwd)
        text_y <- y_line - tick
      }
    } else {
      yvals <- c(data_sub$y[data_sub$x == label1], data_sub$y[data_sub$x == label2])
      yvals <- yvals[!is.na(yvals)]; if (!length(yvals)) next
      y_max <- max(yvals); y_min <- min(yvals); shift_amt <- if (!single_factor && !paired) 6 * tick else 0
      if (y_max > 0) {
        y_line <- y_max + offset + shift_amt
        segments(x1, y_line, x1, y_line - tick, lwd = lwd)
        segments(x2, y_line, x2, y_line - tick, lwd = lwd)
        text_y <- y_line + tick
      } else {
        y_line <- y_min - offset - shift_amt
        segments(x1, y_line, x1, y_line + tick, lwd = lwd)
        segments(x2, y_line, x2, y_line + tick, lwd = lwd)
        text_y <- y_line - tick
      }
    }
    segments(x1, y_line, x2, y_line, lwd = lwd)
    text((x1 + x2)/2, text_y, labels = '*', cex = 1.2)
  }
}


BoxPlot_ver0 <- function(data, wid=0.2, cap=0.05, xlab='', ylab='PSC amplitude (pA)', main='', 
                    xrange=c(0.75,2.25), yrange=c(-400, 0), tick_length=0.2, 
                    x_tick_interval = NULL, y_tick_interval=100, lwd=1, 
                    type=6, amount=0.05, p.cex=0.5, filename='boxplot.svg', 
                    height=2.5, width=4, bg='transparent', alpha=0.6, na.rm=FALSE, save=FALSE){
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # data1 <- if ("s" %in% colnames(data)) data[, !colnames(data) %in% "s"] else data

  WBplot(data=data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main, xrange=xrange, yrange=yrange, 
         tick_length=tick_length, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, 
         lwd=lwd, type=type, na.rm=na.rm)
  
  set.seed(42)
  data$x_jitter <- jitter(data$x, amount=amount)
  
  # Set color with alpha transparency for points
  point_color <- rgb(169/255, 169/255, 169/255, alpha=alpha)  # darkgray with alpha=0.4
  points(data$x_jitter, data$y, pch=19, bg='transparent', col=point_color, lwd=lwd/3, cex=p.cex)

  if ("s" %in% colnames(data)) {

    # Connect data points within subjects with gray dotted lines
    line <- TRUE
    if (line) {
      subjects <- unique(data$s)
      for (subj in subjects) {
        subset_data <- data[data$s == subj, ]
        lines(subset_data$x_jitter, subset_data$y, col='darkgray', lwd=lwd, lty=3)  # lty=3 for dotted line
      }
    }
  }

  if (save) {
    dev.off()
  }
}


BoxPlot2_ver0 <- function(formula, data, wid = 0.2, cap = 0.05,
                     xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), tick_length = 0.2,
                     x_tick_interval = NULL, y_tick_interval = 100,
                     xlabel_angle = NULL,
                     lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
                     filename = 'boxplot.svg', height = 2.5, width = 4,
                     bg = 'transparent', alpha = 0.6, na.rm = FALSE, save = FALSE) {

  # formula
  response   <- as.character(formula[[2]])
  predictors <- all.vars(formula[[3]])
  if (!all(c(response, predictors) %in% names(data))) {
    stop('Response or predictor not found in data.')
  }

  if (length(predictors) == 1) {
    data$x <- factor(data[[predictors[1]]])
  } else {
    data$x <- interaction(data[[predictors[1]]], data[[predictors[2]]], sep=' : ')
  }
  data$y <- data[[response]]

  if (is.null(xrange)) {
    xrange <- range(as.numeric(data$x), na.rm=TRUE) + c(-wid, wid)
  }

  if (save) {
    svg(filename, width=width, height=height, bg=bg)
    on.exit(dev.off(), add=TRUE)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # suppress default x‐labels
  orig_axis <- graphics::axis
  assign('axis',
    function(side, at, labels = TRUE, tcl = NA, ...) {
      if (side == 1) {
        orig_axis(side, at = at, labels = FALSE, tcl = -tick_length, ...)
      } else {
        orig_axis(side, at = at, labels = labels, tcl = -tick_length, ...)
      }
    },
    envir = .GlobalEnv
  )
  on.exit(assign('axis', orig_axis, envir = .GlobalEnv), add = TRUE)

  # draw the boxplot (WBplot calls axis(1) & axis(2) internally)
  WBplot_ver0(data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, main = main,
         xrange = xrange, yrange = yrange, tick_length = tick_length,
         x_tick_interval = x_tick_interval, y_tick_interval = y_tick_interval,
         lwd = lwd, type = type, na.rm = na.rm)

  # restore axis()
  assign('axis', orig_axis, envir = .GlobalEnv)


  # draw custom x‐labels
  x_labels <- levels(data$x)

  # strip off everything from the first ' :'
  labs <- sub('\\s*:\\s*.*$', '', x_labels)

  # draw axis
  at <- seq_along(labs)
  if (is.null(xlabel_angle)) {
    axis(1, at = at, labels = labs, tcl = -tick_length, lwd = lwd)
  } else {
    axis(1, at = at, labels = FALSE, tcl = -tick_length, lwd = lwd)
    usr <- par('usr')
    y0  <- usr[3] - 0.1 * diff(usr[3:4])
    text(x = at, y = y0, labels = labs,
         srt = xlabel_angle, adj = 1, xpd = TRUE)
  }

  # —————— jitter & paired/unpaired block ——————
  set.seed(42)
  data$x_jitter <- jitter(as.numeric(data$x), amount = amount)
  
  # identify paired subjects (>=2 non-NA y’s)
  if ('s' %in% names(data)) {
    counts      <- ave(!is.na(data$y), data$s, FUN = sum)
    data$paired <- counts >= 2
  } else {
    data$paired <- rep(FALSE, nrow(data))
  }
  
  # draw lines for repeated measures (>=2)
  if ('s' %in% names(data)) {
    for (subj in unique(data$s[data$paired])) {
      sd <- subset(data, s == subj & !is.na(y))
      sd <- sd[order(as.numeric(sd$x)), ]
      lines(sd$x_jitter, sd$y,
            col = 'darkgray', lwd = lwd, lty = 3)
    }
  }
  
  # draw all paired points (or, if no 's', all points) in gray
  gray_col <- rgb(169/255,169/255,169/255, alpha = alpha)
  points(data$x_jitter[data$paired | !('s' %in% names(data))],
         data$y    [data$paired | !('s' %in% names(data))],
         pch = 19, col = gray_col, cex = p.cex, lwd = lwd/3)
  
  # draw the strictly unpaired points in red—but only if you really had 's'
  if ('s' %in% names(data)) {
    unpaired <- !data$paired & !is.na(data$y)
    red_col  <- rgb(205/255,92/255,92/255, alpha = alpha)
    points(data$x_jitter[unpaired],
           data$y    [unpaired],
           pch = 19, col = red_col, cex = p.cex, lwd = lwd/3)
  }
  
  if (save) dev.off()
}

BoxPlot2 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), tick_length = 0.2, x_tick_interval = NULL, y_tick_interval = 100,
                     xlabel_angle = NULL, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5, filename = 'boxplot.svg', 
                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, na.rm = FALSE, save = FALSE) {

  response   <- as.character(formula[[2]])
  predictors <- all.vars(formula[[3]])
  if (!all(c(response, predictors) %in% names(data))) stop('Response or predictor not found in data.')
  if (length(predictors) == 1) {
    data$x <- factor(data[[predictors]])
  } else {
    data$x <- interaction(data[[predictors[1]]], data[[predictors[2]]], sep=' : ')
  }
  data$y <- data[[response]]
  data_orig <- data

  if (log_y) {
    data$y[data$y <= 0] <- NA
    if (any(yrange <= 0, na.rm=TRUE)) stop('yrange must be > 0 when log_y=TRUE')
    data$y <- log10(data$y)
    yrange <- log10(yrange)
  }

  if (is.null(xrange)) {
    xrange <- range(as.numeric(data$x), na.rm=TRUE) + c(-wid, wid)
  }

  if (save) {
    svg(filename, width=width, height=height, bg=bg)
    on.exit(dev.off(), add=TRUE)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  orig_axis <- graphics::axis
  if (log_y) {
    assign('axis', function(...) {}, envir=.GlobalEnv)
  } else {
    assign('axis', function(side, at, labels=TRUE, tcl=NA, ...){
      if (side==1) orig_axis(side, at=at, labels=FALSE, tcl=-tick_length, ...)
      else         orig_axis(side, at=at, labels=labels, tcl=-tick_length, ...)
    }, envir=.GlobalEnv)
  }
  on.exit(assign('axis', orig_axis, envir=.GlobalEnv), add=TRUE)

  WBplot(data=data_orig, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main,
         xrange=xrange, yrange=yrange, tick_length=tick_length,
         x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval,
         lwd=lwd, type=type, log_y=log_y, na.rm=na.rm)

  assign('axis', orig_axis, envir=.GlobalEnv)

  labs <- sub('\\s*:\\s*.*$', '', levels(data$x))
  at   <- seq_along(labs)
  if (is.null(xlabel_angle)) {
    axis(1, at=at, labels=labs, tcl=-tick_length, lwd=lwd)
  } else {
    axis(1, at=at, labels=FALSE, tcl=-tick_length, lwd=lwd)
    usr <- par('usr'); y0 <- usr[3] - 0.1*diff(usr[3:4])
    text(x=at, y=y0, labels=labs, srt=xlabel_angle, adj=1, xpd=TRUE)
  }

  if (!('s' %in% names(data)) && 'animal_id' %in% names(data)) {
    data$s <- data$animal_id
  }

  set.seed(42)
  data$x_jitter <- jitter(as.numeric(data$x), amount=amount)
  if ('s' %in% names(data)) {
    counts      <- ave(!is.na(data$y), data$s, FUN=sum)
    data$paired <- counts >= 2
    for (subj in unique(data$s[data$paired])) {
      sd <- subset(data, s==subj & !is.na(y))
      sd <- sd[order(as.numeric(sd$x)), ]
      lines(sd$x_jitter, sd$y, col='darkgray', lwd=lwd, lty=3)
    }
  } else {
    data$paired <- FALSE
  }

  gray_col <- adjustcolor('#A9A9A9', alpha.f=alpha)
  red_col  <- adjustcolor('#CD5C5C', alpha.f=alpha)
  pts <- data$paired | !('s' %in% names(data))
  points(data$x_jitter[pts], data$y[pts], pch=19, col=gray_col, cex=p.cex, lwd=lwd/3)
  if ('s' %in% names(data)) {
    up <- !data$paired & !is.na(data$y)
    points(data$x_jitter[up], data$y[up], pch=19, col=red_col, cex=p.cex, lwd=lwd/3)
  }

  if (log_y) {
    major_tick_len <- -0.2
    minor_tick_len <- -0.1
    d0 <- floor(yrange[1]):ceiling(yrange[2])
    axis(2, at=d0, labels=10^d0, tcl=major_tick_len, las=1)
    minor_vals <- unlist(lapply(d0, function(d) (2:9)*10^d))
    mp <- log10(minor_vals)
    mp <- mp[mp >= yrange[1] & mp <= yrange[2]]
    axis(2, at=mp, labels=NA, tcl=minor_tick_len)
  }

  if (save) dev.off()
}

BoxPlot3_ver0 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), xlabel_angle = NULL, tick_length = 0.2, 
                     x_tick_interval = NULL, y_tick_interval = 100, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, na_rm_subjects = FALSE,
                     test_result, alpha_level = 0.05, group_names = NULL, sig_offset = NULL) {

  
  f_str <- deparse(formula)
  has_error <- grepl('Error', f_str)

  if (has_error) {
    err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var <- strsplit(err_part, '/')[[1]][1]
    subject_var <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula <- as.formula(main_formula_str)
  } else {
    subject_var <- NULL
    main_formula <- formula
  }

  response_var <- all.vars(formula(main_formula))[1]
  predictors <- all.vars(formula(main_formula))[-1]

  if (na_rm_subjects && !is.null(subject_var)) {
    data <- df[ !ave(is.na(data[[response_var]]), data[[subject_var]], FUN = any), ]
  }

  BoxPlot2_ver0(formula = formula, data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, xlabel_angle = xlabel_angle, main = main,
           xrange = xrange, yrange = yrange, tick_length = tick_length, y_tick_interval = y_tick_interval,
           lwd = lwd, type = type, amount = amount, p.cex = p.cex, height = height, width = width,
           bg = bg, na.rm = TRUE)

  has_error <- grepl("Error", deparse(formula))
  if (missing(test_result) || is.null(test_result) || nrow(test_result) == 0) {
    return()
  }

  # Prepare data$x & data$y
  response   <- response_var
  # predictors <- all.vars(formula[[3]])
  data$y     <- data[[response]]

  if (length(predictors) == 1) {
    data$x        <- factor(data[[predictors[1]]])
    single_factor <- TRUE
  } else {
    data$x        <- interaction(
      data[[predictors[1]]],
      data[[predictors[2]]],
      sep = ' : '
    )
    single_factor <- FALSE
  }

  x_labels    <- levels(data$x)
  x_positions <- setNames(seq_along(x_labels), x_labels)

  # compute base offset and tick
  y_span <- diff(range(data$y, na.rm = TRUE))
  offset <- if (is.null(sig_offset)) 0.05 * y_span else sig_offset
  tick   <- 0.25 * offset

  # Loop over tests
  for (i in seq_len(nrow(test_result))) {
    p_adj <- test_result$`p adjusted`[i]
    if (is.na(p_adj) || p_adj >= alpha_level) next

    parts <- strsplit(as.character(test_result$contrast[i]), ' vs ')[[1]]
    if (length(parts) != 2) next
    lev1 <- trimws(parts[1])
    lev2 <- trimws(parts[2])

    # Build the two labels that match interaction()
    if (single_factor) {
      label1 <- lev1
      label2 <- lev2
    } else {
      comp      <- as.character(test_result$comparison[i])
      outer_lev <- sub('.*? ([^ ]+) \\(.*', '\\1', comp)
      if (grepl('\\(paired\\)', comp)) {
        label1 <- paste0(lev1, ' : ', outer_lev)
        label2 <- paste0(lev2, ' : ', outer_lev)
      } else {
        label1 <- paste0(outer_lev, ' : ', lev1)
        label2 <- paste0(outer_lev, ' : ', lev2)
      }
    }

    if (!(label1 %in% x_labels) || !(label2 %in% x_labels)) next
    x1 <- x_positions[label1]
    x2 <- x_positions[label2]

    yvals <- c(data$y[data$x == label1], data$y[data$x == label2])
    yvals <- yvals[!is.na(yvals)]
    if (length(yvals) == 0) next
    y_max <- max(yvals)
    y_min <- min(yvals)

    # shift unpaired sig bars by 6*tick
    if (single_factor) {
      shift_amt <- 0
    } else {
      paired    <- grepl('\\(paired\\)', comp)
      shift_amt <- if (!paired) 6 * tick else 0
    }

    if (y_max > 0) {
      y_line <- y_max + offset + shift_amt
      segments(x1, y_line, x1, y_line - tick, lwd = lwd)
      segments(x2, y_line, x2, y_line - tick, lwd = lwd)
      text_y <- y_line + tick
    } else {
      y_line <- y_min - offset - shift_amt
      segments(x1, y_line, x1, y_line + tick, lwd = lwd)
      segments(x2, y_line, x2, y_line + tick, lwd = lwd)
      text_y <- y_line - tick
    }

    segments(x1, y_line, x2, y_line, lwd = lwd)
    text((x1 + x2) / 2, text_y, labels = '*', cex = 1.2)
  }
}

# BoxPlot3 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
#                      xrange = NULL, yrange = c(-400, 0), xlabel_angle = NULL, tick_length = 0.2, 
#                      x_tick_interval = NULL, y_tick_interval = 100, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
#                      height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, na_rm_subjects = FALSE,
#                      test_result, alpha_level = 0.05, group_names = NULL, sig_offset = NULL) {

#   f_str <- deparse(formula)
#   has_error <- grepl('Error', f_str)

#   if (has_error) {
#     err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
#     subject_var <- strsplit(err_part, '/')[[1]][1]
#     subject_var <- gsub('[[:space:]]', '', subject_var)
#     main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
#     main_formula <- as.formula(main_formula_str)
#   } else {
#     subject_var <- NULL
#     main_formula <- formula
#   }

#   response_var <- all.vars(formula(main_formula))[1]
#   predictors <- all.vars(formula(main_formula))[-1]

#   if (na_rm_subjects && !is.null(subject_var)) {
#     data <- df[ !ave(is.na(data[[response_var]]), data[[subject_var]], FUN = any), ]
#   }

#   BoxPlot2(formula = formula, data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, xlabel_angle = xlabel_angle, main = main,
#            xrange = xrange, yrange = yrange, tick_length = tick_length, y_tick_interval = y_tick_interval,
#            lwd = lwd, type = type, amount = amount, p.cex = p.cex, height = height, width = width, log_y=log_y,
#            bg = bg, na.rm = TRUE)

#   has_error <- grepl("Error", deparse(formula))
#   if (missing(test_result) || is.null(test_result) || nrow(test_result) == 0) {
#     return()
#   }

#   response   <- response_var
#   data$y     <- data[[response]]

#   if (length(predictors) == 1) {
#     data$x        <- factor(data[[predictors[1]]])
#     single_factor <- TRUE
#   } else {
#     data$x        <- interaction(
#       data[[predictors[1]]],
#       data[[predictors[2]]],
#       sep = ' : '
#     )
#     single_factor <- FALSE
#   }

#   x_labels    <- levels(data$x)
#   x_positions <- setNames(seq_along(x_labels), x_labels)

#   y_span <- diff(range(data$y, na.rm = TRUE))
#   offset <- if (is.null(sig_offset)) 0.05 * y_span else sig_offset
#   tick   <- 0.25 * offset

#   for (i in seq_len(nrow(test_result))) {
#     p_adj <- test_result$`p adjusted`[i]
#     if (is.na(p_adj) || p_adj >= alpha_level) next

#     parts <- strsplit(as.character(test_result$contrast[i]), ' vs ')[[1]]
#     if (length(parts) != 2) next
#     lev1 <- trimws(parts[1])
#     lev2 <- trimws(parts[2])

#     if (single_factor) {
#       label1 <- lev1
#       label2 <- lev2
#       paired <- FALSE
#     } else {
#       comp      <- as.character(test_result$comparison[i])
#       outer_lev <- sub('.*? ([^ ]+) \\(.*', '\\1', comp)
#       if (grepl('\\(paired\\)', comp)) {
#         label1 <- paste0(lev1, ' : ', outer_lev)
#         label2 <- paste0(lev2, ' : ', outer_lev)
#         paired <- TRUE
#       } else {
#         label1 <- paste0(outer_lev, ' : ', lev1)
#         label2 <- paste0(outer_lev, ' : ', lev2)
#         paired <- FALSE
#       }
#     }

#     if (!(label1 %in% x_labels) || !(label2 %in% x_labels)) next
#     x1 <- x_positions[label1]
#     x2 <- x_positions[label2]

#     if (log_y) {
#       yvals <- c(data[[response]][data$x == label1], data[[response]][data$x == label2])
#       yvals <- yvals[!is.na(yvals) & yvals > 0]
#       if (length(yvals) == 0) next
#       yvals <- log10(yvals)
#       y_max <- max(yvals)
#       y_min <- min(yvals)
#       y_span_log <- diff(range(log10(data[[response]][data[[response]] > 0]), na.rm = TRUE))
#       offset <- if (is.null(sig_offset)) 0.05 * y_span_log else sig_offset
#       tick   <- 0.25 * offset
#       shift_amt <- if (!single_factor && !paired) 6 * tick else 0
#       if (y_max > 0) {
#         y_line <- y_max + offset + shift_amt
#         segments(x1, y_line, x1, y_line - tick, lwd = lwd)
#         segments(x2, y_line, x2, y_line - tick, lwd = lwd)
#         text_y <- y_line + tick
#       } else {
#         y_line <- y_min - offset - shift_amt
#         segments(x1, y_line, x1, y_line + tick, lwd = lwd)
#         segments(x2, y_line, x2, y_line + tick, lwd = lwd)
#         text_y <- y_line - tick
#       }
#     } else {
#       yvals <- c(data$y[data$x == label1], data$y[data$x == label2])
#       yvals <- yvals[!is.na(yvals)]
#       if (length(yvals) == 0) next
#       y_max <- max(yvals)
#       y_min <- min(yvals)
#       shift_amt <- if (!single_factor && !paired) 6 * tick else 0
#       if (y_max > 0) {
#         y_line <- y_max + offset + shift_amt
#         segments(x1, y_line, x1, y_line - tick, lwd = lwd)
#         segments(x2, y_line, x2, y_line - tick, lwd = lwd)
#         text_y <- y_line + tick
#       } else {
#         y_line <- y_min - offset - shift_amt
#         segments(x1, y_line, x1, y_line + tick, lwd = lwd)
#         segments(x2, y_line, x2, y_line + tick, lwd = lwd)
#         text_y <- y_line - tick
#       }
#     }

#     segments(x1, y_line, x2, y_line, lwd = lwd)
#     text((x1 + x2) / 2, text_y, labels = '*', cex = 1.2)
#   }
# }

BoxPlot3 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), xlabel_angle = NULL, tick_length = 0.2, 
                     x_tick_interval = NULL, y_tick_interval = 100, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, log_y = FALSE, na_rm_subjects = FALSE,
                     test_result, alpha_level = 0.05, group_names = NULL, sig_offset = NULL) {

  f_str <- deparse(formula)
  has_error <- grepl('Error', f_str)

  if (has_error) {
    err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var <- strsplit(err_part, '/')[[1]][1]
    subject_var <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula <- as.formula(main_formula_str)
  } else {
    subject_var <- NULL
    main_formula <- formula
  }

  response_var <- all.vars(formula(main_formula))[1]
  predictors <- all.vars(formula(main_formula))[-1]

  if (na_rm_subjects && !is.null(subject_var)) {
    data <- df[ !ave(is.na(data[[response_var]]), data[[subject_var]], FUN = any), ]
  }

  BoxPlot2(formula = formula, data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, xlabel_angle = xlabel_angle, main = main,
           xrange = xrange, yrange = yrange, tick_length = tick_length, y_tick_interval = y_tick_interval,
           lwd = lwd, type = type, amount = amount, p.cex = p.cex, height = height, width = width, log_y=log_y,
           bg = bg, na.rm = TRUE)

  has_error <- grepl("Error", deparse(formula))
  if (missing(test_result) || is.null(test_result) || nrow(test_result) == 0) {
    return()
  }

  if ('p adjusted' %in% names(test_result)) {
    p_column <- 'p adjusted'
  } else if ('p value' %in% names(test_result)) {
    p_column <- 'p value'
  } else {
    stop("test_result must contain 'p adjusted' or 'p value'")
  }

  response <- response_var
  data$y <- data[[response]]

  if (length(predictors) == 1) {
    data$x <- factor(data[[predictors[1]]])
    single_factor <- TRUE
  } else {
    data$x <- interaction(
      data[[predictors[1]]],
      data[[predictors[2]]],
      sep = ' : '
    )
    single_factor <- FALSE
  }

  x_labels <- levels(data$x)
  x_positions <- setNames(seq_along(x_labels), x_labels)

  y_span <- diff(range(data$y, na.rm = TRUE))
  offset <- if (is.null(sig_offset)) 0.05 * y_span else sig_offset
  tick <- 0.25 * offset

  for (i in seq_len(nrow(test_result))) {
    p_value <- test_result[[p_column]][i]
    if (is.na(p_value) || p_value >= alpha_level) next

    label1 <- NULL
    label2 <- NULL
    paired <- FALSE

    contrast.available <- 'contrast' %in% names(test_result) &&
      !is.na(test_result$contrast[i]) &&
      nzchar(trimws(as.character(test_result$contrast[i])))

    parameter.available <- 'parameter' %in% names(test_result) &&
      !is.na(test_result$parameter[i]) &&
      nzchar(trimws(as.character(test_result$parameter[i])))

    if (contrast.available) {
      parts <- strsplit(as.character(test_result$contrast[i]), ' vs ')[[1]]
      if (length(parts) != 2) next

      lev1 <- trimws(parts[1])
      lev2 <- trimws(parts[2])

      if (single_factor) {
        label1 <- lev1
        label2 <- lev2
      } else {
        if (!'comparison' %in% names(test_result) ||
            is.na(test_result$comparison[i])) next

        comp <- as.character(test_result$comparison[i])
        outer_lev <- sub('.*? ([^ ]+) \\(.*', '\\1', comp)

        if (grepl('\\(paired\\)', comp)) {
          label1 <- paste0(lev1, ' : ', outer_lev)
          label2 <- paste0(lev2, ' : ', outer_lev)
          paired <- TRUE
        } else {
          label1 <- paste0(outer_lev, ' : ', lev1)
          label2 <- paste0(outer_lev, ' : ', lev2)
          paired <- FALSE
        }
      }
    } else if (parameter.available) {
      parameter <- trimws(as.character(test_result$parameter[i]))

      if (single_factor) {
        parts <- strsplit(parameter, ' - ', fixed=TRUE)[[1]]

        if (length(parts) == 2) {
          label1 <- trimws(parts[1])
          label2 <- trimws(parts[2])
        } else if (length(x_labels) == 2) {
          label1 <- x_labels[1]
          label2 <- x_labels[2]
        } else {
          next
        }
      } else {
        separator <- NULL

        if (grepl(' in ', parameter, fixed=TRUE)) {
          separator <- ' in '
          paired <- TRUE
        } else if (grepl(' for ', parameter, fixed=TRUE)) {
          separator <- ' for '
          paired <- FALSE
        }

        if (is.null(separator)) next

        parameter.parts <- strsplit(
          parameter,
          separator,
          fixed=TRUE
        )[[1]]

        if (length(parameter.parts) != 2) next

        parts <- strsplit(
          trimws(parameter.parts[1]),
          ' - ',
          fixed=TRUE
        )[[1]]

        if (length(parts) != 2) next

        lev1 <- trimws(parts[1])
        lev2 <- trimws(parts[2])
        outer_lev <- trimws(parameter.parts[2])

        labels.first <- c(
          paste0(lev1, ' : ', outer_lev),
          paste0(lev2, ' : ', outer_lev)
        )

        labels.second <- c(
          paste0(outer_lev, ' : ', lev1),
          paste0(outer_lev, ' : ', lev2)
        )

        if (all(labels.first %in% x_labels)) {
          label1 <- labels.first[1]
          label2 <- labels.first[2]
        } else if (all(labels.second %in% x_labels)) {
          label1 <- labels.second[1]
          label2 <- labels.second[2]
        } else {
          next
        }
      }
    } else if (single_factor && length(x_labels) == 2) {
      label1 <- x_labels[1]
      label2 <- x_labels[2]
    } else {
      next
    }

    if (!(label1 %in% x_labels) || !(label2 %in% x_labels)) next

    x1 <- x_positions[label1]
    x2 <- x_positions[label2]

    if (log_y) {
      yvals <- c(data[[response]][data$x == label1], data[[response]][data$x == label2])
      yvals <- yvals[!is.na(yvals) & yvals > 0]
      if (length(yvals) == 0) next

      yvals <- log10(yvals)
      y_max <- max(yvals)
      y_min <- min(yvals)
      y_span_log <- diff(range(log10(data[[response]][data[[response]] > 0]), na.rm = TRUE))
      offset <- if (is.null(sig_offset)) 0.05 * y_span_log else sig_offset
      tick <- 0.25 * offset
      shift_amt <- if (!single_factor && !paired) 6 * tick else 0

      if (y_max > 0) {
        y_line <- y_max + offset + shift_amt
        segments(x1, y_line, x1, y_line - tick, lwd = lwd)
        segments(x2, y_line, x2, y_line - tick, lwd = lwd)
        text_y <- y_line + tick
      } else {
        y_line <- y_min - offset - shift_amt
        segments(x1, y_line, x1, y_line + tick, lwd = lwd)
        segments(x2, y_line, x2, y_line + tick, lwd = lwd)
        text_y <- y_line - tick
      }
    } else {
      yvals <- c(data$y[data$x == label1], data$y[data$x == label2])
      yvals <- yvals[!is.na(yvals)]
      if (length(yvals) == 0) next

      y_max <- max(yvals)
      y_min <- min(yvals)
      shift_amt <- if (!single_factor && !paired) 6 * tick else 0

      if (y_max > 0) {
        y_line <- y_max + offset + shift_amt
        segments(x1, y_line, x1, y_line - tick, lwd = lwd)
        segments(x2, y_line, x2, y_line - tick, lwd = lwd)
        text_y <- y_line + tick
      } else {
        y_line <- y_min - offset - shift_amt
        segments(x1, y_line, x1, y_line + tick, lwd = lwd)
        segments(x2, y_line, x2, y_line + tick, lwd = lwd)
        text_y <- y_line - tick
      }
    }

    segments(x1, y_line, x2, y_line, lwd = lwd)
    text((x1 + x2) / 2, text_y, labels = '*', cex = 1.2)
  }
}

scatter_plot <- function(scatter, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, height=4, width=4, main='',
                         xlab = expression(A[fast] * " " * (pA)), ylab = expression(A[slow] * " " * (pA)),
                         colors=c("black", "indianred"), open_symbols=FALSE, lwd=1, p.cex=0.5, filename='scatter.svg', save=FALSE) {
  # Create a basic scatter plot
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Map levels to 1 and 2 alternately
  unique_levels <- unique(scatter$level)
  n <- length(unique_levels)
  scatter$level <- as.numeric(factor(scatter$level, levels = unique_levels, labels = rep(1:n, length.out = length(unique_levels))))

  # Determine plot symbols (only open circles if open_symbols is TRUE)
  pch <- if (open_symbols) 1 else 19
  cols <- hex_palette(n=n, color1=colors[1], color2=colors[2], reverse = FALSE)

  plot(scatter$x, scatter$y, col = cols[scatter$level], 
       pch = pch, cex = p.cex, xlim=xlim, ylim=ylim, xlab = xlab, ylab = ylab, 
       main = main, xaxt='n', yaxt='n', bty='n')

  # Add abline y = x as light gray; use segments as abline function extends beyond axes
  segments(min(xlim), min(ylim), max(xlim), max(ylim), col = "lightgray")

  # Define tick intervals and lengths
  x_ticks <- seq(min(xlim), max(xlim), by=x_tick_interval)
  y_ticks <- seq(min(ylim), max(ylim), by=y_tick_interval)
  tick_length <- -0.2

  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }
}

# Custom Theil-Sen function
theil_sen <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }
  
  n <- length(x)
  slopes <- c()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
  }
  
  slope <- median(slopes, na.rm=TRUE)
  intercept <- median(y - slope * x)
  
  return(list(slope=slope, intercept=intercept))
}

# Custom function to compute z-scores and p-values
compute_regression_stats <- function(estimate, bootstrap_estimates, conf_level=0.95) {
  # Standard error
  std_err <- sd(bootstrap_estimates, na.rm=TRUE)
  
  # Z-score
  z_value <- estimate / std_err
  
  # P-value (two-tailed)
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  
  # Confidence intervals
  lower_ci <- quantile(bootstrap_estimates, (1 - conf_level) / 2, na.rm=TRUE)
  upper_ci <- quantile(bootstrap_estimates, 1 - (1 - conf_level) / 2, na.rm=TRUE)
  
  return(list(
    estimate=estimate,
    std_err=std_err,
    z_value=z_value,
    p_value=p_value,
    lower_ci=lower_ci,
    upper_ci=upper_ci
  ))
}

# Function to format p-values, showing as 0.00000 when below threshold
format_p_value <- function(p, dp=5) {
  if (p < 10^(-dp)) {
    return(formatC(0, format='f', digits=dp))
  } else {
    return(formatC(p, format='f', digits=dp))
  }
}

# Custom Siegel Estimator function
siegel_sen <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }
  
  n <- length(x)
  slopes_per_point <- numeric(n)  # Vector to store median slopes for each point
  
  for (i in 1:n) {
    slopes <- c()
    for (j in 1:n) {
      if (i != j && x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
    slopes_per_point[i] <- median(slopes, na.rm=TRUE)  # Store the median slope for point i
  }
  
  slope <- median(slopes_per_point, na.rm=TRUE)  # Final slope is the median of individual medians
  intercept <- median(y - slope * x, na.rm=TRUE)
  
  return(list(slope=slope, intercept=intercept))
}

theil_sen_with_ci <- function(x, y, n_bootstrap=1e4, conf_level=0.95, seed=42, dp=5, n_points=1000) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }

  # Create a finer sequence for smoother CI plotting
  x1 <- seq(min(x), max(x), length.out=n_points)

  # Calculate Theil-Sen slope and intercept
  n <- length(x)
  slopes <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
  }

  slope <- median(slopes, na.rm=TRUE)
  intercept <- median(y - slope * x)

  # Bootstrap to estimate confidence intervals
  set.seed(seed)  # For reproducibility
  boot_slopes <- numeric(n_bootstrap)
  boot_intercepts <- numeric(n_bootstrap)
  preds <- matrix(NA, nrow=n_bootstrap, ncol=length(x1))

  for (b in 1:n_bootstrap) {
    sample_indices <- sample(1:n, replace=TRUE)
    x_boot <- x[sample_indices]
    y_boot <- y[sample_indices]
    
    slopes_boot <- c()
    if (length(unique(x_boot)) > 1) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          if (x_boot[j] != x_boot[i]) {
            slopes_boot <- c(slopes_boot, (y_boot[j] - y_boot[i]) / (x_boot[j] - x_boot[i]))
          }
        }
      }
      
      boot_slope <- median(slopes_boot, na.rm=TRUE)
      boot_intercept <- median(y_boot - boot_slope * x_boot)
      
      # Store bootstrap results
      boot_slopes[b] <- boot_slope
      boot_intercepts[b] <- boot_intercept
      
      # Predict y values for finer sequence x1
      preds[b, ] <- boot_intercept + boot_slope * x1
    }
  }

  # Remove NAs from bootstrap results
  boot_slopes <- boot_slopes[!is.na(boot_slopes)]
  boot_intercepts <- boot_intercepts[!is.na(boot_intercepts)]

  # Compute confidence intervals for predicted y values at x1
  ci <- apply(preds, 2, quantile, probs=c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2), na.rm=TRUE)

  # Compute regression stats for slope and intercept
  slope_stats <- compute_regression_stats(slope, boot_slopes, conf_level)
  intercept_stats <- compute_regression_stats(intercept, boot_intercepts, conf_level)

  # Dynamic confidence interval labels based on conf_level
  lower_ci_label <- paste0(round((1 - conf_level) / 2 * 100, 1), '%')
  upper_ci_label <- paste0(round((1 + conf_level) / 2 * 100, 1), '%')

  # Format results, rounding to the specified decimal places and applying p-value formatting
  summary_table <- data.frame(
    est=round(c(intercept_stats$estimate, slope_stats$estimate), dp),
    `se`=round(c(intercept_stats$std_err, slope_stats$std_err), dp),
    `z`=round(c(intercept_stats$z_value, slope_stats$z_value), dp),
    `P(>|z|)`=c(format_p_value(intercept_stats$p_value, dp), format_p_value(slope_stats$p_value, dp)),
    check.names=FALSE
  )
  
  # Add dynamically labeled confidence intervals to the table
  summary_table[lower_ci_label] <- round(c(intercept_stats$lower_ci, slope_stats$lower_ci), dp)
  summary_table[upper_ci_label] <- round(c(intercept_stats$upper_ci, slope_stats$upper_ci), dp)
  
  # Set row names
  rownames(summary_table) <- c('(intercept)', 'slope')
  
  return(list(
    summary_table=summary_table,
    lower_ci=ci[1, ],
    upper_ci=ci[2, ],
    predicted_y=intercept + slope * x1,
    x1=x1  # Return x1 for plotting purposes
  ))
}

# Custom Siegel Estimator with confidence intervals and regression stats
siegel_sen_with_ci <- function(x, y, n_bootstrap=1e4, conf_level=0.95, seed=42, dp=5, n_points=1000) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }

  # Create a finer sequence for smoother CI plotting
  x1 <- seq(min(x), max(x), length.out=n_points)

  n <- length(x)
  slopes_per_point <- numeric(n)
  
  # Original Siegel estimate
  for (i in 1:n) {
    slopes <- c()
    for (j in 1:n) {
      if (i != j && x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
    slopes_per_point[i] <- median(slopes, na.rm=TRUE)
  }
  
  slope <- median(slopes_per_point, na.rm=TRUE)
  intercept <- median(y - slope * x)

  # Bootstrap to estimate confidence intervals
  set.seed(seed)  # For reproducibility
  boot_slopes <- numeric(n_bootstrap)
  boot_intercepts <- numeric(n_bootstrap)
  preds <- matrix(NA, nrow=n_bootstrap, ncol=length(x1))

  for (b in 1:n_bootstrap) {
    sample_indices <- sample(1:n, replace=TRUE)
    x_boot <- x[sample_indices]
    y_boot <- y[sample_indices]
    
    slopes_per_point <- numeric(n)
    if (length(unique(x_boot)) > 1) {
      for (i in 1:n) {
        slopes <- c()
        for (j in 1:n) {
          if (i != j && x_boot[j] != x_boot[i]) {
            slopes <- c(slopes, (y_boot[j] - y_boot[i]) / (x_boot[j] - x_boot[i]))
          }
        }
        slopes_per_point[i] <- median(slopes, na.rm=TRUE)
      }
      
      boot_slope <- median(slopes_per_point, na.rm=TRUE)
      boot_intercept <- median(y_boot - boot_slope * x_boot)
      
      # Store bootstrap results
      boot_slopes[b] <- boot_slope
      boot_intercepts[b] <- boot_intercept
      
      # Predict y values for finer sequence x1
      preds[b, ] <- boot_intercept + boot_slope * x1
    }
  }

  # Remove NAs from bootstrap results
  boot_slopes <- boot_slopes[!is.na(boot_slopes)]
  boot_intercepts <- boot_intercepts[!is.na(boot_intercepts)]

  # Compute confidence intervals for predicted y values at x1
  ci <- apply(preds, 2, quantile, probs=c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2), na.rm=TRUE)

  # Compute regression stats for slope and intercept
  slope_stats <- compute_regression_stats(slope, boot_slopes, conf_level)
  intercept_stats <- compute_regression_stats(intercept, boot_intercepts, conf_level)

  # Dynamic confidence interval labels based on conf_level
  lower_ci_label <- paste0(round((1 - conf_level) / 2 * 100, 1), '%')
  upper_ci_label <- paste0(round((1 + conf_level) / 2 * 100, 1), '%')

  # Format results, rounding to the specified decimal places and applying p-value formatting
  summary_table <- data.frame(
    est=round(c(intercept_stats$estimate, slope_stats$estimate), dp),
    `se`=round(c(intercept_stats$std_err, slope_stats$std_err), dp),
    `z`=round(c(intercept_stats$z_value, slope_stats$z_value), dp),
    `P(>|z|)`=c(format_p_value(intercept_stats$p_value, dp), format_p_value(slope_stats$p_value, dp)),
    check.names=FALSE
  )
  
  # Add dynamically labeled confidence intervals to the table
  summary_table[lower_ci_label] <- round(c(intercept_stats$lower_ci, slope_stats$lower_ci), dp)
  summary_table[upper_ci_label] <- round(c(intercept_stats$upper_ci, slope_stats$upper_ci), dp)
  
  # Set row names
  rownames(summary_table) <- c('(intercept)', 'slope')
  
  # Return the summary table and the confidence intervals for predictions
  return(list(
    summary_table=summary_table,
    lower_ci=ci[1, ],
    upper_ci=ci[2, ],
    predicted_y=intercept + slope * x1,
    x1=x1  # Return x1 for plotting purposes
  ))
}

convert_A_to_scatter <- function(A, sign=1) {
  # Get the unique 's' values
  s_vals <- unique(A$s)
  
  # Initialize empty vectors for x and y
  x_vals <- numeric(length(s_vals))
  y_vals <- numeric(length(s_vals))
  
  # Loop through each 's' and assign the correct values to x and y
  for (i in seq_along(s_vals)) {
    x_vals[i] <- A$y[A$s == s_vals[i] & A$x == 1]
    y_vals[i] <- A$y[A$s == s_vals[i] & A$x == 2]
  }
  
  # Create the scatter data frame
  scatter <- data.frame(s=s_vals, x=sign*x_vals, y=sign*y_vals)
  
  return(scatter)
}

convert_to_scatter <- function(A, sign=1) {
  # Split data frame A by 's' and 'x' values to get unique levels
  split_data <- split(A, A$s)
  
  # Create scatter data frame
  scatter <- do.call(rbind, lapply(split_data, function(df) {
    n <- nrow(df) / 2  # Calculate the number of levels
    data.frame(
      s = df$s[1:n],
      level = df$x[1:n],
      x = sign * (df$y[1:n]),
      y = sign *(df$y[(n + 1):(2 * n)])
    )
  }))
  
  return(scatter)
}

create2condition_df <- function(mat, cols = c(1, 10), var_name = 'amplitude', 
  levels = list(condition = c('fast', 'slow'), cell_type = NULL),
  start_id = 1) {

  if (length(cols) != 2) stop('cols must contain exactly 2 column indices')

  n <- nrow(mat)
  total <- 2 * n

  df <- data.frame(
    s = factor(rep(seq_len(n) + start_id - 1, times = 2)),
    value = c(mat[, cols[1]], mat[, cols[2]])
  )

  for (level_name in names(levels)) {
    lv <- levels[[level_name]]
    if (is.null(lv)) next
    if (length(lv) == 1) {
      df[[level_name]] <- factor(rep(lv, total))
    } else if (length(lv) == 2) {
      df[[level_name]] <- factor(rep(lv, each = n), levels = lv)
    } else {
      stop(paste('Level', level_name, 'must be of length 1 or 2'))
    }
  }

  colnames(df)[colnames(df) == 'value'] <- var_name
  df <- df[, c(setdiff(names(df), var_name), var_name)]
  return(df)
}

# ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
#   height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, lwd=1, p.cex=0.5, 
#   filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), 
#   save=FALSE, bg='transparent', dp=3, return.output=FALSE) {

#   # scatter <- convert_A_to_scatter(A=A, sign=sign)
#   scatter <- convert_to_scatter(A=A, sign=sign)

#   # Create the scatter plot
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
#   if ('level' %in% colnames(scatter)) {
#     unique_levels <- unique(scatter$level)
#     n <- length(unique_levels)
#     scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
#   } else {
#     n <- 1  # If 'level' column does not exist
#     scatter$level <- rep(1, dim(scatter)[1])
#   }

#   # Determine plot symbols (only open circles if open_symbols is TRUE)
#   pch <- if (open_symbols) 1 else 19

#   cols <- hex_palette(n=n, color1=colors[1], color2=colors[2], reverse=FALSE)

#   plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
#        xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

#   # Add regression line or non-parametric line based on 'reg' parameter
#   if (!reg) {
#     segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
#   } else {
#     reg.method <- match.arg(reg.method)

#     if (reg.method == 'Siegel') {
#       reg_func <- if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen
#     } else if (reg.method == 'Theil-Sen') {
#       reg_func <- if (plot.CI || return.output) theil_sen_with_ci else theil_sen
#     }

#     # Perform regression and extract results
#     if (plot.CI || return.output) {
#       reg_results <- reg_func(scatter$x, scatter$y, n_bootstrap=reg.CI.settings$nboot, 
#                               conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
      
#       summary_table <- reg_results$summary_table
#       intercept <- summary_table["(intercept)", "est"]
#       slope <- summary_table["slope", "est"]
      
#       # Confidence intervals for predictions
#       lower_ci <- reg_results$lower_ci
#       upper_ci <- reg_results$upper_ci
#       predicted_y <- reg_results$predicted_y
#       x1 <- reg_results$x1
#     } else {
#       out <- reg_func(scatter$x, scatter$y)
#       intercept <- out$intercept
#       slope <- out$slope
#       x1 <- seq(min(scatter$x), max(scatter$x), length.out=1000)  # Fallback for plotting
#       predicted_y <- intercept + slope * x1
#     }

#     # Plot the regression line using x1
#     lines(x1, predicted_y, lwd=lwd, col=reg.color, lty=3)

#     # Plot confidence intervals if requested
#     if (plot.CI) {
#       lines(x1, lower_ci, col=reg.color, lty=2)
#       lines(x1, upper_ci, col=reg.color, lty=2)
#     }
#   }

#   # Define tick intervals and lengths
#   x_ticks <- seq(min(xlim), max(xlim), by=x_tick_interval)
#   y_ticks <- seq(min(ylim), max(ylim), by=y_tick_interval)
  
#   # Customize x-axis
#   axis(1, at=x_ticks, labels=x_ticks, tcl=-tick_length, lwd=lwd)

#   # Customize y-axis with horizontal labels
#   axis(2, at=y_ticks, labels=y_ticks, tcl=-tick_length, las=1, lwd=lwd)

#   if (save) {
#     dev.off()
#   }

#   if (reg && return.output) {
#     return(summary_table)
#   }
# }


# ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
#                         height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, 
#                         lwd=1, p.cex=0.5, filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', 
#                         reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), save=FALSE, 
#                         bg='transparent', dp=3, return.output=FALSE) {
  
#   # Convert A to scatter
#   scatter <- convert_to_scatter(A=A, sign=sign)

#   # Create the scatter plot
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
#   if ('level' %in% colnames(scatter)) {
#     unique_levels <- unique(scatter$level)
#     n <- length(unique_levels)
#     scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
#   } else {
#     n <- 1  # If 'level' column does not exist
#     scatter$level <- rep(1, dim(scatter)[1])
#   }

#   # Determine plot symbols
#   pch <- if (open_symbols) 1 else 19
#   cols <- hex_palette(n=2, color1=colors[1], color2=colors[2], reverse=FALSE)
  
#   # Plot scatter points
#   plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
#        xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

#   # Initialize list to store summary tables for each level
#   summary_tables <- list()

#   # Add regression line or non-parametric line based on 'reg' parameter
#   if (!reg) {
#     segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
#   } else {
#     reg.method <- match.arg(reg.method)
#     levels <- unique(scatter$level)

#     reg_func <- switch(reg.method,
#                        'Siegel' = if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen,
#                        'Theil-Sen' = if (plot.CI || return.output) theil_sen_with_ci else theil_sen)

#     for (level in levels) {
#       # Filter data for the current level
#       level_data <- scatter[scatter$level==level,]
      
#       if (plot.CI || return.output) {
#         # Perform regression with confidence intervals for the current level
#         reg_results <- reg_func(level_data$x, level_data$y, n_bootstrap=reg.CI.settings$nboot, 
#                                 conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
        
#         # Store the summary table for this level
#         summary_tables[[as.character(level)]] <- reg_results$summary_table
        
#         intercept <- reg_results$summary_table["(intercept)", "est"]
#         slope <- reg_results$summary_table["slope", "est"]
        
#         # Confidence intervals for predictions
#         lower_ci <- reg_results$lower_ci
#         upper_ci <- reg_results$upper_ci
#         predicted_y <- reg_results$predicted_y
#         x1 <- reg_results$x1
#       } else {
#         # Perform standard regression without confidence intervals for the current level
#         out <- reg_func(level_data$x, level_data$y)
#         intercept <- out$intercept
#         slope <- out$slope
#         x1 <- seq(min(level_data$x), max(level_data$x), length.out=reg.points)
#         predicted_y <- intercept + slope * x1
#       }
      
#       # Plot the regression line for the current level
#       lines(x1, predicted_y, lwd=lwd, col=cols[as.integer(level)], lty=3)
      
#       # Plot confidence intervals if requested
#       if (plot.CI) {
#         lines(x1, lower_ci, col=cols[as.integer(level)], lty=2)
#         lines(x1, upper_ci, col=cols[as.integer(level)], lty=2)
#       }
#     }
#   }

#   # Customize axes
#   axis(1, at=seq(min(xlim), max(xlim), by=x_tick_interval), tcl=-tick_length, lwd=lwd)
#   axis(2, at=seq(min(ylim), max(ylim), by=y_tick_interval), las=1, tcl=-tick_length, lwd=lwd)

#   if (save) {
#     dev.off()
#   }

#   # Return the list of summary tables if requested
#   if (reg && return.output) {
#     return(summary_tables)
#   }
# }


# ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
#                         height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, 
#                         lwd=1, p.cex=0.5, filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', 
#                         reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), save=FALSE, 
#                         show_pairs=FALSE, bg='transparent', dp=3, return.output=FALSE) {
  
#   # Convert A to scatter
#   scatter <- convert_to_scatter(A=A, sign=sign)

#   # Create the scatter plot
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
#   if ('level' %in% colnames(scatter)) {
#     unique_levels <- unique(scatter$level)
#     n <- length(unique_levels)
#     scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
#   } else {
#     n <- 1  # If 'level' column does not exist
#     scatter$level <- rep(1, dim(scatter)[1])
#   }

#   # Determine plot symbols
#   pch <- if (open_symbols) 1 else 19
#   cols <- hex_palette(n=2, color1=colors[1], color2=colors[2], reverse=FALSE)
  
#   # Plot scatter points
#   plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
#        xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

#   # Initialize list to store summary tables for each level
#   summary_tables <- list()

#   # Add regression line or non-parametric line based on 'reg' parameter
#   if (!reg) {
#     segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
#   } else {
#     reg.method <- match.arg(reg.method)
#     levels <- unique(scatter$level)

#     reg_func <- switch(reg.method,
#                        'Siegel' = if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen,
#                        'Theil-Sen' = if (plot.CI || return.output) theil_sen_with_ci else theil_sen)

#     for (level in levels) {
#       # Filter data for the current level
#       level_data <- scatter[scatter$level==level,]
      
#       if (plot.CI || return.output) {
#         # Perform regression with confidence intervals for the current level
#         reg_results <- reg_func(level_data$x, level_data$y, n_bootstrap=reg.CI.settings$nboot, 
#                                 conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
        
#         # Store the summary table for this level
#         summary_tables[[as.character(level)]] <- reg_results$summary_table
        
#         intercept <- reg_results$summary_table["(intercept)", "est"]
#         slope <- reg_results$summary_table["slope", "est"]
        
#         # Confidence intervals for predictions
#         lower_ci <- reg_results$lower_ci
#         upper_ci <- reg_results$upper_ci
#         predicted_y <- reg_results$predicted_y
#         x1 <- reg_results$x1
#       } else {
#         # Perform standard regression without confidence intervals for the current level
#         out <- reg_func(level_data$x, level_data$y)
#         intercept <- out$intercept
#         slope <- out$slope
#         x1 <- seq(min(level_data$x), max(level_data$x), length.out=reg.points)
#         predicted_y <- intercept + slope * x1
#       }
      
#       # Plot the regression line for the current level
#       lines(x1, predicted_y, lwd=lwd, col=cols[as.integer(level)], lty=3)
      
#       # Plot confidence intervals if requested
#       if (plot.CI) {
#         lines(x1, lower_ci, col=cols[as.integer(level)], lty=2)
#         lines(x1, upper_ci, col=cols[as.integer(level)], lty=2)
#       }
#     }
#   }

#   if (show_pairs && "s" %in% colnames(scatter)) {
#     subjects <- unique(scatter$s)
#     for (subj in subjects) {
#       subj_data <- scatter[scatter$s == subj, ]
#       # Optionally, order the data by x or some grouping variable:
#       # subj_data <- subj_data[order(subj_data$x), ]
#       # Connect them with lines
#       if (nrow(subj_data) > 1) {
#         lines(subj_data$x, subj_data$y, col='darkgray', lty=3, lwd=lwd)
#       }
#     }
#   }

#   # Customize axes
#   axis(1, at=seq(min(xlim), max(xlim), by=x_tick_interval), tcl=-tick_length, lwd=lwd)
#   axis(2, at=seq(min(ylim), max(ylim), by=y_tick_interval), las=1, tcl=-tick_length, lwd=lwd)

#   if (save) {
#     dev.off()
#   }

#   # Return the list of summary tables if requested
#   if (reg && return.output) {
#     return(summary_tables)
#   }
# }
# Define the start and end colors of your palette Slate Blue to Indian Red
hex_palette <- function(n, color1='#6A5ACD', color2='#CD5C5C', reverse = FALSE) {
  
  # Create a sequence of colors
  colors <- colorRampPalette(c(color1, color2))(n)
  
  # reverse colors if reverse=TRUE
  if (reverse) {
    colors <- rev(colors)
  }
  
  return(colors)
}

ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
                        height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, 
                        lwd=1, p.cex=0.5, filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', 
                        reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), save=FALSE, 
                        show_pairs=FALSE, bg='transparent', dp=3, return.output=FALSE, log_xy=FALSE) {
  
  scatter <- convert_to_scatter(A=A, sign=sign)
  
  if (log_xy) {
    if (xlim[1] == 0) xlim[1] <- 0.1
    if (ylim[1] == 0) ylim[1] <- 0.1
    scatter <- scatter[scatter$x > 0 & scatter$y > 0, ]
    scatter$x <- log10(scatter$x)
    scatter$y <- log10(scatter$y)
    xlim_plot <- log10(xlim)
    ylim_plot <- log10(ylim)
  } else {
    xlim_plot <- xlim
    ylim_plot <- ylim
  }
  
  if (save) {
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }
  
  if ('level' %in% colnames(scatter)) {
    unique_levels <- unique(scatter$level)
    n <- length(unique_levels)
    scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
  } else {
    n <- 1
    scatter$level <- rep(1, dim(scatter)[1])
  }
  
  pch <- if (open_symbols) 1 else 19
  cols <- hex_palette(n=2, color1=colors[1], color2=colors[2], reverse=FALSE)
  
  plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim_plot, ylim=ylim_plot, 
       xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')
  
  summary_tables <- list()
  
  if (!reg) {
    if (log_xy) {
      abline(0, 1, lwd=lwd, col=reg.color, lty=3)
    } else {
      segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)
    }
  } else {
    reg.method <- match.arg(reg.method)
    levels <- unique(scatter$level)
    reg_func <- switch(reg.method,
                       'Siegel' = if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen,
                       'Theil-Sen' = if (plot.CI || return.output) theil_sen_with_ci else theil_sen)
    for (level in levels) {
      level_data <- scatter[scatter$level==level,]
      if (plot.CI || return.output) {
        reg_results <- reg_func(level_data$x, level_data$y, n_bootstrap=reg.CI.settings$nboot, 
                                conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
        summary_tables[[as.character(level)]] <- reg_results$summary_table
        intercept <- reg_results$summary_table["(intercept)", "est"]
        slope <- reg_results$summary_table["slope", "est"]
        lower_ci <- reg_results$lower_ci
        upper_ci <- reg_results$upper_ci
        predicted_y <- reg_results$predicted_y
        x1 <- reg_results$x1
      } else {
        out <- reg_func(level_data$x, level_data$y)
        intercept <- out$intercept
        slope <- out$slope
        x1 <- seq(min(level_data$x), max(level_data$x), length.out=reg.points)
        predicted_y <- intercept + slope * x1
      }
      lines(x1, predicted_y, lwd=lwd, col=cols[as.integer(level)], lty=3)
      if (plot.CI) {
        lines(x1, lower_ci, col=cols[as.integer(level)], lty=2)
        lines(x1, upper_ci, col=cols[as.integer(level)], lty=2)
      }
    }
  }
  
  if (show_pairs && "s" %in% colnames(scatter)) {
    subjects <- unique(scatter$s)
    for (subj in subjects) {
      subj_data <- scatter[scatter$s == subj, ]
      if (nrow(subj_data) > 1) {
        lines(subj_data$x, subj_data$y, col='darkgray', lty=3, lwd=lwd)
      }
    }
  }
  
  if (log_xy) {
    at_x <- log10(outer(1:9, 10^(floor(xlim_plot[1]):ceiling(xlim_plot[2])), "*"))
    at_x <- at_x[at_x >= xlim_plot[1] & at_x <= xlim_plot[2]]
    labs_x <- 10^at_x
    major_x <- which(log10(labs_x) %% 1 == 0)
    axis(1, at=at_x[major_x], labels=format(labs_x[major_x], trim=TRUE, scientific=FALSE), lwd=lwd, tcl=-tick_length)
    axis(1, at=at_x, labels=FALSE, lwd=lwd, tcl=-0.5*tick_length)
    at_y <- log10(outer(1:9, 10^(floor(ylim_plot[1]):ceiling(ylim_plot[2])), "*"))
    at_y <- at_y[at_y >= ylim_plot[1] & at_y <= ylim_plot[2]]
    labs_y <- 10^at_y
    major_y <- which(log10(labs_y) %% 1 == 0)
    axis(2, at=at_y[major_y], labels=format(labs_y[major_y], trim=TRUE, scientific=FALSE), las=1, lwd=lwd, tcl=-tick_length)
    axis(2, at=at_y, labels=FALSE, lwd=lwd, tcl=-0.5*tick_length)
  } else {
    x_ticks <- seq(min(xlim), max(xlim), by=x_tick_interval)
    y_ticks <- seq(min(ylim), max(ylim), by=y_tick_interval)
    axis(1, at=x_ticks, labels=x_ticks, tcl=-tick_length, lwd=lwd)
    axis(2, at=y_ticks, labels=y_ticks, las=1, tcl=-tick_length, lwd=lwd)
  }
  
  if (save) {
    dev.off()
  }
  
  if (reg && return.output) {
    return(summary_tables)
  }
}

traces2plot <- function(V, traces, offsets=NULL, color1='#6A5ACD', color2='#CD5C5C', xlim=NULL, ylim=NULL, lwd=1, 
  xbar=100, ybar=50, reverse=TRUE, show_text=FALSE, normalise=FALSE, height=3, width=3, 
  filename='traces.svg', dt=0.1, save=FALSE){
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  n <- length(traces)
  colors <- hex_palette(n=n, color1=color1, color2=color2, reverse=reverse)
  if (reverse) traces = rev(traces)
  if (is.null(offsets)) offsets <- rep(0, n)

  x = 0:dt:(dim(V)[1]-1)*dt

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  if (is.null(ylim)) ylim <- c(0, max(apply(V[, traces],2,max)))

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  y <- V[, traces[1]]
  
  # # Design a low-pass filter
  # fs <- 10  # Sampling frequency (e.g., 1 Hz if data points are 1 second apart)
  # cutoff <- 0.1  # Cutoff frequency (e.g., 0.1 Hz)
  # butter_order <- 5  # Order of the Butterworth filter
  # bf <- butter(butter_order, cutoff / (0.5 * fs), type = "low")

  # # Apply the low-pass filter
  # y_smooth <- filtfilt(bf, y)


  plot(x[idx1:idx2]+offsets[1], y[idx1:idx2], type='l', col=colors[1], xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab = '', ylab = '')

  # Loop through remaining traces and add them to the plot
  for (i in 2:n) {
    y <- V[, traces[i]]
    lines(x[idx1:idx2]+offsets[i], y[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }

  # Define scale bar lengths and ybar position
  # ybar_start <- (max(ylim) - min(ylim)) / 10
  ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20

  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar

  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise){
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }

  # Add labels to the scale bars
  if (show_text){
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start -xbar/4, y = (y_start + y_end) / 2, labels = paste(ybar, 'mV'), adj = c(0.5, 0.5), srt = 90)
  }

  if (save) {
    dev.off()
  }

}


fun_single_example <- function(rawdata, fits, start_time=50, baseline=50, idx=1, model='product2', response_sign_method = c('smooth', 'regression', 'cumsum')){
  
  ind1 <- (start_time - baseline)/dx
  ind2 <- start_time/dx

  y <- rawdata[, idx+1]
  y <- na.omit(y)
  y <- y[ind1:length(y)]
  y <- y - mean(y[ind1:ind2])
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  # sign <- sign_fun(y, direction_method=response_sign_method) 
  # single_example <- data.frame('x'=x, 'y'= sign * y)
  single_example <- data.frame('x'=x, 'y'= y)

  if (model == 'product') {
    fun.2.fit <- product1
  } else if (model == 'product2') {
    fun.2.fit <- product2
  }

  pars <- fits[[idx]]$fits
  if (model == 'product') {
    pars[4] <- pars[4] + baseline
  } else if (model == 'product2') {
    pars[4] <- pars[4] + baseline
    pars[8] <- pars[8] + baseline
  }    
  yfit <- fun.2.fit(pars, x) 
  single_example$y_fit <- yfit

  if (model == 'product2'){
    yfit1 <- product1(pars[1:4], x) 
    yfit2 <- product1(pars[5:8], x) 
    single_example$y_fit1 <- yfit1
    single_example$y_fit2 <- yfit2
  }
  return(single_example)
}

# single_fit_egs <- function(traces, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, height=4, width=2.5, xbar=100, ybar=50, log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
#   if (save) {
#     # Open SVG device
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }
  
#   x <- traces$x
#   y <- traces$y
  
#   fit1 <- traces$yfit1
#   if (identical(func, product2N)){
#     fit2 <- traces$yfit2
#   }
  
#   if (identical(func, product3N)){
#     fit2 <- traces$yfit2
#     fit3 <- traces$yfit3
#   }

#   if (is.null(xlim)) xlim <- c(min(x), max(x))
  
#   if (is.null(ylim)) {
#     if (log_y) {
#       y <- -y
#       fit1 <- -fit1
#       if (identical(func, product2N)){
#         fit2 <- -fit2
#       }
#       if (identical(func, product3N)){
#         fit2 <- -fit2
#         fit3 <- -fit3
#       }

#       # Define custom major tick positions
#       y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
#       log_y_ticks <- log(y_ticks)
#       valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

#       # Set ylim based on valid ticks, using the lowest major tick for the lower bound
#       ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
#     } else {
#       ylim <- c(-max(y, na.rm=TRUE), 0)
#     }
#   }

#   if (log_y) {
#     y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
#     fit1 <- ifelse(fit1 > 0, log(fit1), NA)
#     if (identical(func, product2N)){
#       fit2 <- ifelse(fit2 > 0, log(fit2), NA)
#     }
#     if (identical(func, product3N)){
#       fit2 <- ifelse(fit2 > 0, log(fit2), NA)
#       fit3 <- ifelse(fit3 > 0, log(fit3), NA)
#     }
#   } 

#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))

#   plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

#   if (identical(func, product1N)){
#     fits <- cbind(fit)
#   }else if (identical(func, product2N)){
#     fits <- cbind(fit1, fit2)
#   }else if (identical(func, product3N)){
#     fits <- cbind(fit1, fit2, fit3)
#   }

  
#   # Loop through remaining traces and add them to the plot
#   for (i in 1:dim(fits)[2]) {
#     y_fit <- fits[, i]
#     lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
#   }
  
#   # Define scale bar lengths and ybar position
#   ybar <- ifelse(log_y, exp(1), ybar)
#   ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
#   # Add scale bars at the bottom right
#   x_start <- max(xlim) - xbar - 50
#   y_start <- ybar_start
#   x_end <- x_start + xbar
#   y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
#   # Draw the scale bars
#   segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
#   if (!normalise) {
#     segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
#   }
  
#   # Add labels to the scale bars
#   if (show_text) {
#     text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
#     if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, 'pA')), adj = c(0.5, 0.5), srt = 90)
#   }
  
#   # Add the y-axis only if log_y is TRUE
#   if (log_y) {
#     tick_length <- -0.2
#     minor_tick_length <- -0.1
    
#     # Major tick positions and labels for the log scale
#     y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
#     log_y_ticks <- log(y_ticks)
#     valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

#     # Minor tick positions for the log scale
#     minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
#                        20, 30, 40, 50, 60, 70, 80, 90, 
#                        200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
#     log_minor_y_ticks <- log(minor_y_ticks)
#     valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

#     # Add major ticks
#     axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

#     # Add minor ticks
#     axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
#     # Add y-axis label
#     mtext('PSC amplitude (pA)', side=2, line=2.5)
#   }
  
#   if (save) {
#     dev.off()
#   }
# }

single_fit_egs <- function(traces, sign=-1, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, main='', height=4, width=2.5, 
  xbar=100, ybar=50, scalebar_x_offset=50, ybar_units = 'pA', log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  x <- traces$x
  y <- traces$y
  
  fit1 <- if (identical(func, product1N)) traces$yfit else traces$yfit1
  
  if (identical(func, product2N)){
    fit2 <- traces$yfit2
  }
  
  if (identical(func, product3N)){
    fit2 <- traces$yfit2
    fit3 <- traces$yfit3
  }

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  
  if (is.null(ylim)) {
    if (log_y) {
      y <- sign * y
      fit1 <- sign * fit1
      if (identical(func, product2N)){
        fit2 <- sign * fit2
      }
      if (identical(func, product3N)){
        fit2 <- sign * fit2
        fit3 <- sign * fit3
      }

      # Define custom major tick positions
      y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
      log_y_ticks <- log(y_ticks)
      valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

      # Set ylim based on valid ticks, using the lowest major tick for the lower bound
      ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
    } else {
      ylim <- c(-max(y, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
    fit1 <- ifelse(fit1 > 0, log(pmax(fit1, .Machine$double.eps)), NA)
    if (identical(func, product2N)){
      fit2 <- ifelse(fit2 > 0, log(pmax(fit2, .Machine$double.eps)), NA)
    }
    if (identical(func, product3N)){
      fit2 <- ifelse(fit2 > 0, log(pmax(fit2, .Machine$double.eps)), NA)
      fit3 <- ifelse(fit3 > 0, log(pmax(fit3, .Machine$double.eps)), NA)
    }
    y[ y < 0 ] <- NA
    y[which.min(abs(x - xlim[1]))] <- 0
  } 

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', main=main, xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  if (identical(func, product1N)){
    fits <- cbind(fit1)
  }else if (identical(func, product2N)){
    fits <- cbind(fit1, fit2)
  }else if (identical(func, product3N)){
    fits <- cbind(fit1, fit2, fit3)
  }

  if (log_y) {
    fits[ fits < 0 ] <- NA
    for (col_idx in seq_len(ncol(fits))) {
      col_data <- fits[, col_idx]
      
      # Find the first non-NA value's index
      f <- which(!is.na(col_data))[1]
      
      if (!is.na(f) && f > 0) {
        # Check if row f-1 is NA
        if (is.na(col_data[f - 1])) {
          col_data[f - 1] <- 0
        }
      }
      
      fits[, col_idx] <- col_data
    }
  }

  # Loop through remaining traces and add them to the plot
  for (i in 1:dim(fits)[2]) {
    y_fit <- fits[, i]
    lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }
  
  # Define scale bar lengths and ybar position
  ybar <- ifelse(log_y, exp(1), ybar)
  ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - scalebar_x_offset
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise) {
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, ybar_units)), adj = c(0.5, 0.5), srt = 90)
  }
  
  # Add the y-axis only if log_y is TRUE
  if (log_y) {
    tick_length <- -0.2
    minor_tick_length <- -0.1
    
    # Major tick positions and labels for the log scale
    y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
    log_y_ticks <- log(y_ticks)
    valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

    # Minor tick positions for the log scale
    minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
                       20, 30, 40, 50, 60, 70, 80, 90, 
                       200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
    log_minor_y_ticks <- log(minor_y_ticks)
    valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

    # Add major ticks
    axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

    # Add minor ticks
    axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
    # Add y-axis label
    if (grepl("A$", ybar_units)) {
      mtext(bquote(PSC~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    } else if (grepl("V$", ybar_units)) {
      mtext(bquote(PSP~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    }

  }
  
  if (save) {
    dev.off()
  }
}

# single examples
single_egs <- function(traces, sign=-1, xlim=NULL, ylim=NULL, main='', lwd=1.2, width=5, height=5, log_y=FALSE, 
  bg='transparent', xbar=100, ybar=50, ybar_units = 'pA', colors = c('#A9A9A9', '#D3D3D3'), filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  x <- traces$time
  ymat <- as.matrix(traces[, colnames(traces) != 'time', drop = FALSE])

  if (is.null(xlim)) xlim <- c(min(x), max(x))

  if (is.null(ylim)) {
    if (log_y) {
      ymat <- sign * ymat
      # Define custom major tick positions
      y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
      log_y_ticks <- log(y_ticks)
      valid_ticks <- log_y_ticks[log_y_ticks <= log(max(ymat[ymat > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

      # Set ylim based on valid ticks, using the lowest major tick for the lower bound
      ylim <- c(min(valid_ticks), log(max(ymat[ymat > 0], na.rm=TRUE)))
    } else {
      ylim <- c(min(ymat, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    ymat <- apply(ymat, 2, function(y) {
      y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
      y[y < 0] <- NA
      y[which.min(abs(x - xlim[1]))] <- 0
      y
    })
  }

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  n <- ncol(ymat)
  cols <- colorRampPalette(colors)(n)

  plot(x[idx1:idx2], ymat[idx1:idx2, ncol(ymat)], col = cols[n], xlab = '', ylab = '', xlim = xlim, ylim = ylim,
       type = 'l', bty = 'l', las = 1, axes = FALSE, frame = FALSE, lwd = lwd, main = main)

  # then plot the rest in reverse order
  if (ncol(ymat) > 1){
    for (ii in (ncol(ymat) - 1):1) {
      lines(x[idx1:idx2], ymat[idx1:idx2, ii], col = cols[ii])
    }
  }

  ybar <- ifelse(log_y, exp(1), ybar)
  ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)

  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar

  text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
  text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold", paste(ybar, ybar_units)), adj = c(0.5, 0.5), srt = 90)

  # Add the y-axis only if log_y is TRUE
  if (log_y) {
    tick_length <- -0.2
    minor_tick_length <- -0.1
    
    # Major tick positions and labels for the log scale
    y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
    log_y_ticks <- log(y_ticks)
    valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

    # Minor tick positions for the log scale
    minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
                       20, 30, 40, 50, 60, 70, 80, 90, 
                       200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
    log_minor_y_ticks <- log(minor_y_ticks)
    valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

    # Add major ticks
    axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

    # Add minor ticks
    axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
    # Add y-axis label
    if (grepl("A$", ybar_units)) {
      mtext(bquote(PSC~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    } else if (grepl("V$", ybar_units)) {
      mtext(bquote(PSP~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    }
  }
  
  if (save) {
    dev.off()
  }
}

SingleFitExample <- function(traces, xlim=NULL, ylim=NULL, ylab='PSC amplitude (pA)', tick_length=0.2, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, 
  height=4, width=2.5, xbar=100, ybar=50, log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }
  
  x <- traces$x
  y <- traces$y
  
  fit1 <- traces$yfit1
  if (identical(func, product2N)){
    fit2 <- traces$yfit2
  }
  
  if (identical(func, product3N)){
    fit2 <- traces$yfit2
    fit3 <- traces$yfit3
  }

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  
  if (is.null(ylim)) {
    if (log_y) {
      y <- -y
      fit1 <- -fit1
      if (identical(func, product2N)){
        fit2 <- -fit2
      }
      if (identical(func, product3N)){
        fit2 <- -fit2
        fit3 <- -fit3
      }

      # Define custom major tick positions
      y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
      log_y_ticks <- log(y_ticks)
      valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

      # Set ylim based on valid ticks, using the lowest major tick for the lower bound
      ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
    } else {
      ylim <- c(-max(y, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
    fit1 <- ifelse(fit1 > 0, log(fit1), NA)
    if (identical(func, product2N)){
      fit2 <- ifelse(fit2 > 0, log(fit2), NA)
    }
    if (identical(func, product3N)){
      fit2 <- ifelse(fit2 > 0, log(fit2), NA)
      fit3 <- ifelse(fit3 > 0, log(fit3), NA)
    }
  } 

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  if (identical(func, product1N)){
    fits <- cbind(fit)
  }else if (identical(func, product2N)){
    fits <- cbind(fit1, fit2)
  }else if (identical(func, product3N)){
    fits <- cbind(fit1, fit2, fit3)
  }

  
  # Loop through remaining traces and add them to the plot
  for (i in 1:dim(fits)[2]) {
    y_fit <- fits[, i]
    lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }
  
  # Define scale bar lengths and ybar position
  ybar <- ifelse(log_y, exp(1), ybar)
  ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise) {
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, 'pA')), adj = c(0.5, 0.5), srt = 90)
  }
  
  # Add the y-axis only if log_y is TRUE
  if (log_y) {
    tick_length <- -tick_length
    minor_tick_length <- -tick_length/2
    
    # Major tick positions and labels for the log scale
    y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
    log_y_ticks <- log(y_ticks)
    valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

    # Minor tick positions for the log scale
    minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
                       20, 30, 40, 50, 60, 70, 80, 90, 
                       200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
    log_minor_y_ticks <- log(minor_y_ticks)
    valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

    # Add major ticks
    axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

    # Add minor ticks
    axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
    # Add y-axis label
    mtext(ylab, side=2, line=2.5)
  }
  
  if (save) {
    dev.off()
  }
}


trace_extract <- function(rawdata, start_time=50, baseline=50, idx=1, response_sign_method = c('smooth', 'regression', 'cumsum')){
  ind1 <- (start_time - baseline)/dx
  ind2 <- start_time/dx

  y <- rawdata[, idx+1]
  y <- na.omit(y)
  y <- y[ind1:length(y)]
  y <- y - mean(y[ind1:ind2])
  # x <- 0:dx:(length(y)-1)*dx

  sign <- sign_fun(y, direction_method=response_sign_method) 
  return(sign * y)
}

# Function to process each sheet with user-defined base name
process_sheet <- function(sheet_name, summary, data, dt, stimulation_time, baseline, smooth, base_name, func=product2) {
  # Get the summary list
  summary_list <- get(paste0(base_name, '_summary'), envir = .GlobalEnv)
  summary_list[[sheet_name]] <- summary
  assign(paste0(base_name, '_summary'), summary_list, envir = .GlobalEnv)
  
  # Process fits
  fits <- t(sapply(1:length(summary), function(ii){
    X <- summary[[ii]]$output
    as.vector(t(X))
  }))
  
  # Create new column names by appending 1 and 2 to the original names
  if (identical(func, product1)){
    new_colnames <- rep(colnames(summary[[1]]$output), 2)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')   
  }else if (identical(func, product2)){
    new_colnames <- rep(colnames(summary[[1]]$output), 2)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')
  }else if (identical(func, product3)){
    new_colnames <- rep(colnames(summary[[1]]$output), 3)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2', 'A3')
  }
  
  colnames(fits) <- new_colnames
  rownames(fits) <- 1:dim(fits)[1]
  
  fits_list <- get(paste0(base_name, '_fits'), envir = .GlobalEnv)
  fits_list[[sheet_name]] <- fits
  assign(paste0(base_name, '_fits'), fits_list, envir = .GlobalEnv)
  
  # Process peaks
  peaks <- sapply(2:dim(data[[sheet_name]])[2], function(ii){
    x <- data[[sheet_name]][,1]
    y <- data[[sheet_name]][,ii]
    
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
  })
  
  peaks_list <- get(paste0(base_name, '_peaks'), envir = .GlobalEnv)
  peaks_list[[sheet_name]] <- peaks
  assign(paste0(base_name, '_peaks'), peaks_list, envir = .GlobalEnv)
}



# # for analysis of noise baseline to determoine cut-off frequencies

# FFT_plot <- function(response, response0=NULL, dx=0.1, xlim=c(0, 3000), ylim=c(10^-5, 10^2), xlab='frequency (Hz)', ylab='amplitude (pA)', col='gray', lwd=1, width=6, height=6, filename="plot.svg", save=FALSE, freq4baseline=25, spar=0.85) {
  
#   y <- response

#   N <- length(y)
#   fft_result <- fft(y)
#   fft_magnitude <- Mod(fft_result) / N

#   # frequency axis
#   sampling_rate <- 1e3 / dx
#   frequencies <- seq(0, sampling_rate / 2, length.out = N / 2 + 1)

#   # Plot the frequency spectrum
#   X <- frequencies
#   Y <- fft_magnitude[1:(N / 2 + 1)]

#   # Filter frequencies within the desired range
#   in_range <- X >= xlim[1] & X <= xlim[2]
#   X <- X[in_range]
#   Y <- Y[in_range]

#   # Compute gain in dB if response0 is provided
#   if (!is.null(response0)) {
#     N0 <- length(response0)
#     fft_result0 <- fft(response0)
#     fft_magnitude0 <- Mod(fft_result0) / N0
#     Y0 <- fft_magnitude0[1:(N0 / 2 + 1)][in_range]
#     gain_db <- 20 * log10(Y / Y0)
#     ylab <- 'gain (dB)'
#     Y <- gain_db
#     if (is.null(ylim)) ylim <- range(c(10 * floor(min(Y) / 10), 1, gain_db))  # Adjust ylim for gain in dB
#     log_scale <- ""
#   } else {
#     log_scale <- "y"
#   }

#   if (save) {
#     # Open SVG device
#     svg(file = filename, width = width, height = height, bg = 'transparent')
#   } else {
#     dev.new(width = width, height = height, noRStudioGD = TRUE)
#   }

#   # Plot the FFT or gain
#   plot(X, Y, col = col, type = 'l', bty = 'l', las = 1, lwd = lwd, 
#        main = '', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = log_scale, axes = FALSE, frame = FALSE)

#   # Adding custom ticks and grid
#   axis(1, at = seq(xlim[1], xlim[2], by = 1000), tcl = -0.2)  # x-axis with custom tick length

#   if (is.null(response0)) {
#     # Major ticks (decades) for amplitude scale
#     major_ticks <- 10^(log10(ylim[1]):log10(ylim[2]))
#     # Minor ticks (in-between decades)
#     minor_ticks <- c(2:9 %o% 10^(log10(ylim[1]):log10(ylim[2])))

#     # Custom y-axis labels using expressions
#     major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))

#     # Add major ticks and labels
#     axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las = 1)

#     # Add minor ticks without labels
#     axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
#   } else {
#     # Major ticks for dB scale
#     major_ticks <- seq(10 * floor(min(Y) / 10), 0, by = 10)
#     axis(2, at = major_ticks, tcl = -0.2, las = 1)
#   }

#   # Adding grid
#   grid()

#   # Add a smooth spline interpolation line for the periodogram
#   if (is.null(response0)) {
#     smooth_fit <- smooth.spline(X, Y, spar = spar)
#     lines(smooth_fit, col = 'slateblue', lty = 3, lwd = 2)

#     # Calculate drop point at 1/sqrt(2) based on spline fit
#     smooth_values <- predict(smooth_fit, X)$y
#     drop_point <- 1 / sqrt(2) * mean(smooth_values[X<freq4baseline])
#     crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    
#     # Linear interpolation to find more accurate crossing points
#     if (length(crossings) > 0) {
#       intersection_freqs <- numeric()
#       for (i in crossings) {
#         x1 <- X[i]
#         x2 <- X[i + 1]
#         y1 <- smooth_values[i]
#         y2 <- smooth_values[i + 1]
#         intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
#       }

#       if (length(intersection_freqs) > 0) {
#         median_freq <- median(intersection_freqs)  # Find the median frequency
#         print(paste0('cutoff frequency is ', round(median_freq,2)))
#         # Draw lines from axes to the median point
#         segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
#         segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
#       }
#     }
#   }

#   # Add -3 dB line if response0 is provided
#   if (!is.null(response0)) {
#     # Calculate -3 dB point
#     minus_3db <- -3
#     crossings <- which(diff(sign(Y - minus_3db)) != 0)  # Find all crossing points
    
#     # Linear interpolation to find more accurate crossing points
#     intersection_freqs <- numeric()
#     for (i in crossings) {
#       x1 <- X[i]
#       x2 <- X[i + 1]
#       y1 <- Y[i]
#       y2 <- Y[i + 1]
#       intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (minus_3db - y1) / (y2 - y1))
#     }

#     if (length(intersection_freqs) > 0) {
#       median_freq <- median(intersection_freqs)  # Find the median frequency
#       print(paste0('cutoff frequency is ', round(median_freq,2)))
#       # Draw lines from axes to the median point
#       segments(x0 = xlim[1], y0 = minus_3db, x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
#       segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
#     }
#   }

#   if (save) {
#     dev.off()
#   }
# }


FFT_plot <- function(response, response0 = NULL, dx = 0.1, xlim = c(0, 3000), ylim = c(10^-5, 10^2), 
                     xlab = 'frequency (Hz)', ylab = 'amplitude (pA)', col = 'gray', lwd = 1, 
                     width = 6, height = 6, filename = "plot.svg", save = FALSE, 
                     freq4baseline = 25, spar = 0.85, logx = FALSE) {
  
  y <- response
  N <- length(y)
  fft_result <- fft(y)
  fft_magnitude <- Mod(fft_result) / N

  # frequency axis
  sampling_rate <- 1e3 / dx
  frequencies <- seq(0, sampling_rate / 2, length.out = N / 2 + 1)

  # Plot the frequency spectrum
  X <- frequencies
  Y <- fft_magnitude[1:(N / 2 + 1)]

  # Filter frequencies within the desired range
  if (logx) xlim <- c(max(1, xlim[1]), xlim[2])
  in_range <- X >= xlim[1] & X <= xlim[2]
  X <- X[in_range]
  Y <- Y[in_range]

  # Compute gain in dB if response0 is provided
  if (!is.null(response0)) {
    N0 <- length(response0)
    fft_result0 <- fft(response0)
    fft_magnitude0 <- Mod(fft_result0) / N0
    Y0 <- fft_magnitude0[1:(N0 / 2 + 1)][in_range]
    gain_db <- 20 * log10(Y / Y0)
    ylab <- 'gain (dB)'
    Y <- gain_db
    if (is.null(ylim)) ylim <- range(c(10 * floor(min(Y) / 10), 1, gain_db))  # Adjust ylim for gain in dB
    log_scale <- ""
  } else {
    log_scale <- "y"
  }

  if (logx) {
    # Filter out zero values to avoid log(0)
    nonzero_idx <- X > 0
    X <- X[nonzero_idx]
    Y <- Y[nonzero_idx]
    # X <- log10(X)
    xlim <- c(max(1, xlim[1]), xlim[2])
    log_scale <- if (log_scale == "y") "xy" else "x"
  }

  if (save) {
    # Open SVG device
    svg(file = filename, width = width, height = height, bg = 'transparent')
  } else {
    dev.new(width = width, height = height, noRStudioGD = TRUE)
  }

  # Plot the FFT or gain
  plot(X, Y, col = col, type = 'l', bty = 'l', las = 1, lwd = lwd, 
       main = '', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = log_scale, axes = FALSE, frame = FALSE)

  # Adding custom ticks and grid
  if (logx) {
    major_ticks <- 10^(log10(xlim[1]):log10(xlim[2]))

    # Minor ticks (in-between decades)
    minor_ticks <- c(2:9 %o% 10^(log10(xlim[1]):log10(xlim[2])))
    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))
    # Add major ticks and labels
    axis(1, at = major_ticks, labels = major_labels, tcl = -0.2)
    # Add minor ticks without labels
    axis(1, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    axis(1, at = seq(xlim[1], xlim[2], by = 1000), tcl = -0.2)  # x-axis with custom tick length
  }

  if (is.null(response0)) {
    # Major ticks (decades) for amplitude scale
    major_ticks <- 10^(log10(ylim[1]):log10(ylim[2]))
    # Minor ticks (in-between decades)
    minor_ticks <- c(2:9 %o% 10^(log10(ylim[1]):log10(ylim[2])))
    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))
    # Add major ticks and labels
    axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las = 1)
   # Add minor ticks without labels
    axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    # Major ticks for dB scale
    major_ticks <- seq(10 * floor(min(Y) / 10), 0, by = 10)
    axis(2, at = major_ticks, tcl = -0.2, las = 1)
  }

  # Adding grid
  grid()

  # Add a smooth spline interpolation line for the periodogram
  if (is.null(response0)) {
    smooth_fit <- smooth.spline(X, Y, spar = spar)
    lines(smooth_fit, col = 'slateblue', lty = 3, lwd = 2)

    # Calculate drop point at 1/sqrt(2) based on spline fit
    smooth_values <- predict(smooth_fit, X)$y
    drop_point <- if (logx)  1 / sqrt(2) * mean(smooth_values[X < log10(freq4baseline)]) else 1 / sqrt(2) * mean(smooth_values[X<freq4baseline])
    crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    if (length(crossings) > 0) {
      intersection_freqs <- numeric()
      for (i in crossings) {
        x1 <- X[i]
        x2 <- X[i + 1]
        y1 <- Y[i]     # smooth_values[i]
        y2 <- Y[i + 1] # smooth_values[i + 1]
        intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
      }

      if (length(intersection_freqs) > 0) {
        median_freq <- median(intersection_freqs)  # Find the median frequency
        # Draw lines from axes to the median point
        segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        print(paste0('cutoff frequency is ', round(median_freq, 2)))       
      }
    }
  }

  # Add -3 dB line if response0 is provided
  if (!is.null(response0)) {
    # Calculate -3 dB point
    drop_point <- if (logx)  mean(Y[X < log10(freq4baseline)]) - 3 else 1 / sqrt(2) * mean(Y[X<freq4baseline])
    # if logx and input tends to output then drop_point is same as minus_3db i.e. -3
    crossings <- which(diff(sign(Y - drop_point)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    intersection_freqs <- numeric()
    for (i in crossings) {
      x1 <- X[i]
      x2 <- X[i + 1]
      y1 <- Y[i]
      y2 <- Y[i + 1]
      intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
    }

    if (length(intersection_freqs) > 0) {
      median_freq <- median(intersection_freqs)  # Find the median frequency
      # Draw lines from axes to the median point
      segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      print(paste0('cutoff frequency is ', round(median_freq, 2)))      
    }
  }

  if (save) {
    dev.off()
  }
}


# for power gain_db <- 10 * log10(power / original_power)
periodogram_plot <- function(response, response0=NULL, dx=0.1, xlim=c(0, 3000), ylim=c(10^-10, 10^2), xlab='frequency (Hz)', ylab='power (pA^2/Hz)', col='gray', lwd=1, width=6, height=6, filename='plot.svg', save=FALSE, spar=0.85, freq4baseline=25) {
  
  y <- response
  fs <- 1e3 / dx  # Compute sampling frequency

  # Compute periodogram for the response
  periodogram <- spectrum(y, plot = FALSE)
  freqs <- periodogram$freq * fs
  power <- periodogram$spec

  # Filter frequencies within the desired range
  in_range <- freqs >= xlim[1] & freqs <= xlim[2]
  freqs <- freqs[in_range]
  power <- power[in_range]

  # Compute gain in dB if response0 is provided
  if (!is.null(response0)) {
    original_periodogram <- spectrum(response0, plot = FALSE)
    original_power <- original_periodogram$spec[in_range]
    gain_db <- 10 * log10(power / original_power)
    ylab <- 'gain (dB)'
    power <- gain_db
    if (is.null(ylim)) ylim <- range(c(10*floor(min(power)/10), 1, gain_db))  # Adjust ylim for gain in dB
    log_scale <- ""
  } else {
    log_scale <- "y"
  }

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Plot the periodogram or gain
  plot(freqs, power, col=col, type = 'l', bty='l', las=1, lwd=lwd, 
       main='', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, log=log_scale, axes=FALSE, frame=FALSE)

  # Adding custom ticks and grid
  axis(1, at = seq(xlim[1], xlim[2], by = 500), tcl = -0.2)  # x-axis with custom tick length

  if (is.null(response0)) {
    # Major ticks (decades) for linear scale
    major_ticks <- 10^(floor(log10(ylim[1])):ceiling(log10(ylim[2])))
    # Minor ticks (in-between decades)
    minor_ticks <- c(outer(2:9, 10^(floor(log10(ylim[1])):ceiling(log10(ylim[2])))))

    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))

    # Add major ticks and labels
    axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las=1)

    # Add minor ticks without labels
    axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    # Major ticks for dB scale
    major_ticks <- seq(10*floor(min(power)/10), 0, by=10)
    axis(2, at = major_ticks, tcl = -0.2, las=1)
  }

  # Adding grid
  grid()

  # Add a smooth spline interpolation line for the periodogram
  if (is.null(response0)){
    smooth_fit <- smooth.spline(freqs, power, spar = spar)
    lines(smooth_fit, col='slateblue', lty=3, lwd=2)

    # Calculate drop point at 1/sqrt(2) based on spline fit
    smooth_values <- predict(smooth_fit, freqs)$y
    drop_point <- 1 / 2 * mean(smooth_values[freqs<freq4baseline])
    crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    # Linear interpolation to find more accurate crossing points
    if (length(crossings) > 0) {
      intersection_freqs <- numeric()
      for (i in crossings) {
        x1 <- freqs[i]
        x2 <- freqs[i + 1]
        y1 <- smooth_values[i]
        y2 <- smooth_values[i + 1]
        intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
      }
    if (length(intersection_freqs) > 0) {
        median_freq <- median(intersection_freqs)  # Find the median frequency
        print(paste0('cutoff frequency is ', round(median_freq,2)))
        # Draw lines from axes to the median point
        segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      }
    }
  }
  # Add -3 dB line if response0 is provided
  if (!is.null(response0)) {
    # Calculate -3 dB point
    minus_3db <- -3
    crossings <- which(diff(sign(power - minus_3db)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    intersection_freqs <- numeric()
    for (i in crossings) {
      x1 <- freqs[i]
      x2 <- freqs[i + 1]
      y1 <- power[i]
      y2 <- power[i + 1]
      intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (minus_3db - y1) / (y2 - y1))
    }

    if (length(intersection_freqs) > 0) {
      median_freq <- median(intersection_freqs)  # Find the median frequency
      print(paste0('cutoff frequency is ', round(median_freq,2)))
      # Draw lines from axes to the median point
      segments(x0 = xlim[1], y0 = minus_3db, x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
      segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
    }
  }

  if (save) {
    dev.off()
  }
}

processed_signals_plot <- function(input, output, width=6, height=6, main=''){
  dev.new(width=width, height=height, noRStudioGD=TRUE)
  plot(t * 1e3, input, type = 'l', col = 'darkgray', xlab = '', ylab = '', main = main, axes = FALSE, frame = FALSE, ylim = range(c(input, output)))
  lines(t * 1e3, output, col = 'lightgray')
  lines(c(max(t) * 1e3 - 1000, max(t) * 1e3), c(-3, -3), col = 'black', lwd = 2)
  lines(c(max(t) * 1e3 - 1000, max(t) * 1e3 - 1000), c(-3, -1), col = 'black', lwd = 2)
}


FFT_plot2 <- function(output, input=NULL, t, fs, width=12, height=6, xlab = 'frequency [Hz]',  ylab = 'Vm/Vc', spar = 0.05) {
  compute_fft <- function(signal, N) {
    fft_result <- fft(signal)
    freq <- (0:(length(t) - 1)) * (fs / length(t))
    valid_indices <- freq < 3e3
    list(
      fft = fft_result[valid_indices],
      freq = freq[valid_indices]
    )
  }
  
  add_epsilon <- function(magnitude, N, epsilon = 1e-6) {
    abs(magnitude) / N + epsilon
  }
  
  find_cutoff_frequency <- function(positive_freq, positive_vm_vc_ratio) {
    spline_fit <- smooth.spline(positive_freq, positive_vm_vc_ratio, spar = spar)
    fitted_spline <- predict(spline_fit, positive_freq)$y
    max_response <- max(positive_vm_vc_ratio)
    cutoff_level <- max_response / sqrt(2)
    cutoff_index <- which.min(abs(fitted_spline - cutoff_level))
    positive_freq[cutoff_index]
  }
  
  plot_signals <- function(t, input, output) {
    plot(t * 1e3, input, type = "l", col = "darkgray", xlab = "", ylab = "", main = "", axes = FALSE, frame = FALSE, ylim = range(c(input, output)))
    lines(t * 1e3, output, col = "lightgray")
    lines(c(max(t) * 1e3 - 1000, max(t) * 1e3), c(-3, -3), col = "black", lwd = 2)
    lines(c(max(t) * 1e3 - 1000, max(t) * 1e3 - 1000), c(-3, -1), col = "black", lwd = 2)
  }
  
  plot_frequency_response <- function(freq, vm_vc_ratio, fitted_spline, cutoff_frequency, xlab, ylab) {
    plot(freq, vm_vc_ratio, type = "l", col = "slateblue", log = "x", xlab = xlab, ylab = ylab, main = "", frame = FALSE)
    lines(freq, fitted_spline, col = "indianred", lty = "dotted")
    abline(v = cutoff_frequency, col = "darkgray", lty = "dotted")
  }
  
  if (!is.null(input)) {
    N <- length(input)
    input_fft_res <- compute_fft(input, N)
    output_fft_res <- compute_fft(output, N)
    
    input_fft_magnitude <- add_epsilon(input_fft_res$fft, N)
    output_fft_magnitude <- add_epsilon(output_fft_res$fft, N)
    
    vm_vc_ratio <- output_fft_magnitude / input_fft_magnitude
    positive_freq_indices <- input_fft_res$freq > 0
    positive_freq <- input_fft_res$freq[positive_freq_indices]
    positive_vm_vc_ratio <- vm_vc_ratio[positive_freq_indices]
    
    cutoff_frequency <- find_cutoff_frequency(positive_freq, positive_vm_vc_ratio)
    
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    par(mfrow = c(1, 2))  
    plot_signals(t, input, output)
    plot_frequency_response(positive_freq, positive_vm_vc_ratio, predict(smooth.spline(positive_freq, positive_vm_vc_ratio, spar = spar), positive_freq)$y, cutoff_frequency, xlab, ylab)
    
    print(sprintf("-3 dB cutoff frequency is approximately %.2f Hz", cutoff_frequency))
  } else {
    N <- length(output)
    output_fft_res <- compute_fft(output, N)
    
    output_fft_magnitude <- add_epsilon(output_fft_res$fft, N)
    positive_freq_indices <- output_fft_res$freq > 0
    positive_freq <- output_fft_res$freq[positive_freq_indices]
    positive_output_fft_magnitude <- output_fft_magnitude[positive_freq_indices]
    
    cutoff_frequency <- find_cutoff_frequency(positive_freq, positive_output_fft_magnitude)
    
    par(mfrow = c(1, 2))
    plot_signals(t, output, output)
    plot_frequency_response(positive_freq, positive_output_fft_magnitude, predict(smooth.spline(positive_freq, positive_output_fft_magnitude, spar = spar), positive_freq)$y, cutoff_frequency, xlab, ylab)
    
    print(sprintf("-3 dB cutoff frequency is approximately %.2f Hz", cutoff_frequency))
  }
}


# functions for voltage step analysis

# function to plot step for visualisation
step_plot <- function(x, y, tstep=c(50, 250), xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=FALSE) {
  
  xfit <- x[x<sum(tstep)]
  yfit <- y[x<sum(tstep)]
  
  xlim <- c(10*floor(min(xfit)/10), 10*ceiling(max(xfit)/10))

  dev.new(width=width, height=height, noRStudioGD=TRUE)
  plot(xfit, yfit, col='gray', xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, type='l', bty='l', las=1, lwd=lwd, main='')
  abline(v = tstep[1], col='darkgray', lwd=lwd, lty=3) 
  abline(v = tstep[2], col='darkgray', lwd=lwd, lty=3) 

}

vc_step_fit <- function(response, dt=0.1, tstep=c(50, 250), dV=5, t_interval=25, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
  xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, response_sign_method='smooth', lwd=1.2, width=8, height=4, 
  return.output=FALSE, show.output=TRUE, show.plot=TRUE) {

  method <- match.arg(method)
  fun.2.fit <- function(params,x){
    -params[1]*exp(-params[2]*x) -params[3]*exp(-params[4]*x) + params[1] + params[3] + params[5]
  }
  form.2.fit <- y ~ -a*exp(-b*x) -c*exp(-d*x) + a + c + e

  y <- response
  x <- seq(0, length(y)-1) * dt

  X <- x[x<(tstep[2]-dt)]
  Y <- y[x<(tstep[2]-dt)]

  sign <- sign_fun(Y, direction_method=response_sign_method) 
  Y <- sign * Y

  idx1 <- which.max(Y)
  A <- Y[idx1]

  A.fit <- 0.9 * A
  idx2 <-which.min(abs(Y[idx1:length(Y)] - A.fit)) + idx1 - 1
  xfit <- X[idx2:length(X)]
  yfit <- Y[idx2:length(X)]

  output <- nFIT(response=yfit, dt=dt, func=fun.2.fit, method='LM', lower=NULL, upper=NULL, filter=FALSE, return.output=TRUE, show.output=FALSE, show.plot=FALSE)

  fits <- output$fits
  traces <- output$traces

  n <- idx2-idx1
  x2fit <- c(rev(-seq(1:n)) * dt, traces$x)
  y2fit <- fun.2.fit(fits, x2fit)  

  traces=data.frame(x=x2fit - x2fit[1] + tstep[1], y=sign * Y[idx1:length(X)], yfit=sign * y2fit)
    
  dI <- fits[1] + fits[3] + fits[5]
  Rm <- dV/dI * 1e3
  Rs <- dV/y2fit[1] * 1e3
  Cm <- abs(fits[1])/fits[2]/dV + abs(fits[3])/fits[4]/dV # area under given exponetial of form a * exp(-bx) is a/b pAms ie fC then divide by 5 mV gives pC
  
  # to check area use trapz in pracma package:
  Cm_est <- ( pracma:::trapz(x2fit, y2fit) - dI * (max(x2fit) - min(x2fit)) ) / dV
  
  df_output <- data.frame(Cm=Cm, Cm_est=Cm_est, Rm=Rm, Rs=Rs)
  row.names(df_output) <- ''

  # Display the output
  df_output

  if (show.plot){
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    par(mfrow = c(1, 2))
    
    step_plot2(x=x[x<sum(tstep)], y=y[x<sum(tstep)], x2fit=traces$x, y2fit=traces$yfit, tstep=tstep)

    ylim <- c(5*floor(min(traces$y)/5), 0)
    step_plot2(x=traces$x[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], y=traces$y[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], 
      x2fit=traces$x[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], y2fit=traces$yfit[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], tstep=tstep, ylim=ylim, round2=5)
  }

  if (show.output){
    print(df_output)
  }
  
  if (return.output){
    return(list(pars=df_output, fits=fits))
  }

}

step_plot2 <- function(x, y, x2fit, y2fit, tstep=c(50, 250), xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, lwd=1.2, filter=FALSE, width=5, height=5, tcl = -0.2, round2=10, save=FALSE) {
  xlim <- c(round2*floor(min(x)/round2), round2*ceiling(max(x)/round2))
  ylim <- if (is.null(ylim)) c(round2*floor(min(y2fit)/round2), round2*ceiling(max(y)/round2)) else ylim
  
  plot(x, y, col='darkgray', xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, type='l', bty='l', las=1, lwd=lwd, main='',
       tcl=tcl, xaxs='i', yaxs='i')
  
  lines(x2fit, y2fit, col='indianred', lty=3, lwd=lwd)
}



# Function to wait for a mouse click
wait_for_click <- function() {
  cat("Click on the plot to continue\n")
  locator(1)
}


# Function to generate initial start values
generate_start_values <- function(y, x) {
  xmin <- min(x) * runif(1, 0.9, 1.1)
  Imax <- max(y) * runif(1, 0.9, 1.1)
  tau <- runif(1, 0, 1)
  list(xmin = xmin, Imax = Imax, tau=tau)
}

# Function to fit the model with a fallback to fixing n = 1 if fitting fails
fit_fun <- function(y, x, start_values, lower_bounds, upper_bounds) {
  fit <- tryCatch({
    nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
          start = start_values, 
          lower = lower_bounds, 
          upper = upper_bounds, 
          control = nls.lm.control(maxiter = 500))
  }, error = function(e) {
    cat("Initial fit failed with error:", e$message, "\n")
    NULL
  })
  
  
  return(fit)
}

fits_fun <- function(mat, ylab = "amplitude (pA)", ids = NULL, attempts = 10, seed = 7) {
  # Run the plotting and fitting in lapply with mouse click pause
  set.seed(seed)
  N <- dim(mat)[2]
  
  # Open a graphical device
  dev.new(width = 4, height = 4, noRStudioGD = TRUE)
  
  fits <- lapply(1:N, function(ii) {
    y <- c(rev(mat[, ii]))
    attempt <- 0
    fit <- NULL

    while (is.null(fit) && attempt < attempts) {
      start_values <- generate_start_values(y, x)
      lower_bounds <- c(0, 0, 0)
      upper_bounds <- c(1, 10 * max(y), 1)
      
      fit <- fit_fun(y, x, start_values, lower_bounds, upper_bounds)
      attempt <- attempt + 1
    }

    if (!is.null(fit)) {
      plot(x, y, main = "", xlab = "input (intensity)", ylab = ylab, 
           xlim = c(0, 1), ylim = c(0, max(y) * 1.1), bty = 'n')
      curve(exp_fun(x, coef(fit)[1], coef(fit)[2], coef(fit)[3]), add = TRUE, 
            col = "slateblue", lty = 3, lwd = 2)
    } else {
      plot(x, y, main = "fit fails", xlab = "input (intensity)", ylab = ylab, 
           xlim = c(0, 1), ylim = c(0, max(y) * 1.1), bty = 'n')
    }

    wait_for_click()

    fit
  })
  
  # Clear the final plot
  dev.off()
  names(fits) <- ids

  # Print fit summaries
  lapply(fits, function(fit) {
    if (!is.null(fit)) {
      summary(fit)
    } else {
      cat("fit fails\n")
    }
  })

  # Extract coefficients
  coefficients <- do.call(rbind, lapply(fits, function(fit) {
    if (!is.null(fit)) {
      as.data.frame(t(coef(fit)))
    } else {
      data.frame(xmin = NA, Imax = NA, tau = NA)
    }
  }))

  return(coefficients)
}

exp_fun <- function(x, xmin, Imax, tau)
  Imax * (1 - exp(-(x-xmin)/tau))



IO_plot <- function(boxplot_values, main='', xlab='relative LED intensity', ylab='relative amplitude', col='indianred', 
  xlim=c(0, 1.05), ylim=c(0, 1.05), x_tick_interval=0.2, y_tick_interval=0.2, tick_length=-0.2, xround=1, yround=1, 
  seed=42, fit.attempts=10, maxiter=500, pch=20, lwd=1.2, height=4, width=4, filename='IO.svg', save=FALSE){
  
  y <- abs(boxplot_values$Median)
  x <- boxplot_values$x

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(x, y, main=main, xlab=xlab, ylab=ylab, col=col, pch=pch, xlim=xlim, ylim=ylim, lwd=lwd, xaxt='n', yaxt='n', bty='n')

  # Fit the exponential to median values
  lower_bounds <- c(0, 0, 0)  
  upper_bounds <- c(1, 2*max(y), 10) 

  set.seed(seed)
  
  attempt <- 0
  fit <- NULL
  while (is.null(fit) && attempt < fit.attempts) {
    start_values <- generate_start_values(y, x)
    fit <- nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit model to the Q1 values
  y_q1 <- abs(boxplot_values$Q1)
  upper_bounds <- c(1, 2*max(y_q1), 10) 
  
  attempt <- 0
  fit1 <- NULL
  while (is.null(fit1) && attempt < fit.attempts) {
    start_values <- generate_start_values(y_q1, x)
    fit1 <- nlsLM(y_q1 ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit model to the Q3 values
  y_q3 <- abs(boxplot_values$Q3)
  upper_bounds <- c(1, 2*max(y_q3), 10) 
  
  attempt <- 0
  fit3 <- NULL
  while (is.null(fit3) && attempt < fit.attempts) {
    start_values <- generate_start_values(y_q3, x)
    fit3 <- nlsLM(y_q3 ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Filter x-values for which the y-values of the fitted curve are positive
  x_curve <- seq(0, max(x), length.out = 200)
  
  y_fit <- exp_fun(x_curve, coef(fit)[1], coef(fit)[2], coef(fit)[3])
  x_curve_pos <- x_curve[y_fit > 0]
  y_fit_pos <- y_fit[y_fit > 0]
  lines(x_curve_pos, y_fit_pos, col = col, lty = 1, lwd = lwd)
  
  y_fit1 <- exp_fun(x_curve, coef(fit1)[1], coef(fit1)[2], coef(fit1)[3])
  x_curve_pos1 <- x_curve[y_fit1 > 0]
  y_fit_pos1 <- y_fit1[y_fit1 > 0]
  lines(x_curve_pos1, y_fit_pos1, col = "gray", lty = 3, lwd = lwd)

  y_fit3 <- exp_fun(x_curve, coef(fit3)[1], coef(fit3)[2], coef(fit3)[3])
  x_curve_pos3 <- x_curve[y_fit3 > 0]
  y_fit_pos3 <- y_fit3[y_fit3 > 0]
  lines(x_curve_pos3, y_fit_pos3, col = "gray", lty = 3, lwd = lwd)

  # Define tick intervals and lengths
  x_ticks <- seq(xround*floor(xlim[1]/xround), xround*ceiling(xlim[2]/xround), by=x_tick_interval)
  y_ticks <- seq(yround*floor(ylim[1]/yround), yround*ceiling(ylim[2]/yround), by=y_tick_interval)
  
  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }

}


IO_plot2 <- function(boxplot_values, main='', xlab='relative LED intensity', ylab='relative amplitude', col='indianred', 
  xlim=c(0, 1.05), ylim=c(0, 1.05), x_tick_interval=0.2, y_tick_interval=0.2, tick_length=-0.2, xround=1, yround=1, 
  seed=42, fit.attempts=10, maxiter=500, pch=20, lwd=1.2, height=4, width=4, spar=0.5, filename='IO.svg', save=FALSE){

  y <- abs(boxplot_values$Median)
  x <- boxplot_values$x

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(x, y, main=main, xlab=xlab, ylab=ylab, col=col, pch=pch, xlim=xlim, ylim=ylim, lwd=lwd, xaxt='n', yaxt='n', bty='n')

  # Fit the exponential to median values
  lower_bounds <- c(0, 0, 0)  
  upper_bounds <- c(1, 2*max(y), 1) 

  set.seed(seed)
  
  attempt <- 0
  fit <- NULL
  while (is.null(fit) && attempt < fit.attempts) {
    start_values <- generate_start_values(y, x)
    fit <- nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit spline to Q1 values
  y_q1 <- abs(boxplot_values$Q1)
  spline_fit1 <- smooth.spline(x, y_q1, spar=spar)

  # Fit spline to Q3 values
  y_q3 <- abs(boxplot_values$Q3)
  spline_fit3 <- smooth.spline(x, y_q3, spar=spar)

  # Filter x-values for which the y-values of the fitted exponential curve are positive
  x_curve <- seq(0, max(x), length.out = 200)
  
  y_fit <- exp_fun(x_curve, coef(fit)[1], coef(fit)[2], coef(fit)[3])
  x_curve_pos <- x_curve[y_fit > 0]
  y_fit_pos <- y_fit[y_fit > 0]
  lines(x_curve_pos, y_fit_pos, col = col, lty = 1, lwd = lwd)

  # Filter and plot the spline for Q1 where the y-values are positive
  spline_fit1$y <- spline_fit1$y[spline_fit1$y > 0]
  spline_fit1$x <- spline_fit1$x[1:length(spline_fit1$y)]
  lines(spline_fit1, col="gray", lty=3, lwd=lwd)

  # Filter and plot the spline for Q3 where the y-values are positive
  spline_fit3$y <- spline_fit3$y[spline_fit3$y > 0]
  spline_fit3$x <- spline_fit3$x[1:length(spline_fit3$y)]
  lines(spline_fit3, col="gray", lty=3, lwd=lwd)

  # Define tick intervals and lengths
  x_ticks <- seq(xround*floor(xlim[1]/xround), xround*ceiling(xlim[2]/xround), by=x_tick_interval)
  y_ticks <- seq(yround*floor(ylim[1]/yround), yround*ceiling(ylim[2]/yround), by=y_tick_interval)
  
  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }
}

datasets2list_2 <- function(names, idxs, id='_fits'){
  all_objects <- ls(envir = .GlobalEnv)
  objects <- all_objects[grep(paste0(id, "$"), all_objects)]
  out_list <- list()
  # Loop through each name and id to group the objects
  for (name in names) {
    for (idx in idxs) {
      # Construct the object name pattern
      pattern <- paste0(name, idx, id)
      
      # Find the object that matches the pattern
      matched_object <- objects[objects == pattern]
      
      # If a match is found, add it to the grouped_fits list
      if (length(matched_object) > 0) {
        if (is.null(out_list[[name]])) {
          out_list[[name]] <- list()
        }
        out_list[[name]][[idx]] <- get(matched_object)
      }
    }
  }
  return(out_list)
}



traces_smoothfits <- function(y, fits, dt=0.1, N=1, IEI=50, stimulation_time=150, baseline=50, func=product1N, filter=FALSE, fc=1000, upsample.fit = c(upsample=TRUE, factor=100)){
  
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  upsample <- upsample.fit[['upsample']]
  factor <- upsample.fit[['factor']]

  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

  xfit <- if (upsample) seq(0, dx * (length(yorig) - 1), by = dx/factor) else xorig


  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces1 <- data.frame(x=xorig, y=yorig, yfilter=yfilter)
  traces2 <- data.frame(xfit=xfit)


  if (identical(func, product1N)){
    fits[N+3] <- fits[N+3] + baseline
  } else if (identical(func, product2N)){
    fits[N+3] <- fits[N+3] + baseline; fits[2*N+6] <- fits[2*N+6] + baseline
  } else if (identical(func, product3N)){
    fits[N+3] <- fits[N+3] + baseline; fits[2*N+6] <- fits[2*N+6] + baseline; fits[3*N+9] <- fits[3*N+9] + baseline
  }    
  traces2$yfit <- func(fits,traces2$x+dx, N=N, IEI=IEI)  
  if (identical(func, product2N)){
    traces2$yfit1 <- product1N(fits[1:(N+3)],traces2$x+dx, N=N, IEI=IEI) 
    traces2$yfit2 <- product1N(fits[(N+4):(2*N+6)],traces2$x+dx, N=N, IEI=IEI) 
  } 
  if (identical(func, product3N)){
    traces2$yfit1 <- product1N(fits[1:(N+3)],traces2$x+dx, N=N, IEI=IEI) 
    traces2$yfit2 <- product1N(fits[(N+4):(2*N+6)],traces2$x+dx, N=N, IEI=IEI)
    traces2$yfit3 <- product1N(fits[(2*N+7):(3*N+9)],traces2$x+dx, N=N, IEI=IEI) 
  } 

  return(list(original=traces1, fit=traces2))
}



fit_plot5 <- function(traces, func=product1N, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=4, height=4, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$original$x, traces$original$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main='')
  
  if (filter) {
    lines(traces$original$x, traces$original$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$fit$xfit, traces$fit$yfit, col='#CD5C5C', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2N)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3N)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}


smooth.plots <- function(y, fits, N=1, IEI=50, dt=0.1,  stimulation_time=150, baseline=50, func=product1N, filter=FALSE, fc=1000, upsample.fit = c(upsample=FALSE, factor=100),
  xlab='time (ms)', ylab='', xlim=NULL, ylim=NULL, lwd=1.2, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE){

  traces <- traces_smoothfits(y=y, fits=fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, upsample.fit = upsample.fit)

  fit_plot5(traces=traces, func=func, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, lwd=lwd, filter=filter, width=width, height=height, bg=bg, filename=filename, save=save) 
  
}

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to identify points that do not belong in a cluster
# eps defines the maximum distance between two points for one to be considered in the neighborhood of the other.
# As such, it directly affects how clusters are formed.If too small, many points will be labeled as noise, and the 
# algorithm might fail to form meaningful clusters. If too large, distinct clusters may merge into one big cluster.

# An appropriate eps value is based on the data distribution. Examine the k-Nearest Neighbors (k-NN) Distance Plot
# Define minPts as the minimum points in a cluster. The value k = minPts - 1. Then plot the k-th nearest neighbor 
# distances for all data points in ascending order. Look for the “elbow” in the curve. 
# The distance value at this point is a good candidate for eps.

# DBSCAN is an algorithm that discovers clusters by identifying regions of high density, which are separated by sequentia of lower density. 
# It is effective for finding clusters of varying shapes and sizes in large datasets, even when noise or outliers are present.
# The algorithm operates using two key parameters:
# minPts: The minimum number of points required for a region to be considered a cluster.
# epsilon (eps): Defines the neighborhood radius around a point for identifying neighboring points, often measured using Euclidean distance.

# k-NN distance: In general, for a given value of k, the k-NN distance is the distance from a point to its k-th nearest neighbor

# The kNN distance plot helps identify the “elbow” point, which corresponds to a good eps value
# Use k = minPts - 1; minPts <- 4  is default for 2D data
# Small minPts values (e.g., 3 or 4) work well for data with:
#   Small clusters
#   High density
# Large minPts values (e.g., 10 or 15) work well for:
#   Large datasets
#   Clusters with high variability in density

# dbscan::dbscan performs the following operations
#   1.  Core Points Identification:
#   • For each point in the dataset:
#   • Compute the number of points within a radius of eps (including the point itself).
#   • If this count is ≥ minPts, the point is classified as a core point.
#   2.  Cluster Formation:
#   • Starting from a core point, expand the cluster by adding all points that are directly reachable from the core point:
#   • A point is directly reachable if it lies within the eps radius of the core point.
#   • If this point is also a core point, repeat the process to find more directly reachable points (density reachability).
#   3.  Noise Points Identification:
#   • Any point that:
#   • Is not a core point.
#   • Is not directly reachable from any core point.
#   • These points are classified as noise (or outliers) and assigned a cluster label of 0.


kNNdistplot2 <- function(x, k, minPts,  bty="n", lwd=1, lty=1, axes=FALSE, frame=FALSE, xtick=1, ...) {
    
    if (missing(k) && missing(minPts)) 
        stop("k or minPts need to be specified.")
    if (missing(k)) 
        k <- minPts - 1
    
    # Sort the k-NN distances
    kNNdist <- sort(dbscan::kNNdist(x, k, ...))
    factor <- 10^floor(log10(max(kNNdist)))  # Dynamically determine rounding factor
    ylim <- c(0, factor * ceiling(max(kNNdist) / factor))
    # Dynamically calculate tick intervals based on the unified limit
    yticks <- pretty(ylim, n = 5)
    xlim <- c(1, length(kNNdist))
    # Plot with additional parameters passed using ...
    plot(kNNdist, type="l", ylab=paste0(k, "-NN distance"), 
         xlab="points sorted by distance", xlim=xlim, ylim=ylim, 
         bty=bty, lwd=lwd, lty=lty, axes=axes, frame=frame, ...)
 
  if (!axes){
    axis(1, at=seq(xlim[1], xlim[2], by=xtick), tcl=-0.2, las = 1)  
    axis(2, at=yticks, tcl=-0.2, las = 1)  
  }

}

DBSCAN_analyse <- function(data, minPts = 4, k = NULL, height = 5, width = 10, bg = 'transparent', 
                           eps = NA, filename = 'DBscan.svg', save = FALSE) {
  # Ensure k is derived from minPts if k is not provided
  if (is.null(k)) {
    k <- minPts - 1
  }
  
  # Open a wider plotting window
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  
  # Split the plotting region into two panels
  par(mfrow = c(1, 2))  # 1 row, 2 columns
  
  # Panel 1: kNN plot
  kNNdistplot2(data, k = k)
  
  proceed <- if (is.na(eps)) FALSE else TRUE
  
  if (proceed) {
    # Draw a horizontal line at the specified eps
    abline(h = eps, col = 'indianred', lty = 3)
  } else {
    while (!proceed) {
      # Prompt for eps value
      while (is.na(eps)) {
        cat('\nEnter eps: ')
        eps <- as.numeric(readLines(n = 1))
        if (is.na(eps)) {
          cat('\nInvalid input. Please enter a numeric value.\n')
        }
      }

      # Redraw the kNN plot in the left panel
      par(mfg = c(1, 1))  # Focus back on the first panel
      kNNdistplot2(data, k = k)
      abline(h = eps, col = 'indianred', lty = 3)

      # Ask user if they're happy with the eps
      cat('\nAre you happy with the position of the line (y/n)? ')
      response <- tolower(readLines(n = 1))

      if (response == 'y') {
        proceed <- TRUE
      } else {
        eps <- NA
        cat('\nTry again...\n')
      }
    }
  }

  # Apply DBSCAN on the dataset
  dbscan_result <- dbscan::dbscan(x=data, eps=eps, minPts=minPts)
  
  # Panel 2: DBSCAN plot
  par(mfg = c(1, 2))  # Switch to the second panel
  DBSCAN_plot(data, dbscan_result)

  # Save the plot to an SVG file if save = TRUE
  if (save) {
    # Open SVG device
    svg(file = filename, width = width, height = height, bg = bg)
    # Recreate the plots
    par(mfrow = c(1, 2))
    # Replot kNN plot
    kNNdistplot2(data, k = k)
    abline(h = eps, col = 'indianred', lty = 3)
    # Replot DBSCAN plot
    par(mfg = c(1, 2))
    DBSCAN_plot(data, dbscan_result)
    dev.off()  # Close SVG device
  }
}

DBSCAN_plot <- function(data, dbscan_result) {
  # Ensure `data` is a matrix for proper subsetting
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # Extract column names for axis labels
  xlab <- colnames(data)[1]
  ylab <- colnames(data)[2]
  
  # Determine the maximum value across both axes
  # max_limit <- 100 * ceiling(max(data) / 100)
  factor <- 10^floor(log10(max(data))) 
  max_limit <- factor * ceiling(max(data) / factor)

  lim <- c(0, max_limit)
  
  # Dynamically calculate tick intervals based on the unified limit
  ticks <- pretty(lim, n = 5)
  
  # Plot all data points
  plot(data, col = 'black', pch = 19, cex = 0.75, main = 'DBSCAN clustering results',
       xlim = lim, ylim = lim, bty = 'n', lwd = 1, lty = 1, axes = FALSE, frame = FALSE,
       xlab = xlab, ylab = ylab)
  
  # Highlight noise points (cluster == 0)
  noise_indices <- which(dbscan_result$cluster == 0)
  if (length(noise_indices) > 0) {
    noise_points <- data[noise_indices, , drop = FALSE]
    points(noise_points, col = 'indianred', pch = 19)
    text(noise_points, labels = noise_indices, pos = 4, col = 'indianred', cex = 0.75)
  }
  
  # Add axes with adaptive ticks
  axis(1, at = ticks, tcl = -0.2, las = 1)  # x-axis
  axis(2, at = ticks, tcl = -0.2, las = 1)  # y-axis
}

# DBSCAN_plot <- function(data, dbscan_result, height = height, width = width) {
#   dev.new(width = width, height = height, noRStudioGD = TRUE)
  
#   # Ensure `data` is a matrix for proper subsetting
#   if (!is.matrix(data)) {
#     data <- as.matrix(data)
#   }
  
#   # Extract column names for axis labels
#   xlab <- colnames(data)[1]
#   ylab <- colnames(data)[2]
  
#   # Determine the maximum value across both axes
#   max_limit <- 100 * ceiling(max(data) / 100)
#   lim <- c(0, max_limit)
  
#   # Dynamically calculate tick intervals based on the unified limit
#   ticks <- pretty(lim, n = 5)
  
#   # Plot all data points
#   plot(data, col = 'black', pch = 19, cex = 0.75, main = '', xlim = lim, ylim = lim, 
#        bty = 'n', lwd = 1.2, lty = 1, axes = FALSE, frame = FALSE, 
#        xlab = xlab, ylab = ylab)
  
#   # Identify noise points (cluster == 0)
#   noise_indices <- which(dbscan_result$cluster == 0)
#   if (length(noise_indices) > 0) {
#     noise_points <- data[noise_indices, , drop = FALSE]
    
#     # Highlight noise points
#     points(noise_points, col = "indianred", pch = 19)
    
#     # Add labels to the right of noise points using row identifiers
#     text(noise_points, labels = noise_indices, pos = 4, col = "indianred", cex = 0.75)
#   }
  
#   # Add axes with adaptive ticks
#   axis(1, at = ticks, tcl = -0.2)  # x-axis
#   axis(2, at = ticks, tcl = -0.2)  # y-axis
# }

# DBSCAN_analyse <- function(data, height=5, width=5) {

#   # Use kNNdistplot to select eps
#   dev.new(width=width, height=height, noRStudioGD=TRUE)
#   kNNdistplot2(data, k=2)

#   eps <- NA
#   proceed <- FALSE

#   # Loop until the user is happy with the abline position
#   while (!proceed) {
#     # Prompt for eps value
#     while (is.na(eps)) {
#       cat('\nEnter eps: ')
#       eps <- as.numeric(readLines(n=1))
#       if (is.na(eps)) {
#         cat('\nInvalid input. Please enter a numeric value.\n')
#       }
#     }

#     # Get current x-axis limits
#     x_limits <- par("usr")[1:2]  # This gets the x-axis limits from the plot (xmin and xmax)

#     # Draw the abline with the current eps value, restricted to x-axis limits
#     segments(x0=x_limits[1], y0=eps, x1=x_limits[2], y1=eps, col="indianred", lty=3)

#     # Ask user if they're happy with the line
#     cat('\nAre you happy with the position of the line (y/n)? ')
#     response <- tolower(readLines(n=1))

#     # Check if user is happy
#     if (response == 'y') {
#       proceed <- TRUE
#     } else {
#       # Reset eps and prompt again
#       eps <- NA
#       cat('\ntry again...\n')
#     }
#   }

#   dev.off()

#   # Apply DBSCAN on the dataset
#   dbscan_result <- dbscan::dbscan(data, eps=eps, minPts=3)

#   DBSCAN_plot(data, dbscan_result, height=height, width=width)
# }

# Function to remove loaded objects
remove_loaded_objects <- function(loaded_objects) {
  for (dataset in names(loaded_objects)) {
    objects_to_remove <- loaded_objects[[dataset]]
    rm(list = objects_to_remove, envir = .GlobalEnv)
  }
  # Remove the loaded_objects list itself
  rm(loaded_objects, envir = .GlobalEnv)
}

create_test_output <- function(parameter, test_result) {
  # Create an empty data frame for mixed types
  output_matrix <- data.frame(
    parameter = parameter,  # Add the parameter name
    test = as.character(test_result$method),
    alternative = as.character(test_result$alternative),
    W = as.numeric(test_result$statistic),
    p.value = as.numeric(test_result$p.value),
    stringsAsFactors = FALSE  # Avoid factors for character columns
  )
  
  return(output_matrix)
}

create_art_output <- function(formula, data, parameter=NULL, dp=5) {
  
  model <- NULL  # Clear any old model

  temp_model <- tryCatch(
    {
      # Use `try(..., silent = TRUE)` to suppress messages
      res <- try(art(formula = formula, data = data), silent = TRUE)
      
      # 2) Check if `res` is a `try-error`; if so, raise an error for `tryCatch` to handle
      if (inherits(res, "try-error")) {
        # Extract the error message from the `try-error` object
        error_msg <- attr(res, "condition")$message
        stop(error_msg)
      }
      
      # If successful, `res` is the fitted 'art' model
      res
    },
    error = function(e) {
      e
    }
  )

  if (!inherits(temp_model, "error")) {
    # Assign the successful model to 'model'
    model <- temp_model
    anova_res <- anova(model)
  }else{
    # Return a single-row data frame with the error
    return(data.frame(
      parameter = parameter,
      test      = 'Analysis of Variance of Aligned Rank Transformed Data', 
      factor    = '',
      df        = '',
      dfe       = '',
      `F value` = '',
      `Pr(>F)`  = '',
      Note      = gsub("\\.$", "", error_msg),   # Store the error text here
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  
  anova_res <- anova(model)
  
  # Figure out if it's an anova.art object or a regular anova
  if ('anova.art' %in% class(anova_res)) {
    header_text <- 'Analysis of Variance of Aligned Rank Transformed Data'
  } else {
    header_text <- 'ANOVA'
  }
  
  # Convert the ANOVA result to a data frame
  anova_df <- as.data.frame(anova_res)
  factor_names <- rownames(anova_df)
  
  # Build the output data frame with multiple rows for each factor
  # and only the first row containing the parameter & test name
  output_matrix <- data.frame(
    parameter = c(parameter, rep("", nrow(anova_df) - 1)),  # Only first row has 'parameter'
    test      = c(header_text, rep("", nrow(anova_df) - 1)), # Only first row has 'test'
    factor    = factor_names,                                # e.g., celltype, condition, etc.
    df        = round(anova_df$Df, digits=dp),               # Degrees of freedom
    dfe       = round(anova_df$Df.res, digits=dp),           # Residual degrees of freedom
    `F value` = round(anova_df$`F value`, digits=dp),        # F statistic
    `Pr(>F)`  = round(anova_df$`Pr(>F)`, digits=dp),         # p-value
    Note      = "",                                          # Blank note for successful runs
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(output_matrix)
}

# if ms and pA then output would be fC so 1e3 corrects to pC
trap_fun <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / 1e3
}

fit_plot2 <- function(traces, func=product2, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, main='', bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$x, traces$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main=main)
  lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}

# single_egs <- function(x, y, sign=-1, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, height=4, width=2.5, xbar=100, ybar=50,  color='#4C77BB', filename='trace1.svg', bg='transparent', save=FALSE) {
  
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   if (is.null(ylim)) ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)

#   if (is.null(xlim)) xlim <- c(min(x), max(x))
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))

#   plot(x[idx1:idx2], y[idx1:idx2], type='l', col=color, xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

#   #  scale bar lengths and ybar position
#   ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
  
#   # Add scale bars at the bottom right
#   x_start <- max(xlim) - xbar - 50
#   y_start <- ybar_start
#   x_end <- x_start + xbar
#   y_end <- y_start + ybar
  
#   # Draw the scale bars
#   segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black')
#   segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black')
  
#   # Add labels to the scale bars
#   if (show_text) {
#     text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
#     text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  labels = paste(ybar, 'pA'), adj = c(0.5, 0.5), srt = 90)
#   }
    
#   if (save) {
#     dev.off()
#   }
# }

# assumes stimulation occurs at end of baseline
# use tmax to limit fit (with very small but clear responses, sometimes trace can return to baseline BUT fluctuations in trace make measurement inaccurate)
charge_transfer_fun <- function(x, y, tmax=NULL, fc=300, dt=0.1, baseline=40, filter=TRUE, width=5, height=5, bg='transparent', filename='trace.svg', showplot=FALSE, save=FALSE){

  idx1 <- baseline / dt
  y <- y - mean(y[1:idx1])

  if (filter) {
    fs=1 / dt * 1000; bf <- butter(2, fc / (fs / 2), type='low')
    yfilter <- as.numeric(signal::filter(bf, y))
  } else {
    ind=1
    yfilter=y
  }

  if (showplot){
    traces <- data.frame(x=x, y=y, yfilter=yfilter)
    fit_plot2(traces, width=width, height=height, bg=bg, filename=filename, save=save)
  }

  idx2 <- if (is.null(tmax)) length(x) else tmax / dt

  trap_fun(x[idx1:idx2], yfilter[idx1:idx2])

}

charge_fun <- function(data_list, condition='Control', fc=300, dt=0.1, tmax=NULL, baseline=10, filter=TRUE, showplot=TRUE) {
  sapply(data_list, function(df) {
    if (condition %in% colnames(df)) {
      x <- df$time
      y <- df[[condition]]  # Select column dynamically
      auc <- charge_transfer_fun(x, y, tmax=tmax, fc=fc, dt=dt, baseline=baseline, filter=filter, showplot=showplot)  # Compute charge transfer
      return(auc)
    } else {
      return(NA)  # Return NA if column is missing
    }
  })
}

# Function to calculate the half-width of y(t) = A * (exp(-t / tau2) - exp(-t / tau1))
half_width <- function(A, tau1, tau2, limit=100) {
  # Define the response function y(t)
  y <- function(t) {
    A * (exp(-t / tau2) - exp(-t / tau1))
  }
  
  # Find the peak value of y(t) and the corresponding time (t_peak)
  opt <- optimize(y, interval = c(0, 100), maximum = TRUE)
  t_peak <- opt$maximum
  y_max <- opt$objective
  
  # Define the target half-maximum value
  y_half_max <- y_max / 2
  
  # Define a function for the difference from half-maximum
  half_max_eq <- function(t) {
    y(t) - y_half_max
  }
  
  # Solve for t1 (before the peak) where y(t) = y_max / 2
  t1 <- uniroot(half_max_eq, interval = c(0, t_peak))$root
  
  # Solve for t2 (after the peak) where y(t) = y_max / 2
  t2 <- uniroot(half_max_eq, interval = c(t_peak, limit))$root
  
  # Calculate the half-width
  half_width <- t2 - t1
  
  # Return results as a numeric vector
  c(t1 = t1, t2 = t2, half_width = half_width)
}


amplifier_gain <- function(dataset = NULL, headstage_gain = 0.5, additional_gain = NULL,
                           AD_range = c(-10, 10), AD_bits = 16,
                           dp = 3, tol = 1e-3, VClamp = TRUE) {
  
  if (!is.null(dataset) && !is.null(additional_gain)) {
    amplifier_gain3(dataset = dataset,
                    headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    tol = tol,
                    dp = dp,
                    VClamp = VClamp)
    
  } else if (!is.null(dataset)) {
    amplifier_gain1(dataset = dataset,
                    headstage_gain = headstage_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    dp = dp,
                    tol = tol,
                    VClamp = VClamp)
    
  } else {
    amplifier_gain2(headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    VClamp = VClamp)
  }
}

dpA_fun <- function(dataset, tol=1e-2){
  sapply(1:dim(dataset)[2], function(ii){
    vec <- diff(sort(unique(dataset[,ii])))
    vec <- vec[vec>tol]
    min(vec)
    }
  )
}

amplifier_gain1 <- function(dataset, headstage_gain=0.5, AD_range=c(-10, 10), AD_bits=16, dp=3, tol=1e-3, VClamp=TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  
  if (VClamp) {
    dpA <- dpA_fun(dataset=dataset, tol=tol)
    recording_range <- dpA * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain mV/pA' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dpA,
      check.names = FALSE
    )
    
  } else {
    dV <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    recording_range <- dV * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain V/V' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dV,
      check.names = FALSE
    )
  }
  
  return(output)
}

amplifier_gain2 <- function(headstage_gain=0.5, additional_gain=20, AD_range=c(-10, 10), AD_bits=16, VClamp=TRUE) {
  
  if (length(additional_gain) == 1 && length(headstage_gain) > 1) {
    additional_gain <- rep(additional_gain, length(headstage_gain))
  }
  
  if (length(headstage_gain) != length(additional_gain)) {
    stop("if 'additional_gain' is a vector, it must be the same length as 'headstage_gain'")
  }
  
  min_A_D <- AD_range[1]
  max_A_D <- AD_range[2]
  digitiser_range <- abs(diff(AD_range))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits
  
  if (VClamp) {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain mV/pA' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain mV/pA' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      check.names = FALSE
    )
  } else {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain V/V' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain V/V' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      check.names = FALSE
    )
  }
  return(output)
}

amplifier_gain3 <- function(dataset, headstage_gain = 0.5, additional_gain = 20,
                            AD_range = c(-10, 10), AD_bits = 16,
                            tol = 1e-3, dp=3, VClamp = TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits  # actual theoretical digitisation
  
  if (VClamp) {
    dpA_actual <- dpA_fun(dataset = dataset, tol = tol)
    n <- dUnit / dpA_actual 
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain mV/pA' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
    
  } else {
    dV_actual <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    n <- dUnit / dV_actual
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain V/V' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
  }
  
  return(output)
}


fit_limit <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, 
                      y_abline=0.1, height=4, width=4, show_plot=FALSE) { 
  # Calculate peak (unused in the remainder but may be important elsewhere)
  peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, 
                   baseline=baseline, smooth=smooth)
  
  ind1 <- as.integer((stimulation_time - baseline)/dt)
  ind2 <- as.integer(stimulation_time/dt)
  y2plot <- y - mean(y[ind1:ind2])
  
  Y <- y2plot[ind1:length(y2plot)]
  X <- seq(0, dt * (length(Y) - 1), by = dt)

  out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
  A_abline <- out[1]
  avg_t.abline <- out[2]
  avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline

  if (show_plot) {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    
    plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty = 3)
    
    # Get the left and bottom boundaries of the plot
    left_axis <- par("usr")[1]
    bottom_axis <- par("usr")[3]
    
    # Draw horizontal dotted line from the left to avg_t.abline at height A_abline
    lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
    # Draw vertical dotted line from avg_t.abline down to the bottom of the plot
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
    # Determine index for labeling and add labels
    ind3 <- as.integer(avg_t.abline/dt)
    text(x = max(X[ind1:ind3]) * 1.05, y = A_abline * 1.2, 
         labels = paste0(y_abline * 100, " %"), pos = 4, cex = 0.6)
    text(x = max(X[ind1:ind3]) * 1.05, y = bottom_axis * 0.95, 
         labels = paste0(avg_t.abline, " ms"), pos = 4, cex = 0.6)
    
    # Prompt user for the range of x to use for nFIT
    x_limit <- NA
    while (is.na(x_limit)) {
      cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
      x_limit <- as.numeric(readLines(n = 1))
      if (is.na(x_limit)) {
        cat('\nInvalid input. Please enter a numeric value.\n')
      }
    }
    dev.off()
  } else {
    x_limit <- avg_t.abline
  }
  
  return(x_limit)
}

smooth_moving_avg <- function(y, n = 5) {
  y_length <- length(y)
  result <- rep(NA, y_length)
  
  for (i in 1:y_length) {
    # Determine the start and end indices for the window
    start_idx <- max(1, i - floor(n / 2))
    end_idx <- min(y_length, i + floor(n / 2))
    
    # Calculate the mean for the current window
    result[i] <- mean(y[start_idx:end_idx], na.rm = TRUE)
  }
  
  return(result)
}

downsample_fun <- function(data, ds) {
  if (is.vector(data)) {
    data[seq(1, length(data), by = ds)]
  } else {
    data[seq(1, nrow(data), by = ds), , drop = FALSE]
  }
}


# determine_tmax2
determine_tmax2 <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, lwd=1.2, cex=0.6,
  tmax=NULL, y_abline=0.1, xbar=50, ybar=50, xbar_lab='ms', ybar_lab='pA') {
  if (is.null(tmax)) {
    # Calculate peak information (assumes peak.fun and abline_fun are defined)
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
    
    ind1 <- as.integer((stimulation_time - baseline) / dt)
    ind2 <- as.integer(stimulation_time / dt)
    y2plot <- y - mean(y[ind1:ind2])
    
    # Prepare data for plotting
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by=dt)
    
    out <- abline_fun(X, Y, N=N, y_abline=y_abline)
    A_abline <- out[1]
    avg_t.abline <- if (is.na(out[2])) max(X) else out[2]
    
    # Draw the main plot (without axes)
    plot(X, Y, col='indianred', type='l', axes=FALSE, xlab='', ylab='', lwd=lwd, main='', bty='n')
    
    # draw ablines and add text
    usr <- par('usr')
    left_axis <- usr[1]
    bottom_axis <- usr[3]
    lines(c(min(X), max(X)), c(0, 0), col='black', lwd=lwd, lty=3)
    lines(c(min(X), avg_t.abline), c(A_abline, A_abline), col='black', lwd=lwd, lty=3)
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col='black', lwd=lwd, lty=3)
    ind3 <- as.integer(avg_t.abline / dt)
    text(x=max(X[ind1:ind3]) * 1.05, y=A_abline * 1.2, labels=paste0(y_abline * 100, ' %'), pos=4, cex=cex, font=2)
    text(x=max(X[ind1:ind3]) * 1.05, y=bottom_axis * 0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=cex, font=2)
    stim_index <- round(baseline / dt) + 1
    if (stim_index > length(X)) stim_index <- length(X)
    points(X[stim_index], Y[stim_index], pch=8, col='black', cex=1)
    x_offset <- 0.02 * diff(range(X))
    text(x=X[stim_index] + x_offset, y=Y[stim_index], labels='stim', pos=4, col='black', cex=cex, font=2)
    x_limit <- avg_t.abline

  } else {
    x_limit <- tmax
  }
  
  x_limit <- x_limit + stimulation_time - baseline
    
  return(x_limit)
}

# fit_plot3
fit_plot3 <- function(traces, func=product2, lwd=1.2, cex=0.6, filter=FALSE, xbar=50, ybar=50, 
  xbar_lab='ms', ybar_lab='pA') {
  plot(traces$x, traces$y, col='gray', type='l', axes=FALSE, xlab='', ylab='',
       bty='n', lwd=lwd)
  if (filter && !is.null(traces$yfilter)) {
    lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  }
  lines(traces$x, traces$yfit, col='indianred', lty=3, lwd=2 * lwd)
  if (identical(func, product2) || identical(func, product2N)) {
    lines(traces$x, traces$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }
  if (identical(func, product3) || identical(func, product3N)) {
    lines(traces$x, traces$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }
  if (!is.null(traces$bl)) {
    abline(v=traces$bl, col='black', lwd=lwd, lty=3)
  }
  
  # scale bars
  usr <- par('usr')
  x_range <- usr[1:2]
  y_range <- usr[3:4]
  ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
  x_start <- x_range[2] - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar
  
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black')
  segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black')
  text(x=(x_start + x_end) / 2, y=y_start - ybar / 20, 
       labels=paste(xbar, xbar_lab), adj=c(0.5, 1), cex=cex)
  text(x=x_start - xbar / 4, y=(y_start + y_end) / 2, 
       labels=paste(ybar, ybar_lab), srt=90, adj=c(0.5, 0.5), cex=cex)
}

# drawPlot2
drawPlot2 <- function(traces, func=product2N, lwd=1.2, cex=1, filter=FALSE, xbar=50, ybar=50, 
                      xbar_lab='ms', ybar_lab='pA') {
  fit_plot3(traces=traces, func=func, lwd=lwd, cex=cex, filter=filter,
            xbar=xbar, ybar=ybar, xbar_lab=xbar_lab, ybar_lab=ybar_lab)
}



# MCwilcox <- function(formula, df, alternative = 'two.sided',
#                      exact = NULL, na_rm_subjects = TRUE, p_adjust = 'holm') {
#   f_str <- deparse(formula)
#   has_error <- grepl('Error', f_str)

#   if (has_error) {
#     err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
#     subject_var <- strsplit(err_part, '/')[[1]][1]
#     subject_var <- gsub('[[:space:]]', '', subject_var)
#     main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
#     main_formula <- as.formula(main_formula_str)
#   } else {
#     subject_var <- NULL
#     main_formula <- formula
#   }

#   response_var <- all.vars(formula(main_formula))[1]
#   predictors <- all.vars(formula(main_formula))[-1]

#   if (na_rm_subjects && !is.null(subject_var)) {
#     df <- df[ !ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
#   }

#   if (length(predictors) < 1) {
#     stop('Formula must contain at least one predictor for comparisons')
#   }

#   paired_var <- predictors[1]
#   unpaired_var <- if (length(predictors) > 1) predictors[2] else NULL
#   results <- list()

#   ### Paired ###
#   if (!is.null(subject_var) && !is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       levels_pair <- levels(df[[paired_var]])
#     } else {
#       levels_pair <- sort(unique(df[[paired_var]]))
#     }

#     for (lev in levels_pair) {
#       subset_df <- df[df[[paired_var]] == lev, ]
#       if (is.factor(subset_df[[unpaired_var]])) {
#         levels_unpair <- levels(subset_df[[unpaired_var]])
#       } else {
#         levels_unpair <- sort(unique(subset_df[[unpaired_var]]))
#       }
#       if (length(levels_unpair) < 2) next

#       for (i in seq_len(length(levels_unpair) - 1)) {
#         lev1 <- levels_unpair[i]
#         lev2 <- levels_unpair[i + 1]
#         d1 <- subset_df[subset_df[[unpaired_var]] == lev1, ]
#         d2 <- subset_df[subset_df[[unpaired_var]] == lev2, ]
#         common_subj <- intersect(d1[[subject_var]], d2[[subject_var]])
#         d1 <- d1[d1[[subject_var]] %in% common_subj, ]
#         d2 <- d2[d2[[subject_var]] %in% common_subj, ]
#         d1 <- d1[order(d1[[subject_var]]), ]
#         d2 <- d2[order(d2[[subject_var]]), ]
#         y1 <- d1[[response_var]]
#         y2 <- d2[[response_var]]

#         if (length(y1) > 0 && length(y1) == length(y2)) {
#           test <- wilcox.test(y1, y2, paired = TRUE, alternative = alternative, exact = exact)
#           stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#           stat_value <- as.numeric(test$statistic)
#           results[[length(results) + 1]] <- data.frame(
#             comparison   = paste('within', paired_var, lev, '(paired)'),
#             contrast     = paste(unpaired_var, ':', lev1, 'vs', lev2),
#             n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
#             test         = test$method,
#             alternative  = test$alternative,
#             `test stat`  = stat_name,
#             stat         = stat_value,
#             `p value`    = test$p.value,
#             family       = 'paired',
#             stringsAsFactors = FALSE,
#             check.names = FALSE
#           )
#         }
#       }
#     }
#   }

#   ### Unpaired across levels of paired_var within unpaired_var
#   if (!is.null(unpaired_var)) {
#     if (is.factor(df[[unpaired_var]])) {
#       levels_unpair_all <- levels(df[[unpaired_var]])
#     } else {
#       levels_unpair_all <- sort(unique(df[[unpaired_var]]))
#     }

#     for (lev in levels_unpair_all) {
#       subset_df <- df[df[[unpaired_var]] == lev, ]
#       if (is.factor(subset_df[[paired_var]])) {
#         groups_pair <- levels(subset_df[[paired_var]])
#       } else {
#         groups_pair <- sort(unique(subset_df[[paired_var]]))
#       }
#       if (length(groups_pair) < 2) next
#       d1 <- subset_df[subset_df[[paired_var]] == groups_pair[1], ]
#       d2 <- subset_df[subset_df[[paired_var]] == groups_pair[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('within', unpaired_var, lev, '(unpaired)'),
#         contrast   = paste(paired_var, ':', groups_pair[1], 'vs', groups_pair[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   ### Unpaired across levels of unpaired_var within paired_var, if no subjects
#   if (is.null(subject_var) && !is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       levels_pair_all <- levels(df[[paired_var]])
#     } else {
#       levels_pair_all <- sort(unique(df[[paired_var]]))
#     }

#     for (lev in levels_pair_all) {
#       subset_df <- df[df[[paired_var]] == lev, ]
#       if (is.factor(subset_df[[unpaired_var]])) {
#         groups_unpair <- levels(subset_df[[unpaired_var]])
#       } else {
#         groups_unpair <- sort(unique(subset_df[[unpaired_var]]))
#       }
#       if (length(groups_unpair) < 2) next
#       d1 <- subset_df[subset_df[[unpaired_var]] == groups_unpair[1], ]
#       d2 <- subset_df[subset_df[[unpaired_var]] == groups_unpair[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('within', paired_var, lev, '(unpaired)'),
#         contrast   = paste(unpaired_var, ':', groups_unpair[1], 'vs', groups_unpair[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   ### Single-predictor unpaired comparison
#   if (is.null(subject_var) && is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       groups <- levels(df[[paired_var]])
#     } else {
#       groups <- sort(unique(df[[paired_var]]))
#     }
#     if (length(groups) >= 2) {
#       d1 <- df[df[[paired_var]] == groups[1], ]
#       d2 <- df[df[[paired_var]] == groups[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('between', paired_var),
#         contrast   = paste(groups[1], 'vs', groups[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   out <- if (length(results) == 1) results[[1]] else do.call(rbind, results)

#   out$`p adjusted` <- NA
#   for (fam in unique(out$family)) {
#     idx <- which(out$family == fam)
#     out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
#   }

#   out$family <- factor(out$family, levels = c("paired", "unpaired"))
#   out <- out[order(
#     out$family,
#     sub("within (\\w+).*", "\\1", out$comparison),
#     suppressWarnings(as.numeric(sub(".*within \\w+ (\\w+) \\(.*", "\\1", out$comparison)))
#   ), ]
#   out$family <- NULL
#   return(out)
# }

# # Control of family-wise error rate (FWER) or false discovery rate (FDR) occurs within each logically grouped family of tests
# # All paired comparisons are adjusted relative to one another
# # All unpaired comparisons are adjusted relative to one another
# MCwilcox <- function(formula, df, alternative = 'two.sided',
#                      exact = NULL, na_rm_subjects = TRUE,
#                      p_adjust = 'holm') {
  
#   f_str     <- deparse(formula)
#   has_error <- grepl('Error', f_str)
  
#   if (has_error) {
#     err_part         <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
#     subject_var      <- strsplit(err_part, '/')[[1]][1]
#     subject_var      <- gsub('[[:space:]]', '', subject_var)
#     main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
#     main_formula     <- as.formula(main_formula_str)
#   } else {
#     subject_var  <- NULL
#     main_formula <- formula
#   }
  
#   response_var <- all.vars(formula(main_formula))[1]
#   predictors   <- all.vars(formula(main_formula))[-1]
  
#   # drop subjects with any missing response, if requested
#   if (na_rm_subjects && !is.null(subject_var)) {
#     df <- df[!ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
#   }
  
#   results <- list()
  
#   # Helper: is this variable within-subject? (does any subject have >1 value?)
#   is_within_subject <- function(var, subject_var, df) {
#     tab <- tapply(df[[var]], df[[subject_var]], function(x) length(unique(x)))
#     any(tab > 1)
#   }
  
#   # Single predictor
#   if (length(predictors) == 1) {
#     wvar <- predictors[1]
#     levs <- if (is.factor(df[[wvar]])) levels(df[[wvar]]) else sort(unique(df[[wvar]]))
#     if (length(levs) < 2) stop("Need at least two levels of ", wvar)
#     for (i in seq_len(length(levs) - 1)) {
#       a <- levs[i]; b <- levs[i + 1]
#       d1 <- df[df[[wvar]] == a, ]
#       d2 <- df[df[[wvar]] == b, ]
#       if (!is.null(subject_var)) {
#         common <- intersect(d1[[subject_var]], d2[[subject_var]])
#         d1 <- d1[d1[[subject_var]] %in% common, ]
#         d2 <- d2[d2[[subject_var]] %in% common, ]
#         d1 <- d1[order(d1[[subject_var]]), ]
#         d2 <- d2[order(d2[[subject_var]]), ]
#         y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
#         paired <- TRUE
#         n_out <- min(sum(!is.na(y1)), sum(!is.na(y2)))
#       } else {
#         y1 <- d1[[response_var]]
#         y2 <- d2[[response_var]]
#         paired <- FALSE
#         n_out <- paste(sum(!is.na(y1)), 'vs', sum(!is.na(y2)))
#       }
#       if (length(y1) > 0 && length(y2) > 0 && (!paired || length(y1) == length(y2))) {
#         test <- wilcox.test(y1, y2, paired = paired,
#                             alternative = alternative, exact = exact)
#         snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#         snt <- as.numeric(test$statistic)
#         results[[length(results) + 1]] <- data.frame(
#           parameter    = response_var,
#           comparison   = paste('within', wvar, if (paired) '(paired)' else '(unpaired)'),
#           contrast     = paste(a, 'vs', b),
#           n            = n_out,
#           test         = test$method,
#           alternative  = test$alternative,
#           `test stat`  = snm,
#           stat         = snt,
#           `p value`    = test$p.value,
#           family       = if (paired) 'paired' else 'unpaired',
#           stringsAsFactors = FALSE,
#           check.names  = FALSE
#         )
#       }
#     }
#     out <- do.call(rbind, results)
#     out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
#     out$family <- NULL
#     return(out)
#   }
  
#   # Two predictors: run BOTH ways
#   if (length(predictors) == 2) {
#     for (i in 1:2) {
#       pv <- predictors[i]
#       uv <- predictors[3 - i]
      
#       lv_p <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
#       for (lev in lv_p) {
#         subdf <- df[df[[pv]] == lev, ]
#         lv_u  <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
#         if (length(lv_u) < 2) next
#         for (j in seq_len(length(lv_u) - 1)) {
#           a <- lv_u[j]; b <- lv_u[j + 1]
#           d1 <- subdf[subdf[[uv]] == a, ]
#           d2 <- subdf[subdf[[uv]] == b, ]
          
#           if (!is.null(subject_var) && is_within_subject(uv, subject_var, subdf)) {
#             # Paired test
#             cm <- intersect(d1[[subject_var]], d2[[subject_var]])
#             d1_ <- d1[d1[[subject_var]] %in% cm, ]
#             d2_ <- d2[d2[[subject_var]] %in% cm, ]
#             d1_ <- d1_[order(d1_[[subject_var]]), ]
#             d2_ <- d2_[order(d2_[[subject_var]]), ]
#             y1 <- d1_[[response_var]]; y2 <- d2_[[response_var]]
#             paired <- TRUE
#             n_out <- min(sum(!is.na(y1)), sum(!is.na(y2)))
#           } else {
#             # Unpaired
#             y1 <- d1[[response_var]]
#             y2 <- d2[[response_var]]
#             paired <- FALSE
#             n_out <- paste(sum(!is.na(y1)), 'vs', sum(!is.na(y2)))
#           }
#           if (length(y1) > 0 && length(y2) > 0 && (!paired || length(y1) == length(y2))) {
#             test <- wilcox.test(y1, y2, paired = paired,
#                                 alternative = alternative, exact = exact)
#             snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#             snt <- as.numeric(test$statistic)
#             results[[length(results) + 1]] <- data.frame(
#               parameter    = response_var,
#               comparison   = paste('within', pv, lev, if (paired) '(paired)' else '(unpaired)'),
#               contrast     = paste(a, 'vs', b),
#               n            = n_out,
#               test         = test$method,
#               alternative  = test$alternative,
#               `test stat`  = snm,
#               stat         = snt,
#               `p value`    = test$p.value,
#               family       = if (paired) 'paired' else 'unpaired',
#               stringsAsFactors = FALSE,
#               check.names  = FALSE
#             )
#           }
#         }
#       }
#     }
    
#     out <- do.call(rbind, results)
#     # Adjust p-values by family
#     out$`p adjusted` <- NA
#     for (fam in unique(out$family)) {
#       idx <- which(out$family == fam)
#       out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
#     }
#     out$family <- NULL
#     rownames(out) <- seq_len(nrow(out))
#     return(out)
#   }
  
#   stop('Formula must contain one or two predictors.')
# }

# Control of family-wise error rate (FWER) or false discovery rate (FDR) occurs within each logically grouped family of tests
# All paired comparisons are adjusted relative to one another
# All unpaired comparisons are adjusted relative to one another

MCwilcox <- function(formula, df, alternative = 'two.sided',
                     exact = NULL, na_rm_subjects = TRUE,
                     p_adjust = 'holm') {
  
  f_str     <- deparse(formula)
  has_error <- grepl('Error', f_str)
  
  if (has_error) {
    err_part         <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var      <- strsplit(err_part, '/')[[1]][1]
    subject_var      <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula     <- as.formula(main_formula_str)
  } else {
    subject_var  <- NULL
    main_formula <- formula
  }
  
  response_var <- all.vars(formula(main_formula))[1]
  predictors   <- all.vars(formula(main_formula))[-1]
  
  # drop subjects with any missing response, if requested
  if (na_rm_subjects && !is.null(subject_var)) {
    df <- df[!ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
  }
  
  results <- list()
  
  # Helper: is variable within-subject? (does any subject have >1 value?)
  is_within_subject <- function(var, subject_var, df) {
    tab <- tapply(df[[var]], df[[subject_var]], function(x) length(unique(x)))
    if (is.null(tab)) return(FALSE)
    if (length(tab) == 0) return(FALSE)
    tab <- tab[!is.na(tab)]
    any(tab > 1)
  }
  
  # Single predictor
  if (length(predictors) == 1) {
    wvar <- predictors[1]
    levs <- if (is.factor(df[[wvar]])) levels(df[[wvar]]) else sort(unique(df[[wvar]]))
    if (length(levs) < 2) stop("Need at least two levels of ", wvar)
    for (i in seq_len(length(levs) - 1)) {
      a <- levs[i]; b <- levs[i + 1]
      d1 <- df[df[[wvar]] == a, ]
      d2 <- df[df[[wvar]] == b, ]
      if (!is.null(subject_var)) {
        common <- intersect(d1[[subject_var]], d2[[subject_var]])
        d1 <- d1[d1[[subject_var]] %in% common, ]
        d2 <- d2[d2[[subject_var]] %in% common, ]
        d1 <- d1[order(d1[[subject_var]]), ]
        d2 <- d2[order(d2[[subject_var]]), ]
        y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
        paired <- TRUE
        n_out <- min(sum(!is.na(y1)), sum(!is.na(y2)))
      } else {
        y1 <- d1[[response_var]]
        y2 <- d2[[response_var]]
        paired <- FALSE
        n_out <- paste(sum(!is.na(y1)), 'vs', sum(!is.na(y2)))
      }
      if (length(y1) > 0 && length(y2) > 0 && (!paired || length(y1) == length(y2))) {
        test <- wilcox.test(y1, y2, paired = paired,
                            alternative = alternative, exact = exact)
        snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
        snt <- as.numeric(test$statistic)
        results[[length(results) + 1]] <- data.frame(
          parameter    = response_var,
          comparison   = paste('within', wvar, if (paired) '(paired)' else '(unpaired)'),
          contrast     = paste(a, 'vs', b),
          n            = n_out,
          test         = test$method,
          alternative  = test$alternative,
          `test stat`  = snm,
          stat         = snt,
          `p value`    = test$p.value,
          family       = if (paired) 'paired' else 'unpaired',
          stringsAsFactors = FALSE,
          check.names  = FALSE
        )
      }
    }
    out <- do.call(rbind, results)
    out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
    out$family <- NULL
    return(out)
  }
  
  # Two predictors: run BOTH ways
  if (length(predictors) == 2) {
    for (i in 1:2) {
      pv <- predictors[i]
      uv <- predictors[3 - i]
      
      lv_p <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
      for (lev in lv_p) {
        subdf <- df[df[[pv]] == lev, ]
        lv_u  <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
        if (length(lv_u) < 2) next
        for (j in seq_len(length(lv_u) - 1)) {
          a <- lv_u[j]; b <- lv_u[j + 1]
          d1 <- subdf[subdf[[uv]] == a, ]
          d2 <- subdf[subdf[[uv]] == b, ]
          
          if (!is.null(subject_var) && is_within_subject(uv, subject_var, subdf)) {
            # Paired test
            cm <- intersect(d1[[subject_var]], d2[[subject_var]])
            d1_ <- d1[d1[[subject_var]] %in% cm, ]
            d2_ <- d2[d2[[subject_var]] %in% cm, ]
            d1_ <- d1_[order(d1_[[subject_var]]), ]
            d2_ <- d2_[order(d2_[[subject_var]]), ]
            y1 <- d1_[[response_var]]; y2 <- d2_[[response_var]]
            paired <- TRUE
            n_out <- min(sum(!is.na(y1)), sum(!is.na(y2)))
          } else {
            # Unpaired
            y1 <- d1[[response_var]]
            y2 <- d2[[response_var]]
            paired <- FALSE
            n_out <- paste(sum(!is.na(y1)), 'vs', sum(!is.na(y2)))
          }
          if (length(y1) > 0 && length(y2) > 0 && (!paired || length(y1) == length(y2))) {
            test <- wilcox.test(y1, y2, paired = paired,
                                alternative = alternative, exact = exact)
            snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
            snt <- as.numeric(test$statistic)
            results[[length(results) + 1]] <- data.frame(
              parameter    = response_var,
              comparison   = paste('within', pv, lev, if (paired) '(paired)' else '(unpaired)'),
              contrast     = paste(a, 'vs', b),
              n            = n_out,
              test         = test$method,
              alternative  = test$alternative,
              `test stat`  = snm,
              stat         = snt,
              `p value`    = test$p.value,
              family       = if (paired) 'paired' else 'unpaired',
              stringsAsFactors = FALSE,
              check.names  = FALSE
            )
          }
        }
      }
    }
    
    out <- do.call(rbind, results)
    # Adjust p-values by family
    out$`p adjusted` <- NA
    for (fam in unique(out$family)) {
      idx <- which(out$family == fam)
      out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
    }
    family.order <- match(out$family, c('unpaired', 'paired'))
    out <- out[order(family.order),,drop=FALSE]
    out$family <- NULL
    rownames(out) <- seq_len(nrow(out))
    return(out)
  }
  
  stop('Formula must contain one or two predictors.')
}


MCtest.summary <- function(out, formula=formula.MCtest, longdata=longdata.MCtest,
    alternative='two.sided', var.method='unequal') {
  test <- as.data.frame(out$test)
  formula.text <- paste(deparse(formula), collapse=' ')
  has.error <- grepl('Error\\(', formula.text)

  if (has.error) {
    error.text <- sub('.*Error\\((.*)\\).*', '\\1', formula.text)
    subject <- trimws(strsplit(error.text, '/', fixed=TRUE)[[1]][1])
  } else {
    error.text <- ''
    subject <- NULL
  }

  formula.variables <- all.vars(formula)
  response <- formula.variables[1]
  predictors <- setdiff(formula.variables[-1], subject)
  repeated <- predictors[vapply(predictors, function(x) grepl(x, error.text, fixed=TRUE), logical(1))]

  comparison <- contrast <- n <- test.method <- character(nrow(test))

  for (i in seq_len(nrow(test))) {
    row.label <- rownames(test)[i]
    data.subset <- longdata

    if (length(predictors)==1) {
      comparison.variable <- predictors[1]
      contrast.text <- row.label
      contrast.levels <- strsplit(contrast.text, ' - ', fixed=TRUE)[[1]]
    } else {
      row.parts <- strsplit(row.label, ': ', fixed=TRUE)[[1]]
      condition.label <- predictors[startsWith(row.parts[1], paste0(predictors, ' '))][1]
      condition.level <- substring(row.parts[1], nchar(condition.label)+2)
      comparison.label <- predictors[startsWith(row.parts[2], paste0(predictors, ' '))][1]
      contrast.text <- substring(row.parts[2], nchar(comparison.label)+2)
      contrast.levels <- strsplit(contrast.text, ' - ', fixed=TRUE)[[1]]
      condition.matches <- predictors[vapply(predictors, function(x) condition.level %in% as.character(longdata[[x]]), logical(1))]
      comparison.matches <- predictors[vapply(predictors, function(x) all(contrast.levels %in% as.character(longdata[[x]])), logical(1))]
      condition.variable <- if (length(condition.matches)==1) condition.matches else condition.label
      comparison.variable <- if (length(comparison.matches)==1) comparison.matches else comparison.label
      data.subset <- data.subset[as.character(data.subset[[condition.variable]])==condition.level,,drop=FALSE]
    }

    paired <- comparison.variable %in% repeated
    comparison[i] <- paste(
      'within',
      if (length(predictors)==1) comparison.variable else paste(condition.variable, condition.level),
      if (paired) '(paired)' else '(unpaired)'
      )
    contrast[i] <- paste(contrast.levels, collapse=' vs ')

    if (paired) {
      subject.1 <- unique(as.character(data.subset[[subject]][
        as.character(data.subset[[comparison.variable]])==contrast.levels[1] & is.finite(data.subset[[response]])
        ]))
      subject.2 <- unique(as.character(data.subset[[subject]][
        as.character(data.subset[[comparison.variable]])==contrast.levels[2] & is.finite(data.subset[[response]])
        ]))
      n[i] <- length(intersect(subject.1, subject.2))
      test.method[i] <- 'Paired t-test'
    } else {
      n.1 <- sum(as.character(data.subset[[comparison.variable]])==contrast.levels[1] & is.finite(data.subset[[response]]))
      n.2 <- sum(as.character(data.subset[[comparison.variable]])==contrast.levels[2] & is.finite(data.subset[[response]]))
      n[i] <- paste(n.1, 'vs', n.2)
      test.method[i] <- if (var.method=='unequal') 'Welch Two Sample t-test' else 'Two Sample t-test'
    }
  }

  test.column <- function(columns) {
    column <- columns[columns %in% names(test)][1]
    if (length(column)) unname(test[[column]]) else rep(NA_real_, nrow(test))
  }

  output <- data.frame(
    parameter=rep(response, nrow(test)),
    comparison=comparison,
    contrast=contrast,
    n=n,
    test=test.method,
    alternative=rep(alternative, nrow(test)),
    'test stat'=rep('t', nrow(test)),
    df=test.column('df'),
    stat=test.column('t'),
    estimate=test.column('diff'),
    se=test.column('se'),
    Cohens.d=test.column(c('Cohens.d', 'Cohens.d.z')),
    Hedges.g=test.column(c('Hedges.g', 'Hedges.g.z')),
    'p value'=test.column('p'),
    'p adjusted'=unname(out$sig.tests[,'p.adj']),
    check.names=FALSE
    )

  rownames(output) <- NULL
  output
}

MCttest <- function(formula, df, alternative='two.sided', na_rm_subjects=TRUE,
  p_adjust='Sidak-Holm', contr=NULL, var.equal=FALSE){

  if (!is.logical(var.equal) || length(var.equal)!=1 || is.na(var.equal)){
    stop('var.equal must be TRUE or FALSE')
  }

  formula.string <- paste(deparse(formula), collapse='')
  has.error <- grepl('Error', formula.string)

  if (has.error){
    error.part <- sub('.*Error\\((.*)\\).*', '\\1', formula.string)
    subject.name <- strsplit(error.part, '/')[[1]][1]
    subject.name <- gsub('[[:space:]]', '', subject.name)
    main.formula.string <- sub('\\+\\s*Error\\(.*\\)', '', formula.string)
    main.formula <- as.formula(main.formula.string)
  } else {
    subject.name <- NULL
    main.formula <- formula
  }

  formula.variables <- all.vars(main.formula)
  response.name <- formula.variables[1]
  predictor.names <- formula.variables[-1]

  if (!(length(predictor.names) %in% c(1,2))){
    stop('formula must contain one or two predictors')
  }

  required.names <- c(response.name, predictor.names, subject.name)

  if (!all(required.names %in% names(df))){
    stop('formula or subject variables are missing from df')
  }

  if (na_rm_subjects && !is.null(subject.name)){
    remove.subject <- ave(is.na(df[[response.name]]), df[[subject.name]], FUN=any)
    df <- df[!remove.subject,,drop=FALSE]
  }

  analysis.data <- df[,required.names,drop=FALSE]
  analysis.data <- analysis.data[complete.cases(analysis.data),,drop=FALSE]

  if (!nrow(analysis.data)){
    stop('no complete observations are available')
  }

  for (predictor in predictor.names){
    analysis.data[[predictor]] <- droplevels(factor(analysis.data[[predictor]]))

    if (nlevels(analysis.data[[predictor]])<2){
      stop('each predictor must contain at least two levels')
    }
  }

  if (!is.null(subject.name)){
    analysis.data[[subject.name]] <- droplevels(factor(analysis.data[[subject.name]]))
  }

  if (length(predictor.names)==1){
    cell <- analysis.data[[predictor.names]]
  } else {
    cell <- do.call(interaction, c(analysis.data[predictor.names],
      list(sep=':', drop=TRUE, lex.order=FALSE)))
  }

  cell <- droplevels(factor(cell))
  cell.levels <- levels(cell)
  cell.grid <- analysis.data[match(cell.levels, as.character(cell)),predictor.names,drop=FALSE]
  rownames(cell.grid) <- cell.levels

  if (is.null(contr)){

    result.contrasts <- list()

    if (length(predictor.names)==1){
      predictor.levels <- levels(analysis.data[[predictor.names]])

      for (i in seq_len(length(predictor.levels)-1)){
        contrast.row <- rep(0, length(cell.levels))
        contrast.row[match(predictor.levels[i], cell.levels)] <- -1
        contrast.row[match(predictor.levels[i+1], cell.levels)] <- 1
        result.contrasts[[length(result.contrasts)+1]] <- contrast.row
      }

    } else {

      predictor.order <- predictor.names

      if (!is.null(subject.name)){
        subject.rows <- split(seq_len(nrow(analysis.data)), analysis.data[[subject.name]])

        varies.within.subject <- vapply(predictor.names, function(predictor){
          any(vapply(subject.rows, function(i){
            length(unique(analysis.data[[predictor]][i]))>1
          }, logical(1)))
        }, logical(1))

        if (sum(varies.within.subject)==1){
          predictor.order <- c(predictor.names[!varies.within.subject],
            predictor.names[varies.within.subject])
        }
      }

      for (target.name in predictor.order){

        conditioning.name <- setdiff(predictor.names, target.name)
        conditioning.levels <- levels(analysis.data[[conditioning.name]])
        target.levels <- levels(analysis.data[[target.name]])

        for (conditioning.level in conditioning.levels){
          for (i in seq_len(length(target.levels)-1)){

            negative.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i]
              ]

            positive.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i+1]
              ]

            if (length(negative.cell)==1 && length(positive.cell)==1){
              contrast.row <- rep(0, length(cell.levels))
              contrast.row[match(negative.cell, cell.levels)] <- -1
              contrast.row[match(positive.cell, cell.levels)] <- 1
              result.contrasts[[length(result.contrasts)+1]] <- contrast.row
            }
          }
        }
      }
    }

    contr <- do.call(rbind, result.contrasts)
    colnames(contr) <- cell.levels

  } else {

    if (is.numeric(contr) && is.null(dim(contr))){
      contr <- matrix(contr, nrow=1)
    }

    if (!is.matrix(contr) || !is.numeric(contr) || !nrow(contr) ||
      anyNA(contr) || any(!is.finite(contr)) || ncol(contr)!=length(cell.levels)){
      stop('contr must be a finite numeric matrix with one column per factor-level combination')
    }

    if (!is.null(colnames(contr))){
      if (!setequal(colnames(contr), cell.levels)){
        stop(paste('contr column names must match:', paste(cell.levels, collapse=', ')))
      }

      contr <- contr[,cell.levels,drop=FALSE]
    } else {
      colnames(contr) <- cell.levels
    }
  }

  pairwise.contrast <- apply(contr, 1, function(x){
    nonzero <- sort(as.numeric(x[x!=0]))
    length(nonzero)==2 && all(nonzero==c(-1,1))
  })

  if (any(!pairwise.contrast)){
    stop('each row of contr must contain one -1, one 1 and zeros elsewhere')
  }

  results <- vector('list', nrow(contr))

  for (i in seq_len(nrow(contr))){

    negative.cell <- cell.levels[contr[i,]==-1]
    positive.cell <- cell.levels[contr[i,]==1]
    keep <- cell %in% c(negative.cell, positive.cell)

    test.data <- analysis.data[keep,,drop=FALSE]

    if (!is.null(subject.name)){
      test.data[[subject.name]] <- droplevels(test.data[[subject.name]])
    }

    test.data$.cell <- factor(as.character(cell[keep]), levels=c(negative.cell, positive.cell))

    negative.data <- test.data[test.data$.cell==negative.cell,,drop=FALSE]
    positive.data <- test.data[test.data$.cell==positive.cell,,drop=FALSE]

    paired <- FALSE

    if (!is.null(subject.name)){
      negative.subjects <- unique(as.character(negative.data[[subject.name]]))
      positive.subjects <- unique(as.character(positive.data[[subject.name]]))
      common.subjects <- intersect(negative.subjects, positive.subjects)

      if (setequal(negative.subjects, positive.subjects)){
        paired <- TRUE
      } else if (length(common.subjects)){
        stop('a contrast cannot contain a mixture of paired and unpaired subjects')
      }
    }

    if (paired){
      pair.table <- table(test.data[[subject.name]], test.data$.cell)

      if (any(pair.table!=1)){
        stop('paired comparisons require one complete observation per subject and contrast level')
      }

      negative.data <- negative.data[order(negative.data[[subject.name]]),,drop=FALSE]
      positive.data <- positive.data[order(positive.data[[subject.name]]),,drop=FALSE]
      n.output <- nrow(positive.data)
    } else {
      if (!is.null(subject.name) &&
        (anyDuplicated(negative.data[[subject.name]]) || anyDuplicated(positive.data[[subject.name]]))){
        stop('unpaired comparisons require one observation per subject and contrast level')
      }

      n.output <- paste(nrow(positive.data), 'vs', nrow(negative.data))
    }

    positive.values <- positive.data[[response.name]]
    negative.values <- negative.data[[response.name]]

    test <- t.test(positive.values, negative.values, paired=paired,
      var.equal=if (paired) TRUE else var.equal,
      alternative=alternative)

    estimate <- if (paired){
      unname(test$estimate)
    } else {
      unname(test$estimate[1]-test$estimate[2])
    }

    statistic <- unname(test$statistic)
    df.test <- unname(test$parameter)
    se <- unname(test$stderr)

    if (paired){
      Cohens.d <- statistic/sqrt(n.output)
      df.effect <- df.test
    } else {
      n.positive <- length(positive.values)
      n.negative <- length(negative.values)
      pooled.sd <- sqrt(
        ((n.positive-1)*var(positive.values)+(n.negative-1)*var(negative.values))/
          (n.positive+n.negative-2)
        )
      Cohens.d <- estimate/pooled.sd
      df.effect <- n.positive+n.negative-2
    }

    j.df <- if (df.effect<=3e2){
      gamma(df.effect/2)/(sqrt(df.effect/2)*gamma((df.effect-1)/2))
    } else {
      1-3/(4*df.effect-1)
    }

    Hedges.g <- Cohens.d*j.df
    n.output <- as.character(n.output)
    differing <- predictor.names[
      vapply(predictor.names, function(x){
        as.character(cell.grid[negative.cell,x])!=as.character(cell.grid[positive.cell,x])
      }, logical(1))
      ]

    same <- setdiff(predictor.names, differing)

    if (length(predictor.names)==1){
      parameter <- response.name
      comparison <- paste('within', predictor.names, if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(as.character(cell.grid[positive.cell,differing]), 'vs',
        as.character(cell.grid[negative.cell,differing]))
    } else if (length(differing)==1){
      parameter <- response.name
      comparison <- paste('within', same, as.character(cell.grid[positive.cell,same]),
        if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(as.character(cell.grid[positive.cell,differing]), 'vs',
        as.character(cell.grid[negative.cell,differing]))
    } else {
      parameter <- response.name
      comparison <- paste('between cells', if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(positive.cell, 'vs', negative.cell)
    }

    results[[i]] <- data.frame(
      parameter=parameter,
      comparison=comparison,
      contrast=contrast.name,
      n=n.output,
      test=test$method,
      alternative=test$alternative,
      'test stat'=names(test$statistic),
      df=df.test,
      stat=statistic,
      estimate=estimate,
      se=se,
      Cohens.d=Cohens.d,
      Hedges.g=Hedges.g,
      'p value'=test$p.value,
      family=if (paired) 'paired' else 'unpaired',
      check.names=FALSE
      )
  }

  stats_summary <- do.call(rbind, results)
  rownames(stats_summary) <- NULL

  stats_summary$'p adjusted' <- NA_real_

  for (family in unique(stats_summary$family)){
    index <- which(stats_summary$family==family)
    adjusted <- p.adjust.ff(stats_summary$'p value'[index], method=p_adjust, alpha=0.05)
    stats_summary$'p adjusted'[index] <- adjusted[,'p.adj']
  }

  stats_summary$family <- NULL

  stats_summary
}

# MCwilcox <- function(formula, df, alternative = 'two.sided',
#                      exact = NULL, na_rm_subjects = TRUE,
#                      p_adjust = 'holm') {
  
#   f_str     <- deparse(formula)
#   has_error <- grepl('Error', f_str)

#   if (has_error) {
#     err_part         <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
#     subject_var      <- strsplit(err_part, '/')[[1]][1]
#     subject_var      <- gsub('[[:space:]]', '', subject_var)
#     main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
#     main_formula     <- as.formula(main_formula_str)
#   } else {
#     subject_var  <- NULL
#     main_formula <- formula
#   }

#   response_var <- all.vars(formula(main_formula))[1]
#   predictors   <- all.vars(formula(main_formula))[-1]

#   # drop subjects with any missing response, if requested
#   if (na_rm_subjects && !is.null(subject_var)) {
#     df <- df[!ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
#   }

#   results <- list()

#   # single predictor + Error(subject) → paired only
#   if (!is.null(subject_var) && length(predictors) == 1) {
#     wvar <- predictors[1]
#     levs <- if (is.factor(df[[wvar]])) levels(df[[wvar]]) else sort(unique(df[[wvar]]))
#     if (length(levs) < 2) {
#       stop("Need at least two levels of ", wvar, " for paired comparisons")
#     }
#     for (i in seq_len(length(levs)-1)) {
#       a <- levs[i]; b <- levs[i+1]
#       d1 <- df[df[[wvar]]==a, ]
#       d2 <- df[df[[wvar]]==b, ]
#       common <- intersect(d1[[subject_var]], d2[[subject_var]])
#       d1 <- d1[d1[[subject_var]] %in% common, ]
#       d2 <- d2[d2[[subject_var]] %in% common, ]
#       d1 <- d1[order(d1[[subject_var]]), ]
#       d2 <- d2[order(d2[[subject_var]]), ]
#       y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
#       if (length(y1)>0 && length(y1)==length(y2)) {
#         test <- wilcox.test(y1, y2, paired = TRUE,
#                             alternative = alternative, exact = exact)
#         snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#         snt <- as.numeric(test$statistic)
#         results[[length(results)+1]] <- data.frame(
#           parameter    = response_var,
#           comparison   = paste('within', wvar, '(paired)'),
#           contrast     = paste(a, 'vs', b),
#           n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
#           test         = test$method,
#           alternative  = test$alternative,
#           `test stat`  = snm,
#           stat         = snt,
#           `p value`    = test$p.value,
#           family       = 'paired',
#           stringsAsFactors = FALSE,
#           check.names  = FALSE
#         )
#       }
#     }
#     out <- do.call(rbind, results)
#     out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
#     out$family <- NULL
#     return(out)
#   }

#   if (is.null(subject_var) && length(predictors) == 1) {
#     uvar <- predictors[1]
#     levs <- if (is.factor(df[[uvar]])) levels(df[[uvar]]) else sort(unique(df[[uvar]]))
#     if (length(levs) < 2) {
#       stop("Need at least two levels of ", uvar, " for unpaired comparisons")
#     }
#     for (i in seq_len(length(levs)-1)) {
#       a <- levs[i]; b <- levs[i+1]
#       d1 <- df[df[[uvar]]==a, response_var]
#       d2 <- df[df[[uvar]]==b, response_var]
#       test <- wilcox.test(d1, d2, paired=FALSE,
#                           alternative = alternative, exact = exact)
#       snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       snt <- as.numeric(test$statistic)
#       results[[length(results)+1]] <- data.frame(
#         parameter    = response_var,
#         comparison   = paste('within', uvar, '(unpaired)'),
#         contrast     = paste(a, 'vs', b),
#         n            = paste(sum(!is.na(d1)), 'vs', sum(!is.na(d2))),
#         test         = test$method,
#         alternative  = test$alternative,
#         `test stat`  = snm,
#         stat         = snt,
#         `p value`    = test$p.value,
#         family       = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names  = FALSE
#       )
#     }
#     out <- do.call(rbind, results)
#     out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
#     out$family <- NULL
#     return(out)
#   }

#   if (length(predictors) < 2) {
#     stop('Formula must contain at least one predictor (single unpaired) or two (for paired+unpaired)')
#   }

#   pv <- predictors[1]; uv <- predictors[2]

#   # paired within pv if Error(subject) present
#   if (!is.null(subject_var)) {
#     lv_p <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
#     for (lev in lv_p) {
#       subdf <- df[df[[pv]]==lev, ]
#       lv_u  <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
#       if (length(lv_u)<2) next
#       for (i in seq_len(length(lv_u)-1)) {
#         a <- lv_u[i]; b <- lv_u[i+1]
#         d1 <- subdf[subdf[[uv]]==a, ]
#         d2 <- subdf[subdf[[uv]]==b, ]
#         cm <- intersect(d1[[subject_var]], d2[[subject_var]])
#         d1 <- d1[d1[[subject_var]] %in% cm, ]
#         d2 <- d2[d2[[subject_var]] %in% cm, ]
#         d1 <- d1[order(d1[[subject_var]]), ]
#         d2 <- d2[order(d2[[subject_var]]), ]
#         y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
#         if (length(y1)>0 && length(y1)==length(y2)) {
#           test <- wilcox.test(y1, y2, paired = TRUE,
#                               alternative = alternative, exact = exact)
#           snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#           snt <- as.numeric(test$statistic)
#           results[[length(results)+1]] <- data.frame(
#             parameter    = response_var,
#             comparison   = paste('within', pv, lev, '(paired)'),
#             contrast     = paste(a, 'vs', b),
#             n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
#             test         = test$method,
#             alternative  = test$alternative,
#             `test stat`  = snm,
#             stat         = snt,
#             `p value`    = test$p.value,
#             family       = 'paired',
#             stringsAsFactors = FALSE,
#             check.names  = FALSE
#           )
#         }
#       }
#     }
#   }

#   # unpaired across UV
#   lv_u_all <- if (is.factor(df[[uv]])) levels(df[[uv]]) else sort(unique(df[[uv]]))
#   for (lev in lv_u_all) {
#     subdf <- df[df[[uv]]==lev, ]
#     gp    <- if (is.factor(subdf[[pv]])) levels(subdf[[pv]]) else sort(unique(subdf[[pv]]))
#     if (length(gp)<2) next
#     d1 <- subdf[subdf[[pv]]==gp[1], ]
#     d2 <- subdf[subdf[[pv]]==gp[2], ]
#     test <- wilcox.test(d1[[response_var]], d2[[response_var]],
#                         paired=FALSE,
#                         alternative = alternative, exact = exact)
#     snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#     snt <- as.numeric(test$statistic)
#     results[[length(results)+1]] <- data.frame(
#       parameter    = response_var,
#       comparison   = paste('within', uv, lev, '(unpaired)'),
#       contrast     = paste(gp[1], 'vs', gp[2]),
#       n            = paste(sum(!is.na(d1[[response_var]])),
#                            'vs',
#                            sum(!is.na(d2[[response_var]]))),
#       test         = test$method,
#       alternative  = test$alternative,
#       `test stat`  = snm,
#       stat         = snt,
#       `p value`    = test$p.value,
#       family       = 'unpaired',
#       stringsAsFactors = FALSE,
#       check.names  = FALSE
#     )
#   }

#   # if no subject_var also do reverse unpaired within pv
#   if (is.null(subject_var)) {
#     lv_p_all <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
#     for (lev in lv_p_all) {
#       subdf <- df[df[[pv]]==lev, ]
#       gu    <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
#       if (length(gu)<2) next
#       d1 <- subdf[subdf[[uv]]==gu[1], ]
#       d2 <- subdf[subdf[[uv]]==gu[2], ]
#       test <- wilcox.test(d1[[response_var]], d2[[response_var]],
#                           paired=FALSE,
#                           alternative = alternative, exact = exact)
#       snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       snt <- as.numeric(test$statistic)
#       results[[length(results)+1]] <- data.frame(
#         parameter    = response_var,
#         comparison   = paste('within', pv, lev, '(unpaired)'),
#         contrast     = paste(gu[1], 'vs', gu[2]),
#         n            = paste(sum(!is.na(d1[[response_var]])),
#                              'vs',
#                              sum(!is.na(d2[[response_var]]))),
#         test         = test$method,
#         alternative  = test$alternative,
#         `test stat`  = snm,
#         stat         = snt,
#         `p value`    = test$p.value,
#         family       = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names  = FALSE
#       )
#     }
#   }

#   # output
#   out <- do.call(rbind, results)
#   out$`p adjusted` <- NA
#   for (fam in unique(out$family)) {
#     idx <- which(out$family == fam)
#     out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
#   }

#   # reorder as before
#   out$family <- factor(out$family, levels = c("paired","unpaired"))
#   out <- out[ order(
#     out$family,
#     sub("within (\\w+).*", "\\1", out$comparison),
#     suppressWarnings(as.numeric(sub(".*within \\w+ (\\w+) \\(.*", "\\1", out$comparison)))
#   ), ]
#   out <- out[ order(out$comparison), ]

#   # drop the family column
#   out$family <- NULL
#   rownames(out) <- seq_len(nrow(out))
#   return(out)
# }


create_test_output <- function(parameter, test_result) {
  stat_name <- names(test_result$statistic)
  
  output_matrix <- data.frame(
    parameter = parameter,
    test = as.character(test_result$method),
    alternative = as.character(test_result$alternative),
    statistic = as.numeric(test_result$statistic),
    'p value' = as.numeric(test_result$p.value),
    stringsAsFactors = FALSE,
    check.names = FALSE    
  )
  
  # Rename the statistic column to W or V as appropriate
  colnames(output_matrix)[colnames(output_matrix) == 'statistic'] <- stat_name
  
  return(output_matrix)
}

create_df <- function(mat, levels = list(), var_name = 'amplitude', start_id = 1) {
  n <- nrow(mat)
  total <- n * ncol(mat)
  
  df <- data.frame(
    s = factor(rep(seq_len(n) + start_id - 1, times = ncol(mat)))
  )
  
  for (name in names(levels)) {
    label <- levels[[name]]
    if (is.character(label) && length(label) > 1) {
      df[[name]] <- factor(
        rep(seq_along(label), each = n),
        levels = seq_along(label),
        labels = label
      )
    } else {
      df[[name]] <- factor(rep(label, total))
    }
  }
  
  df[[var_name]] <- c(mat)
  return(df)
}

sequential_wilcox <- function(formula, data, 
                              alternative = 'two.sided', exact = NULL, na.rm = TRUE, 
                              p_adjust = 'holm', group_names = NULL) {
  
  f_str     <- deparse(formula)
  has_error <- grepl('Error', f_str)

  if (has_error) {
    err_part         <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var      <- strsplit(err_part, '/')[[1]][1]
    subject_var      <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula     <- as.formula(main_formula_str)
    paired           <- TRUE
  } else {
    main_formula <- formula
    paired       <- FALSE
  }

  vars      <- all.vars(main_formula)
  response  <- vars[1]
  group_var <- vars[2]

  levels <- sort(unique(data[[group_var]]))

  if (!is.null(group_names)) {
    if (length(group_names) != length(levels)) {
      stop('Number of group_names must match the number of unique levels in group')
    }
    names(levels) <- group_names
  } else {
    names(levels) <- as.character(levels)
  }

  comparisons <- lapply(seq_along(levels[-length(levels)]), function(ii) {
    group1 <- data[[response]][data[[group_var]] == levels[ii]]
    group2 <- data[[response]][data[[group_var]] == levels[ii + 1]]

    if (paired && !is.null(subject_var)) {
      subj1 <- data[[subject_var]][data[[group_var]] == levels[ii]]
      subj2 <- data[[subject_var]][data[[group_var]] == levels[ii + 1]]
      common <- intersect(subj1, subj2)
      group1 <- group1[subj1 %in% common]
      group2 <- group2[subj2 %in% common]
    }

    test <- wilcox.test(group1, group2, paired = paired, alternative = alternative, exact = exact, na.rm = na.rm)
    
    stat_name  <- if (!is.null(names(test$statistic))) names(test$statistic)
    stat_value <- as.numeric(test$statistic)

    df <- data.frame(
      parameter    = response,
      contrast     = paste(names(levels)[ii], 'vs', names(levels)[ii + 1]),
      test         = test$method,
      alternative  = test$alternative,
      `test stat`  = stat_name,
      stat         = stat_value,
      `p value`    = test$p.value,
      stringsAsFactors = FALSE,
      check.names  = FALSE
    )
    return(df)
  })

  results <- do.call(rbind, comparisons)
  results$`p adjusted` <- p.adjust(results$`p value`, method = p_adjust)

  return(results)
}



sequential_wilcox_orig <- function(data, groups, var, paired=TRUE, 
                              alternative='two.sided', exact=NULL, na.rm=TRUE, 
                              p_adjust = 'holm', group_names=NULL) {
  
  levels <- sort(unique(data[[groups]]))

  if (!is.null(group_names)) {
    if (length(group_names) != length(levels)) {
      stop('number of group_names should match the number of unique groups in group column')
    }
    names(levels) <- group_names
  } else {
    names(levels) <- as.character(levels)
  }

  comparisons <- lapply(seq_along(levels[-length(levels)]), function(ii) {
    group1 <- data[[var]][data[[groups]] == levels[ii]]
    group2 <- data[[var]][data[[groups]] == levels[ii + 1]]
    
    test <- wilcox.test(group1, group2, paired=paired, alternative=alternative, exact=exact, na.rm=na.rm)   
    
    stat_name <- if (!is.null(names(test$statistic))) names(test$statistic)
    stat_value <- as.numeric(test$statistic)

    df <- data.frame(
      parameter = var,
      contrast = paste(names(levels)[ii], 'vs', names(levels)[ii + 1]),
      test = test$method,
      alternative = test$alternative,
      'test stat' = stat_name,
      'stat' = stat_value,
      'p value' = test$p.value,
      stringsAsFactors = FALSE
    )

    return(df)
  })

  results <- do.call(rbind, comparisons)

  if ('p.value' %in% colnames(results)) {
    colnames(results)[colnames(results) == 'p.value'] <- 'p value'
  }

  results$'p adjusted' <- p.adjust(results$'p value', method = p_adjust)

  return(results)
}


save_graph <- function(svg_path, filename='graph1.svg', width=6, height=4, bg="transparent") {
  old_wd <- getwd()
  setwd(svg_path)
  dev.copy(svg, file=filename, width=width, height=height, bg=bg)
  dev.off()
  setwd(old_wd)
}


# 
analyseABFtk <- function() {

  cexVar <- if (Sys.info()["sysname"] == "Darwin") tclVar("0.6") else tclVar("1.0")

  experimentVar  <- tclVar("voltage clamp")
  unitVar        <- tclVar("")
  dataColVar     <- tclVar("")
  dtVar          <- tclVar("")
  ntracesVar     <- tclVar("")

  is.tkwin <- function(widget) {
    tryCatch({
      tclvalue(tkwinfo("exists", widget)) == "1"
    }, error = function(e) FALSE)
  }


  extract_metadata <- function(abf_dataset) {
    list(
      path                  = abf_dataset$path,
      formatVersion         = abf_dataset$formatVersion,
      channelNames          = abf_dataset$channelNames,
      channelUnits          = abf_dataset$channelUnits,
      samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
      header                = abf_dataset$header,
      tags                  = abf_dataset$tags,
      sections              = abf_dataset$sections
    )
  }

  choose_data_column <- function(channelUnits, experiment) {
    if (experiment == 'voltage clamp') {
      idx <- grep('A', channelUnits, ignore.case = TRUE)
    } else if (experiment == 'current clamp') {
      idx <- grep('V', channelUnits, ignore.case = TRUE)
    } else {
      idx <- integer(0)
    }
    if (length(idx) > 0) return(idx[1])
    else return(NA)
  }

  check_consistency <- function(metadata) {
    dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
    traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
    expType <- tclvalue(experimentVar)
    unit_values <- sapply(metadata, function(meta) {
      col_idx <- choose_data_column(meta$channelUnits, expType)
      if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
    })
    dt_good <- (length(unique(dt_values)) == 1)
    traces_good <- (length(unique(traces_values)) == 1)
    unit_good <- (length(unique(unit_values)) == 1)
    if (dt_good && traces_good && unit_good) {
      return('Data is consistent')
    } else {
      error_msgs <- c()
      if (!dt_good) error_msgs <- c(error_msgs, paste('Inconsistent dt values:', paste(dt_values, collapse = ', ')))
      if (!unit_good) error_msgs <- c(error_msgs, paste('Inconsistent Units:', paste(unit_values, collapse = ', ')))
      if (!traces_good) error_msgs <- c(error_msgs, paste('Inconsistent Traces:', paste(traces_values, collapse = ', ')))
      return(paste(error_msgs, collapse = '; '))
    }
  }

  egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
                       xbar = 100, ybar = 50, color = '#4C77BB', show_bar = FALSE, cex = 0.6) {
    if (is.null(ylim))
      ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
    if (is.null(xlim))
      xlim <- c(min(x), max(x))
    idx1 <- which.min(abs(x - xlim[1]))
    idx2 <- which.min(abs(x - xlim[2]))
    plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color, xlim = xlim, ylim = ylim, bty = 'n', 
      lwd = lwd, lty = 1, axes = FALSE, frame = FALSE, xlab = '', ylab = '')
    if (show_bar) {
      ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
      x_start <- max(xlim) - xbar - 50
      y_start <- ybar_start
      x_end <- x_start + xbar
      y_end <- y_start + ybar
      segments(x_start, y_start, x_end, y_start, lwd = lwd, col = 'black')
      segments(x_start, y_start, x_start, y_end, lwd = lwd, col = 'black')
      if (show_text) {
        text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), 
             adj = c(0.5, 1), cex = cex)
        text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = paste(ybar, 'pA'), 
             adj = c(0.5, 0.5), srt = 90, cex = cex)
      }
    }
  }

  load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
    abf_path <- if (is.null(abf_path)) getwd() else abf_path
    setwd(abf_path)
    N <- length(abf_files)
    datasets <- lapply(seq_len(N), function(ii) readABF(abf_files[ii]))
    names(datasets) <- abf_files
    metadata <- lapply(datasets, extract_metadata)
    return(list(datasets = datasets, metadata = metadata))
  }

  download_data <- function() {
    if (is.null(averaged_data) || length(averaged_data) == 0) {
      tkinsert(consoleText, 'end', 'No averaged data available.\n')
      tkyview.moveto(consoleText, 1.0)
      return()
    }

    df <- as.data.frame(do.call(cbind, averaged_data))
    colnames(df) <- as.character(seq_len(length(averaged_data)))

    download_folder <- tclvalue(folderPathVar)
    if (nchar(download_folder) == 0) {
      tkinsert(consoleText, 'end', 'No folder selected for download.\n')
      tkyview.moveto(consoleText, 1.0)
      return()
    }

    file_path <- tclvalue(tkgetSaveFile(
      initialdir = download_folder,
      defaultextension = ".csv",
      initialfile = "average.csv"
    ))

    if (nchar(file_path) == 0) {
      return()
    }

    # only reached if file_path is valid
    write.csv(df, file = file_path, row.names = FALSE)
    tkinsert(consoleText, 'end', paste0('Saved to: ', file_path, '\n'))
    tkyview.moveto(consoleText, 1.0)
  }

  abf_averages <- function(datasets, baseline = 100, stimulation_time = 350, traces2average = NULL, dataCol = 1, ylim = NULL, xlim = NULL, 
    color = 'darkgray', xbar = 100, ybar = 50, width = 5.25, height = 2.75, cex=0.6, save = FALSE, plotIt = TRUE) {
    
    N <- length(datasets)
    sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
    responses <- lapply(seq_len(N), function(iii) {
      sapply(seq_along(datasets[[iii]]$data), function(ii) {
        datasets[[iii]]$data[[ii]][, dataCol]
      })
    })
    names(responses) <- names(datasets)
    baseline2zero <- function(y, dt, stim, baseline) {
      idx_baseline <- round(baseline / dt)
      idx_start    <- round((stim - baseline) / dt) + 1
      y0 <- y - mean(y[1:idx_baseline])
      y0[idx_start:length(y0)]
    }
    responses0 <- lapply(seq_len(N), function(iii) {
      sapply(seq_len(ncol(responses[[iii]])), function(jj) {
        baseline2zero(responses[[iii]][, jj],
                      dt = sampling_intervals[iii],
                      stim = stimulation_time,
                      baseline = baseline)
      })
    })
    names(responses0) <- names(responses)
    responses0_mean <- if(is.null(traces2average)) {
      lapply(seq_len(N), function(iii) apply(responses0[[iii]], 1, mean))
    } else {
      lapply(seq_len(N), function(iii)
        apply(responses0[[iii]][, traces2average[[iii]], drop = FALSE], 1, mean))
    }
    time <- lapply(seq_len(N), function(iii) {
      dt_val <- sampling_intervals[iii]
      stim_time <- stimulation_time
      base_val <- baseline
      seq(from = stim_time - base_val, by = dt_val, length.out = length(responses0_mean[[iii]]))
    })
    if(plotIt){
      par(mfrow = c(1, N))
      show_bar <- rep(FALSE, N)
      if (N > 0) show_bar[N] <- TRUE
      for(ii in seq_len(N)) {
        egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = 'darkgray',
                   show_bar = FALSE, cex=cex, show_text = FALSE)
      }
    }
    return(list(raw_data = responses,
                baseline_corrected_data = responses0,
                baseline_corrected_mean_data = responses0_mean,
                datasets = datasets))
  }

  combine_abf_data <- function(result) {
    master_abf <- list()
    master_abf$data <- list()
    master_abf$source_files <- c()
    master_abf$samplingIntervalInSec <- result$datasets[[1]]$samplingIntervalInSec
    for(i in seq_along(result$datasets)) {
      ds <- result$datasets[[i]]
      n_traces <- length(ds$data)
      master_abf$data <- c(master_abf$data, ds$data)
      master_abf$source_files <- c(master_abf$source_files, rep(names(result$datasets)[i], n_traces))
    }
    return(master_abf)
  }

  smart_axis_limits <- function(vec, n_steps = 5) {
    rng <- range(vec)
    spread <- diff(rng)
    
    # Pick a base that's a nice round number (1, 2, 5, 10, 20, 50, 100, etc.)
    raw_step <- spread / n_steps
    base <- 10^floor(log10(raw_step))
    
    # Refine to nicer step (1, 2, or 5 × 10^n)
    nice_steps <- c(1, 2, 5, 10)
    best_step <- base * nice_steps[which.min(abs(nice_steps * base - raw_step))]
    
    lower <- floor(rng[1] / best_step) * best_step
    upper <- ceiling(rng[2] / best_step) * best_step
    c(lower, upper)
  }

  # global variables
  master_abf <<- NULL      # Will hold either a concatenated master object or the original structure.
  averaged_data <<- NULL   # Will hold the averaged (baseline_corrected_mean) data.
  traces2average <<- list()  # Used in separate mode.
  # For concatenated (master) mode:
  current_trace <<- 1      
  total_traces <<- 0       
  current_group_selected <<- integer(0)  
  groups_list <<- list()  
  # For separate (non-concatenated) mode:
  current_dataset <<- 1  

  # review functions
  # common plot settings
  tk_par_settings <- function() {
    par(mar = c(4, 5, 2, 1) + 0.1,
        mgp = c(2.5, 0.5, 0),
        tcl = -0.2)
  }

  # global plot size settings
if (Sys.info()["sysname"] == "Darwin") {
  graph_width  <<- 500
  graph_height <<- 300
  graph_hscale <<- 1.0
  graph_vscale <<- 1.0
} else {
  graph_width  <<- 750
  graph_height <<- 450
  graph_hscale <<- 1.5
  graph_vscale <<- 1.5
}

#### review_master_recordings: shows each trace with Accept/Reject/Undo buttons below ####
review_master_recordings <- function() {
  if (is.null(master_abf)) {
    tkinsert(consoleText, 'end', "No master ABF data available. Please load data first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  total_traces <<- length(master_abf$data)
  current_trace <<- 1L
  current_group_selected <<- integer(0)
  groups_list <<- list()

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))

  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(reviewFrame, 0, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 1, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 2, weight = 1)
  tkgrid.rowconfigure(reviewFrame,    1, weight = 1)

  infoLabel <<- tklabel(reviewFrame, text = paste('Trace', current_trace, 'of', total_traces))
  tkgrid(infoLabel, row = 0, column = 1, pady = 5)

  plotWrapper <- tkframe(reviewFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  reviewPlot <<- tkrplot(plotWrapper, fun = function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)

    mat <- master_abf$data[[current_trace]]
    dt  <- master_abf$samplingIntervalInSec * 1000
    time <- seq(0, by = dt, length.out = nrow(mat))
    dc <- as.numeric(tclvalue(dataColVar))
    if (is.na(dc) || dc < 1 || dc > ncol(mat)) dc <- 1
    y  <- mat[, dc]
    stim_time <- as.numeric(tclvalue(stimTimeVar))
    stim_y    <- y[which.min(abs(time - stim_time))]

    plot(time, y, type = 'l', col = 'darkgray',
         xlab = 'time (ms)', ylab = tclvalue(unitVar),
         xlim = smart_axis_limits(time),
         ylim = smart_axis_limits(y),
         axes = FALSE, bty = 'l')
    axis(1); axis(2, las = 1)
    points(stim_time, stim_y, pch = 8, col = 'black')
    text(stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
  }, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(reviewPlot, row = 0, column = 0, sticky = 'nsew')

  redraw_console_master <- function() {
    tkdelete(consoleText, '1.0', 'end')
    if (length(current_group_selected) == 0) {
      tkinsert(consoleText, 'end', 'No traces selected.\n')
    } else {
      tkinsert(consoleText, 'end',
        paste0('Selected traces: ', paste(current_group_selected, collapse = ', '), '\n'))
    }
    tkyview.moveto(consoleText, 1.0)
  }

  move_next_master <- function() {
    if (current_trace < total_traces) {
      current_trace <<- current_trace + 1L
      tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
      tkrreplot(reviewPlot)
    } else {
      tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.\n')
      tkconfigure(acceptButton, state = 'disabled')
      tkconfigure(rejectButton, state = 'disabled')
      tkconfigure(undoButton,   state = 'disabled')
      tkyview.moveto(consoleText, 1.0)
    }
  }

  navBar <- tkframe(reviewFrame)
  tkgrid(navBar, row = 2, column = 1, pady = 5)

  acceptButton <<- tkbutton(navBar, text = 'Accept', command = function() {
    current_group_selected <<- c(current_group_selected, current_trace)
    redraw_console_master()
    move_next_master()
  })
  rejectButton <<- tkbutton(navBar, text = 'Reject', command = move_next_master)
  undoButton   <<- tkbutton(navBar, text = 'Undo',   command = function() {
    if (length(current_group_selected) > 0)
      current_group_selected <<- head(current_group_selected, -1)
    if (current_trace > 1) {
      current_trace <<- current_trace - 1L
      tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
      tkrreplot(reviewPlot)
    }
    redraw_console_master()
  })
  averageGroupButton <<- tkbutton(navBar, text = 'Add Selected Group', command = function() {
    if (length(current_group_selected) > 0) {
      groups_list[[length(groups_list) + 1]] <<- current_group_selected
      current_group_selected <<- integer(0)
      redraw_console_master()
    }
  })
  selectionCompleteButton <<- tkbutton(navBar, text = 'Selection Complete', command = function() {
    tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.\n')
    tkyview.moveto(consoleText, 1.0)
  })

  tkgrid(acceptButton, row = 0, column = 0, padx = 5)
  tkgrid(rejectButton, row = 0, column = 1, padx = 5)
  tkgrid(undoButton,   row = 0, column = 2, padx = 5)
  tkgrid(averageGroupButton, row = 1, column = 0, padx = 5, pady=5)
  tkgrid(selectionCompleteButton, row = 1, column = 1, padx = 5, pady=5)
}

review_recordings <- function() {
  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  if (!exists('abf_analysis_result', envir = .GlobalEnv)) {
    tkinsert(consoleText, 'end', "No analysis result available for review.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  result         <- get('abf_analysis_result', envir = .GlobalEnv)
  datasets       <- result$datasets
  traces2average <<- lapply(datasets, function(x) integer(0))
  current_dataset <<- 1L
  current_trace   <<- 1L

  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(reviewFrame, 0, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 1, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 2, weight = 1)
  tkgrid.rowconfigure(reviewFrame,    1, weight = 1)

  infoLabel <<- tklabel(reviewFrame, text = paste(names(datasets)[1], 'trace', 1))
  tkgrid(infoLabel, row = 0, column = 1, pady = 5)

  plotWrapper <- tkframe(reviewFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  reviewPlot <<- tkrplot(plotWrapper, fun = function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)

    ds    <- datasets[[current_dataset]]
    fname <- names(datasets)[current_dataset]
    tkconfigure(infoLabel, text = paste(fname, 'trace', current_trace))

    if (current_trace > length(ds$data)) {
      plot.new()
      text(0.5, 0.5, paste('No more recordings in', fname))
    } else {
      mat <- ds$data[[current_trace]]
      dt  <- ds$samplingIntervalInSec * 1000
      time<- seq(0, by = dt, length.out = nrow(mat))
      dc  <- as.numeric(tclvalue(dataColVar))
      if (is.na(dc) || dc < 1 || dc > ncol(mat)) dc <- 1
      y   <- mat[, dc]
      stim_time <- as.numeric(tclvalue(stimTimeVar))
      stim_y    <- y[which.min(abs(time - stim_time))]

      plot(time, y, type = 'l', col = 'darkgray',
           xlab = 'time (ms)', ylab = tclvalue(unitVar),
           xlim = smart_axis_limits(time),
           ylim = smart_axis_limits(y),
           axes = FALSE, bty = 'l')
      axis(1); axis(2, las = 1)
      points(stim_time, stim_y, pch = 8, col = 'black')
      text(  stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
    }
  }, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(reviewPlot, row = 0, column = 0, sticky = 'nsew')

  redraw_console_recordings <- function() {
    tkdelete(consoleText, '1.0', 'end')
    approved <- traces2average[[current_dataset]]
    msg <- if (length(approved) == 0) {
      paste0('No approved traces for ', names(datasets)[current_dataset])
    } else {
      paste0('Approved traces for ', names(datasets)[current_dataset], ': ',
             paste(approved, collapse = ', '))
    }
    tkinsert(consoleText, 'end', paste0(msg, '\n'))
    tkyview.moveto(consoleText, 1.0)
  }

  move_next_recordings <- function() {
    ds <- datasets[[current_dataset]]
    if (current_trace < length(ds$data)) {
      current_trace <<- current_trace + 1L
      tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
      tkrreplot(reviewPlot)
    } else {
      if (current_dataset < length(datasets)) {
        current_dataset <<- current_dataset + 1L
        current_trace   <<- 1L
        tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
        tkrreplot(reviewPlot)
      } else {
        tkinsert(consoleText, '1.0', 'end')
        tkinsert(consoleText, 'end', 'Review complete: Approved recordings stored.\n')
        tkconfigure(acceptButton, state = 'disabled')
        tkconfigure(rejectButton, state = 'disabled')
        tkconfigure(undoButton,   state = 'disabled')
        tkyview.moveto(consoleText, 1.0)
      }
    }
  }

  navBar <- tkframe(reviewFrame)
  tkgrid(navBar, row = 2, column = 1, pady = 5)
  acceptButton <<- tkbutton(navBar, text = 'Accept', command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    redraw_console_recordings()
    move_next_recordings()
  })
  rejectButton <<- tkbutton(navBar, text = 'Reject', command = move_next_recordings)
  undoButton   <<- tkbutton(navBar, text = 'Undo',   command = function() {
    if (length(traces2average[[current_dataset]]) > 0)
      traces2average[[current_dataset]] <<- head(traces2average[[current_dataset]], -1)
    if (current_trace > 1) {
      current_trace <<- current_trace - 1L
      tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
      try({
        if (exists('reviewPlot', inherits = TRUE)) {
          widget_id <- as.character(reviewPlot$ID)
          if (tcl('winfo', 'exists', widget_id) == '1') tkrreplot(reviewPlot)
        }
      }, silent = TRUE)

    }
    redraw_console_recordings()
  })
  tkgrid(acceptButton, row = 0, column = 0, padx = 5)
  tkgrid(rejectButton, row = 0, column = 1, padx = 5)
  tkgrid(undoButton,   row = 0, column = 2, padx = 5)
}

 



# averaging Functions
# function to average selected groups for concatenated mode.
average_selected_groups <- function() {
  if (length(groups_list) == 0) {
    tkinsert(consoleText, 'end', "No groups available for averaging. Please select groups first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  
  dt_val      <- master_abf$samplingIntervalInSec * 1000
  stim_time   <- as.numeric(tclvalue(stimTimeVar))
  base_val    <- as.numeric(tclvalue(baselineVar))
  data_column <- as.numeric(tclvalue(dataColVar))
  if (is.na(data_column) || data_column < 1) data_column <- 1

  baseline2zero <- function(y, dt, stim, baseline) {
    idx_baseline <- round(baseline / dt)
    idx_start    <- round((stim - baseline) / dt) + 1
    y0 <- y - mean(y[1:idx_baseline])
    y0[idx_start:length(y0)]
  }

  averaged_data <<- lapply(groups_list, function(indices) {
    mats <- lapply(indices, function(i)
      baseline2zero(master_abf$data[[i]][, data_column],
                    dt_val, stim_time, base_val))
    rowMeans(do.call(cbind, mats))
  })

  current_avg_index <<- 1

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  tkgrid.columnconfigure(plotPanel, 0, weight = 1)
  tkgrid.rowconfigure(plotPanel, 0, weight = 1)

  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(avgFrame, 0, weight = 1)
  tkgrid.columnconfigure(avgFrame, 1, weight = 1)
  tkgrid.columnconfigure(avgFrame, 2, weight = 1)
  tkgrid.rowconfigure(avgFrame, 1, weight = 1)

  plotWrapper <- tkframe(avgFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  drawAvgPlot <- function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)
    
    y <- averaged_data[[current_avg_index]]
    dt_val <- master_abf$samplingIntervalInSec * 1000
    time <- seq(from = stim_time - base_val, by = dt_val, length.out = length(y))
    
    plot(time, y, type = 'l', col = 'darkgray', xlab = 'time (ms)', ylab = tclvalue(unitVar),
         xlim = smart_axis_limits(time), ylim = smart_axis_limits(y),
         axes = FALSE, bty = 'l')
    axis(1); axis(2, las = 1)
    stim_y <- y[which.min(abs(time - stim_time))]
    points(stim_time, stim_y, pch = 8, col = 'black')
    text(stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
  }

  avgPlot <<- tkrplot(plotWrapper, fun = drawAvgPlot, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(avgPlot, row = 0, column = 0, sticky = 'nsew')

  navFrame <- tkframe(avgFrame)
  tkgrid(navFrame, row = 2, column = 1, pady = 5)
  navLabel   <- tklabel(navFrame, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
  nextButton <- tkbutton(navFrame, text = 'Next', command = function() {
    current_avg_index <<- if (current_avg_index < length(averaged_data)) current_avg_index + 1 else 1
    tkconfigure(navLabel, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
    tkrreplot(avgPlot)
  })
  
  tkgrid(navLabel,   row = 0, column = 0, padx = 5)
  tkgrid(nextButton, row = 0, column = 1, padx = 5)

  tkdelete(consoleText, '1.0', 'end')
  tkinsert(consoleText, 'end', 'Averaging complete. Check the updated plot.')
  tkyview.moveto(consoleText, 1.0)
}

averageApprovedTraces_sep <- function() {
  if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
    tkinsert(consoleText, 'end', "No approved traces available. Please review recordings first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  result <- abf_averages(
    datasets         = abf_analysis_result$datasets,
    traces2average   = traces2average,
    baseline         = as.numeric(tclvalue(baselineVar)),
    stimulation_time = as.numeric(tclvalue(stimTimeVar)),
    dataCol          = as.numeric(tclvalue(dataColVar)),
    color            = 'darkgray',
    xbar             = as.numeric(tclvalue(xbarVar)),
    ybar             = as.numeric(tclvalue(ybarVar)),
    cex              = as.numeric(tclvalue(cexVar)),
    plotIt           = FALSE
  )
  averaged_data <<- result$baseline_corrected_mean_data
  datasets       <- result$datasets
  current_avg_index <<- 1

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  tkgrid.columnconfigure(plotPanel, 0, weight = 1)
  tkgrid.rowconfigure(   plotPanel, 0, weight = 1)

  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(avgFrame, 0, weight = 1)
  tkgrid.columnconfigure(avgFrame, 1, weight = 1)
  tkgrid.columnconfigure(avgFrame, 2, weight = 1)
  tkgrid.rowconfigure(   avgFrame, 1, weight = 1)

  plotWrapper <- tkframe(avgFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  drawSingleAvg <- function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)
    y    <- averaged_data[[current_avg_index]]
    dt_val <- datasets[[current_avg_index]]$samplingIntervalInSec * 1000
    time   <- seq(from = as.numeric(tclvalue(stimTimeVar)) - as.numeric(tclvalue(baselineVar)),
                  by   = dt_val,
                  length.out = length(y))

    egs_plot(x = time, y = y,
             color     = 'darkgray',
             show_bar  = TRUE,
             show_text = TRUE,
             xbar      = as.numeric(tclvalue(xbarVar)),
             ybar      = as.numeric(tclvalue(ybarVar)),
             xlim      = smart_axis_limits(time),
             ylim      = smart_axis_limits(y),
             cex       = cex)

    stim_y <- y[which.min(abs(time - as.numeric(tclvalue(stimTimeVar))))]
    points(as.numeric(tclvalue(stimTimeVar)), stim_y, pch = 8, col = 'black')
    text(  as.numeric(tclvalue(stimTimeVar)), stim_y, labels = 'stim', pos = 3, cex = cex)
  }

  avgPlot <<- tkrplot(plotWrapper, fun = drawSingleAvg,
                      hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(avgPlot, row = 0, column = 0, sticky = 'nsew')

  navFrame <- tkframe(avgFrame)
  tkgrid(navFrame, row = 2, column = 1, pady = 5)
  navLabel   <- tklabel(navFrame, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
  nextButton <- tkbutton(navFrame, text = 'Next', command = function() {
    current_avg_index <<- if (current_avg_index < length(averaged_data)) current_avg_index + 1 else 1
    tkconfigure(navLabel, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
    tkrreplot(avgPlot)
  })
  tkgrid(navLabel,   row = 0, column = 0, padx = 5)
  tkgrid(nextButton, row = 0, column = 1, padx = 5)

  tkdelete(consoleText, '1.0', 'end')
  tkinsert(consoleText, 'end', 'Separate-mode averaging complete. Check the updated plot.')
  tkyview.moveto(consoleText, 1.0)
}



  # UI Setup
ABF_analysis_tk <- function() {
    tt <- tktoplevel()
    tkwm.title(tt, 'ABF Analysis')
    
    if (.Platform$OS.type == "windows") {
      hscale <- 2
      vscale <- 2
    } else {
      dpi    <- as.numeric(tclvalue(tcl('winfo','pixels', tt, '1i')))
      w_in   <- 7;  h_in  <- 7
      hscale <- (w_in * dpi) / 480
      vscale <- (h_in * dpi) / 480
    }

    sidebarFrame <- tkframe(tt)
    mainFrame   <- tkframe(tt)
    tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
    tkgrid(mainFrame,   row = 0, column = 1, sticky = 'nsew')
    tkgrid.rowconfigure(tt, 0, weight = 1)
    tkgrid.columnconfigure(tt, 1, weight = 1)

    plotPanel <<- mainFrame

    ## --- folder selector ---
    folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
    tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
    folderPathVar <<- tclVar('')
    folderEntry <- ttkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
    tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
    browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function(){
      folderPath <- tclvalue(tkchooseDirectory())
      if (nchar(folderPath) > 0) {
        tclvalue(folderPathVar) <<- folderPath
        abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
        if (length(abf_list) == 0) {
          tkinsert(consoleText, 'end', 'No ABF files found in the selected folder.\n')
          tkyview.moveto(consoleText, 1.0)
        } else {
          tkdelete(abfListBox, 0, 'end')
          for (f in abf_list) tkinsert(abfListBox, 'end', f)
          firstFilePath <- file.path(folderPath, abf_list[1])
          ds <- readABF(firstFilePath)
          dummy_result <- list(metadata = list(extract_metadata(ds)))
          updateAdditionalParams(dummy_result)
        }
      }
    })
    tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)

    abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
    tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
    abfListBox <<- tklistbox(sidebarFrame, height = 5, selectmode = 'multiple')
    tkgrid(abfListBox, row = 2, column = 0, columnspan = 2, sticky = 'we', padx = 5, pady = 3)

    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

    paramFrame <- tkframe(sidebarFrame)

    ## --- parameters in paramFrame ---
    tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
    experimentCombo <- ttkcombobox(paramFrame,
      textvariable = experimentVar,
      values       = c('voltage clamp','current clamp'),
      width        = 15
    )
    tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
    unitEntry <- ttkentry(paramFrame,
      textvariable = unitVar,
      width        = 10
    )
    tcl(unitEntry, 'state', 'readonly')
    tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
    dataColEntry <- ttkentry(paramFrame, textvariable = dataColVar, width = 10)
    tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
    tkbind(dataColEntry, '<FocusOut>', function(...) {
      dc <- as.numeric(tclvalue(dataColVar))
      if (!exists('abf_analysis_result', envir = .GlobalEnv)) return()
      cu <- abf_analysis_result$datasets[[1]]$channelUnits
      if (!is.na(dc) && dc>=1 && dc<=length(cu)) tclvalue(unitVar) <<- cu[dc]
    })
    tkbind(dataColEntry,'<Return>',function(...) try(tcl("focus",""),silent=TRUE))

    tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
    dtEntry <- ttkentry(paramFrame, textvariable = dtVar, width = 10)
    tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
    tkbind(dtEntry,'<Return>',function(...){})

    tkgrid(tklabel(paramFrame, text = '# traces:'), row = 4, column = 0, sticky = 'w')
    ntracesEntry <- ttkentry(paramFrame, textvariable = ntracesVar, width = 10)
    tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
    tkbind(ntracesEntry,'<Return>',function(...){})

    tkgrid(paramFrame, row = 3, column = 0, columnspan = 2, sticky = 'we', pady = 3)
    tkgrid.columnconfigure(paramFrame, 0, weight = 1)
    tkgrid.columnconfigure(paramFrame, 1, weight = 1)

    ## --- rest of the sidebar ---
    baselineVar <<- tclVar('100')
    stimTimeVar <<- tclVar('150')
    xbarVar     <<- tclVar('100')
    ybarVar     <<- tclVar('50')
    concatMode  <<- tclVar('0')

    tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
    baselineEntry <- ttkentry(sidebarFrame, textvariable = baselineVar, width = 10)
    tkgrid(baselineEntry, row = 4, column = 1, sticky = 'w')
    tkbind(baselineEntry,'<Return>',function(...){})

    tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
    stimTimeEntry <- ttkentry(sidebarFrame, textvariable = stimTimeVar, width = 10)
    tkgrid(stimTimeEntry, row = 5, column = 1, sticky = 'w')
    tkbind(stimTimeEntry,'<Return>',function(...){})

    tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
    xbarEntry <- ttkentry(sidebarFrame, textvariable = xbarVar, width = 10)
    tkgrid(xbarEntry, row = 6, column = 1, sticky = 'w')
    tkbind(xbarEntry, '<FocusOut>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })
    tkbind(xbarEntry, '<Return>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })

    tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
    ybarEntry <- ttkentry(sidebarFrame, textvariable = ybarVar, width = 10)
    tkgrid(ybarEntry, row = 7, column = 1, sticky = 'w')
    tkbind(ybarEntry, '<FocusOut>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })
    tkbind(ybarEntry, '<Return>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })

    tkgrid(tklabel(sidebarFrame, text = 'Text scale (cex):'), row = 8, column = 0, sticky = 'w')
    cexEntry <- ttkentry(sidebarFrame, textvariable = cexVar, width = 10)
    tkgrid(cexEntry, row = 8, column = 1, sticky = 'w')
    tkbind(cexEntry, '<FocusOut>', function(...) {
      try({
        if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
        if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
      }, silent = TRUE)
    })

    tkbind(cexEntry, '<Return>', function(...) {
      try({
        if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
        if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
      }, silent = TRUE)
    })

    concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode,
                                  text = 'Concatenate Imported ABFs')
    tkgrid(concatButton, row = 9, column = 0, columnspan = 3)

    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

    consoleText <<- tktext(sidebarFrame, height = 5)
    tkgrid(consoleText, row = 10, column = 0, columnspan = 3,
           sticky = 'we', padx = 10, pady = 5)

    updateAdditionalParams <<- function(result) {
      if (!is.null(result) && length(result$metadata) >= 1) {
        meta1 <- result$metadata[[1]]
        tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
        if (as.character(tclvalue(concatMode)) != '1') {
          if (!is.null(meta1$header$lActualEpisodes))
            tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
          else
            tclvalue(ntracesVar) <<- 'N/A'
        }
        expType <- tclvalue(experimentVar)
        col_idx <- choose_data_column(meta1$channelUnits, expType)
        if (!is.na(col_idx)) {
          tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
          tclvalue(dataColVar) <<- as.character(col_idx)
        } else {
          tclvalue(unitVar) <<- 'N/A'
          tclvalue(dataColVar) <<- 'N/A'
        }
      }
    }
    tkbind(experimentCombo, '<<ComboboxSelected>>', function(widget, ...) {
      if (exists('abf_analysis_result', envir = .GlobalEnv)) {
        updateAdditionalParams(get('abf_analysis_result', envir = .GlobalEnv))
      }
    })

    runAnalysis <<- function() {
        # clear previous right‑panel widgets (graphs, metadata, table)
        children <- as.character(tkwinfo('children', plotPanel))
        if (length(children) > 0) {
          sapply(children, function(ch) tcl("destroy", ch))
        }

        folderPath <- tclvalue(folderPathVar)
        if (nchar(folderPath) == 0) {
          tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
          tkyview.moveto(consoleText, 1.0)
          return()
        }

      folderPath <- tclvalue(folderPathVar)
      if (nchar(folderPath) == 0) {
        tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      selIndices <- as.integer(tkcurselection(abfListBox))
      allFiles    <- as.character(tkget(abfListBox, 0, 'end'))
      abf_files   <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
      if (length(abf_files) == 0) {
        tkinsert(consoleText, 'end', 'No ABF files selected.\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      result <- tryCatch({
        load_abf_data(abf_files = abf_files, abf_path = folderPath)
      }, error = function(e) {
        tkinsert(consoleText, 'end', paste0('Error during data loading: ', e$message, '\n'))
        tkyview.moveto(consoleText, 1.0)
        NULL
      })
      if (!is.null(result)) {
        tkdelete(consoleText, '1.0', 'end')
        tkinsert(consoleText, 'end', paste0('Data loaded. Processed ', length(abf_files), ' file(s).\n'))
        tkyview.moveto(consoleText, 1.0)
        assign('abf_analysis_result', result, envir = .GlobalEnv)
        updateAdditionalParams(result)

        # display metadata
        meta1 <- result$metadata[[1]]
        first <- result$datasets[[1]]$data[[1]]
        length_sweep <- nrow(first)
        metaText <- paste(
          paste0("Format version: ", meta1$formatVersion),
          paste0("Sampling interval: ", meta1$samplingIntervalInSec, " s"),
          paste0("Channel names: ", paste(meta1$channelNames, collapse = " ")),
          paste0("Channel units: ", paste(meta1$channelUnits, collapse = " ")),
          paste0("Number of sweeps: ", meta1$header$lActualEpisodes),
          paste0("Length of first sweep: ", length_sweep),
          paste0("Path: ", meta1$path),
          sep = "\n"
        )
        kids <- as.character(tkwinfo('children', plotPanel))
        for (k in kids) tryCatch(tkdestroy(.Tk.ID[[k]]), error = function(e) {}, silent = TRUE)
        metaFrame  <- tkframe(plotPanel)
        tkgrid(metaFrame, row = 0, column = 0, sticky = 'w', pady = 2)
        metaLabel  <- tklabel(metaFrame, text = metaText, justify = 'left')
        tkgrid(metaLabel)

        # display first 10 rows of first trace
        out <- first[1:10, ]
        colnames(out) <- meta1$channelUnits
        rownames(out) <- seq(nrow(out))
        tableFrame  <- tkframe(plotPanel)
        tkgrid(tableFrame, row = 1, column = 0, sticky = 'nsew')
        textWidget  <<- tktext(tableFrame, width = 50, height = 11, wrap = 'none')
        tkgrid(textWidget, row = 0, column = 0)
        for (line in capture.output(print(out))) {
          tkinsert(textWidget, 'end', paste0(line, '\n'))
        }

        cons_msg <- check_consistency(result$metadata)
        if (cons_msg == 'Data is consistent') {
          tkinsert(consoleText, 'end', paste0(cons_msg, '\n'))
          tkyview.moveto(consoleText, 1.0)
        } else {
          tkinsert(consoleText, 'end', paste0('ERROR: ', cons_msg, '\n'))
          tkyview.moveto(consoleText, 1.0)
        }

        if (as.character(tclvalue(concatMode)) == '1') {
          master_abf <<- combine_abf_data(result)
          tclvalue(ntracesVar) <<- as.character(length(master_abf$data))
        } else {
          master_abf <<- result
        }
        tkconfigure(runAnalysisButton, text = 'Load Data')
      }
    }

    runAnalysisButton        <<- tkbutton(sidebarFrame, text = 'Load Data',               command = runAnalysis)
    reviewButton             <<- tkbutton(sidebarFrame, text = 'Review Recordings',        command = function() {
                                  if (as.character(tclvalue(concatMode)) == '1') review_master_recordings()
                                  else review_recordings()
                                })
    avgApprovedTracesButton  <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', command = function() {
                                  if (as.character(tclvalue(concatMode)) == '1') average_selected_groups()
                                  else averageApprovedTraces_sep()
                                })
    tkDownloadBtn            <<- tkbutton(sidebarFrame, text = 'Download Data',            command = download_data)

    tkgrid(runAnalysisButton,       row = 11, column = 0, columnspan = 3, pady = 5)
    tkgrid(reviewButton,            row = 12, column = 0, columnspan = 3, pady = 5)
    tkgrid(avgApprovedTracesButton, row = 13, column = 0, columnspan = 3, pady = 5)
    tkgrid(tkDownloadBtn,           row = 14, column = 0, columnspan = 3, pady = 5)

    tkfocus(tt)
    tkwait.window(tt)
  }

  # launch UI
  ABF_analysis_tk()
}



# UI Setup


# ABF_analysis_tk <- function() {
#     tt <- tktoplevel()
#     tkwm.title(tt, 'ABF Analysis')
    
#     if (.Platform$OS.type == "windows") {
#       hscale <- 2
#       vscale <- 2
#     } else {
#       dpi    <- as.numeric(tclvalue(tcl('winfo','pixels', tt, '1i')))
#       w_in   <- 7;  h_in  <- 7
#       hscale <- (w_in * dpi) / 480
#       vscale <- (h_in * dpi) / 480
#     }

#     sidebarFrame <- tkframe(tt)
#     mainFrame   <- tkframe(tt)
#     tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
#     tkgrid(mainFrame,   row = 0, column = 1, sticky = 'nsew')
#     tkgrid.rowconfigure(tt, 0, weight = 1)
#     tkgrid.columnconfigure(tt, 1, weight = 1)

#     plotPanel <<- mainFrame

#     ## --- folder selector ---
#     folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
#     tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
#     folderPathVar <<- tclVar('')
#     folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#     tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
#     browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function(){
#       folderPath <- tclvalue(tkchooseDirectory())
#       if (nchar(folderPath) > 0) {
#         tclvalue(folderPathVar) <<- folderPath
#         abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
#         if (length(abf_list) == 0) {
#           tkinsert(consoleText, 'end', 'No ABF files found in the selected folder.\n')
#           tkyview.moveto(consoleText, 1.0)
#         } else {
#           tkdelete(abfListBox, 0, 'end')
#           for (f in abf_list) tkinsert(abfListBox, 'end', f)
#           firstFilePath <- file.path(folderPath, abf_list[1])
#           ds <- readABF(firstFilePath)
#           dummy_result <- list(metadata = list(extract_metadata(ds)))
#           updateAdditionalParams(dummy_result)
#         }
#       }
#     })
#     tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)

#     abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
#     tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
#     abfListBox <<- tklistbox(sidebarFrame, height = 5, selectmode = 'multiple')
#     tkgrid(abfListBox, row = 2, column = 0, columnspan = 2, sticky = 'we', padx = 5, pady = 3)

#     tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
#     tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
#     tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

#     paramFrame <- tkframe(sidebarFrame)

#     ## --- parameters in paramFrame ---
#     tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
#     experimentCombo <- ttkcombobox(paramFrame,
#       textvariable = experimentVar,
#       values       = c('voltage clamp','current clamp'),
#       width        = 15
#     )
#     tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')

#     tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
#     unitEntry <- tkentry(paramFrame,
#       textvariable = unitVar,
#       width        = 10,
#       state        = 'readonly'
#     )
#     tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')

#     tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
#     dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10)
#     tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
#     tkbind(dataColEntry, '<FocusOut>', function(...) {
#       dc <- as.numeric(tclvalue(dataColVar))
#       if (!exists('abf_analysis_result', envir = .GlobalEnv)) return()
#       cu <- abf_analysis_result$datasets[[1]]$channelUnits
#       if (!is.na(dc) && dc>=1 && dc<=length(cu)) tclvalue(unitVar) <<- cu[dc]
#     })
#     tkbind(dataColEntry,'<Return>',function(...) try(tcl("focus",""),silent=TRUE))

#     tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
#     dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10)
#     tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
#     tkbind(dtEntry,'<Return>',function(...){})

#     tkgrid(tklabel(paramFrame, text = '# traces:'), row = 4, column = 0, sticky = 'w')
#     ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10)
#     tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
#     tkbind(ntracesEntry,'<Return>',function(...){})

#     tkgrid(paramFrame, row = 3, column = 0, columnspan = 2, sticky = 'we', pady = 3)
#     tkgrid.columnconfigure(paramFrame, 0, weight = 1)
#     tkgrid.columnconfigure(paramFrame, 1, weight = 1)

#     ## --- rest of the sidebar ---
#     baselineVar <<- tclVar('100')
#     stimTimeVar <<- tclVar('150')
#     xbarVar     <<- tclVar('100')
#     ybarVar     <<- tclVar('50')
#     concatMode  <<- tclVar('0')

#     tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
#     baselineEntry <- tkentry(sidebarFrame, textvariable = baselineVar, width = 10)
#     tkgrid(baselineEntry, row = 4, column = 1, sticky = 'w')
#     tkbind(baselineEntry,'<Return>',function(...){})

#     tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
#     stimTimeEntry <- tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10)
#     tkgrid(stimTimeEntry, row = 5, column = 1, sticky = 'w')
#     tkbind(stimTimeEntry,'<Return>',function(...){})

#     tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
#     xbarEntry <- tkentry(sidebarFrame, textvariable = xbarVar, width = 10)
#     tkgrid(xbarEntry, row = 6, column = 1, sticky = 'w')
#     tkbind(xbarEntry, '<FocusOut>', function(...) {
#       try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
#     })
#     tkbind(xbarEntry, '<Return>', function(...) {
#       try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
#     })

#     tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
#     ybarEntry <- tkentry(sidebarFrame, textvariable = ybarVar, width = 10)
#     tkgrid(ybarEntry, row = 7, column = 1, sticky = 'w')
#     tkbind(ybarEntry, '<FocusOut>', function(...) {
#       try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
#     })
#     tkbind(ybarEntry, '<Return>', function(...) {
#       try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
#     })

#     tkgrid(tklabel(sidebarFrame, text = 'Text scale (cex):'), row = 8, column = 0, sticky = 'w')
#     cexEntry <- tkentry(sidebarFrame, textvariable = cexVar, width = 10)
#     tkgrid(cexEntry, row = 8, column = 1, sticky = 'w')
#     tkbind(cexEntry, '<FocusOut>', function(...) {
#       try({
#         if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
#         if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
#       }, silent = TRUE)
#     })

#     tkbind(cexEntry, '<Return>', function(...) {
#       try({
#         if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
#         if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
#       }, silent = TRUE)
#     })

#     concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode,
#                                   text = 'Concatenate Imported ABFs')
#     tkgrid(concatButton, row = 9, column = 0, columnspan = 3)

#     tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
#     tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
#     tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

#     consoleText <<- tktext(sidebarFrame, height = 5)
#     tkgrid(consoleText, row = 10, column = 0, columnspan = 3,
#            sticky = 'we', padx = 10, pady = 5)

#     updateAdditionalParams <<- function(result) {
#       if (!is.null(result) && length(result$metadata) >= 1) {
#         meta1 <- result$metadata[[1]]
#         tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#         if (as.character(tclvalue(concatMode)) != '1') {
#           if (!is.null(meta1$header$lActualEpisodes))
#             tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#           else
#             tclvalue(ntracesVar) <<- 'N/A'
#         }
#         expType <- tclvalue(experimentVar)
#         col_idx <- choose_data_column(meta1$channelUnits, expType)
#         if (!is.na(col_idx)) {
#           tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#           tclvalue(dataColVar) <<- as.character(col_idx)
#         } else {
#           tclvalue(unitVar) <<- 'N/A'
#           tclvalue(dataColVar) <<- 'N/A'
#         }
#       }
#     }
#     tkbind(experimentCombo, '<<ComboboxSelected>>', function(widget, ...) {
#       if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#         updateAdditionalParams(get('abf_analysis_result', envir = .GlobalEnv))
#       }
#     })

#     runAnalysis <<- function() {
#         # clear previous right‑panel widgets (graphs, metadata, table)
#         children <- as.character(tkwinfo('children', plotPanel))
#         if (length(children) > 0) {
#           sapply(children, function(ch) tcl("destroy", ch))
#         }

#         folderPath <- tclvalue(folderPathVar)
#         if (nchar(folderPath) == 0) {
#           tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
#           tkyview.moveto(consoleText, 1.0)
#           return()
#         }

#       folderPath <- tclvalue(folderPathVar)
#       if (nchar(folderPath) == 0) {
#         tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
#         tkyview.moveto(consoleText, 1.0)
#         return()
#       }
#       selIndices <- as.integer(tkcurselection(abfListBox))
#       allFiles    <- as.character(tkget(abfListBox, 0, 'end'))
#       abf_files   <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#       if (length(abf_files) == 0) {
#         tkinsert(consoleText, 'end', 'No ABF files selected.\n')
#         tkyview.moveto(consoleText, 1.0)
#         return()
#       }
#       result <- tryCatch({
#         load_abf_data(abf_files = abf_files, abf_path = folderPath)
#       }, error = function(e) {
#         tkinsert(consoleText, 'end', paste0('Error during data loading: ', e$message, '\n'))
#         tkyview.moveto(consoleText, 1.0)
#         NULL
#       })
#       if (!is.null(result)) {
#         tkdelete(consoleText, '1.0', 'end')
#         tkinsert(consoleText, 'end', paste0('Data loaded. Processed ', length(abf_files), ' file(s).\n'))
#         tkyview.moveto(consoleText, 1.0)
#         assign('abf_analysis_result', result, envir = .GlobalEnv)
#         updateAdditionalParams(result)

#         # display metadata
#         meta1 <- result$metadata[[1]]
#         first <- result$datasets[[1]]$data[[1]]
#         length_sweep <- nrow(first)
#         metaText <- paste(
#           paste0("Format version: ", meta1$formatVersion),
#           paste0("Sampling interval: ", meta1$samplingIntervalInSec, " s"),
#           paste0("Channel names: ", paste(meta1$channelNames, collapse = " ")),
#           paste0("Channel units: ", paste(meta1$channelUnits, collapse = " ")),
#           paste0("Number of sweeps: ", meta1$header$lActualEpisodes),
#           paste0("Length of first sweep: ", length_sweep),
#           paste0("Path: ", meta1$path),
#           sep = "\n"
#         )
#         kids <- as.character(tkwinfo('children', plotPanel))
#         for (k in kids) tryCatch(tkdestroy(.Tk.ID[[k]]), error = function(e) {}, silent = TRUE)
#         metaFrame  <- tkframe(plotPanel)
#         tkgrid(metaFrame, row = 0, column = 0, sticky = 'w', pady = 2)
#         metaLabel  <- tklabel(metaFrame, text = metaText, justify = 'left')
#         tkgrid(metaLabel)

#         # display first 10 rows of first trace
#         out <- first[1:10, ]
#         colnames(out) <- meta1$channelUnits
#         rownames(out) <- seq(nrow(out))
#         tableFrame  <- tkframe(plotPanel)
#         tkgrid(tableFrame, row = 1, column = 0, sticky = 'nsew')
#         textWidget  <<- tktext(tableFrame, width = 50, height = 11, wrap = 'none')
#         tkgrid(textWidget, row = 0, column = 0)
#         for (line in capture.output(print(out))) {
#           tkinsert(textWidget, 'end', paste0(line, '\n'))
#         }

#         cons_msg <- check_consistency(result$metadata)
#         if (cons_msg == 'Data is consistent') {
#           tkinsert(consoleText, 'end', paste0(cons_msg, '\n'))
#           tkyview.moveto(consoleText, 1.0)
#         } else {
#           tkinsert(consoleText, 'end', paste0('ERROR: ', cons_msg, '\n'))
#           tkyview.moveto(consoleText, 1.0)
#         }

#         if (as.character(tclvalue(concatMode)) == '1') {
#           master_abf <<- combine_abf_data(result)
#           tclvalue(ntracesVar) <<- as.character(length(master_abf$data))
#         } else {
#           master_abf <<- result
#         }
#         tkconfigure(runAnalysisButton, text = 'Load Data')
#       }
#     }

#     runAnalysisButton        <<- tkbutton(sidebarFrame, text = 'Load Data',               command = runAnalysis)
#     reviewButton             <<- tkbutton(sidebarFrame, text = 'Review Recordings',        command = function() {
#                                   if (as.character(tclvalue(concatMode)) == '1') review_master_recordings()
#                                   else review_recordings()
#                                 })
#     avgApprovedTracesButton  <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', command = function() {
#                                   if (as.character(tclvalue(concatMode)) == '1') average_selected_groups()
#                                   else averageApprovedTraces_sep()
#                                 })
#     tkDownloadBtn            <<- tkbutton(sidebarFrame, text = 'Download Data',            command = download_data)

#     tkgrid(runAnalysisButton,       row = 11, column = 0, columnspan = 3, pady = 5)
#     tkgrid(reviewButton,            row = 12, column = 0, columnspan = 3, pady = 5)
#     tkgrid(avgApprovedTracesButton, row = 13, column = 0, columnspan = 3, pady = 5)
#     tkgrid(tkDownloadBtn,           row = 14, column = 0, columnspan = 3, pady = 5)

#     tkfocus(tt)
#     tkwait.window(tt)
#   }

#   # launch UI
#   ABF_analysis_tk()
# }  



############################################################################################
peak.fun2 <- function(y, dt, stimulation_time, baseline, detection_window=200, smooth=5){
  
  idx1 <- (stimulation_time - baseline) / dt
  idx2 <- baseline / dt
  idx3 <-  detection_window / dt

  y1 <- y[idx1:length(y)]
  bl <- mean(y1[0:idx2])

  y1 <- y1 - bl
  peak <- moving_avg(y1[idx2:(idx2 + idx3)], n = smooth)

  y2 <- y1[idx2:length(y1)]
  x <- (1:length(y2))*dt
  area <- trap_fun(x, y2)

  return(list(response=y1, peak=peak, charge=area, baseline=bl))
}

combine_abf_headers <- function(headers) {
  fields   <- names(headers[[1]])
  combined <- list()

  concat_sweeps <- function(v1, v2) {
    # round and coerce to integer
    v1_int <- as.integer(round(v1))
    v2_int <- as.integer(round(v2))

    # infer the uniform sweep interval
    interval <- v1_int[2] - v1_int[1]

    # compute the offset (last start of first recording)
    offset <- max(v1_int)

    # shift *all* v2 starts by offset+interval multiplicatively so
    # first is offset+interval, second is offset+2*interval, etc.
    v2_shifted <- offset + seq_along(v2_int) * interval

    # combine lengths: length(v1_int) + length(v2_int)
    c(v1_int, as.integer(v2_shifted))
  }

  for (f in fields) {
    vals <- lapply(headers, `[[`, f)

    # Scalar numeric fields: sum if not all identical
    if (all(sapply(vals, is.numeric)) && all(sapply(vals, length) == 1)) {
      if (all(vapply(vals[-1], identical, logical(1), vals[[1]]))) {
        combined[[f]] <- vals[[1]]
      } else {
        combined[[f]] <- sum(unlist(vals))
      }

    # sweepStartInPts: special continuous concatenation
    } else if (f == "sweepStartInPts") {
      v1 <- vals[[1]]
      v2 <- vals[[2]]
      combined[[f]] <- concat_sweeps(v1, v2)

    # Other numeric vectors (e.g. telegraph gains): keep identical vector
    } else if (all(sapply(vals, is.numeric)) &&
               all(vapply(vals, length, integer(1)) > 1) &&
               all(vapply(vals[-1], identical, logical(1), vals[[1]]))) {
      combined[[f]] <- vals[[1]]

    # 4) Fallback: store the list of all values
    } else {
      combined[[f]] <- vals
    }
  }

  combined
}


readABFs <- function(abf_files){

  if (length(abf_files)==1){
    out <- readABF(abf_files[1])
  }else{
    out_list <- lapply(abf_files , readABF)
    names(out_list) <- abf_files
    
    # create new abf style output
    master_data <- unlist(lapply(out_list, `[[`, "data"),recursive = FALSE)

    # confirm out_list entries share the same basic metadata
    same_samplingIntervalInSec <- all(sapply(out_list, function(x) 
      identical(x$samplingIntervalInSec, out_list[[1]]$samplingIntervalInSec)
    ))
    same_channelNames <- all(sapply(out_list, function(x) 
      identical(x$channelNames, out_list[[1]]$channelNames)
    ))
    same_channelUnits <- all(sapply(out_list, function(x) 
      identical(x$channelUnits, out_list[[1]]$channelUnits)
    ))
    same_nTelegraphEnable <- all(sapply(out_list, function(x) 
      identical(x$header$nTelegraphEnable, out_list[[1]]$header$nTelegraphEnable)
    ))
    same_fTelegraphAdditGain <- all(sapply(out_list, function(x) 
      identical(x$header$fTelegraphAdditGain, out_list[[1]]$header$fTelegraphAdditGain)
    ))
    same_fInstrumentScaleFactor <- all(sapply(out_list, function(x) 
      identical(x$header$fInstrumentScaleFactor, out_list[[1]]$header$fInstrumentScaleFactor)
    ))

    # If all checks pass, build the master object
    if (all(
      same_channelNames,
      same_channelUnits,
      same_nTelegraphEnable,
      same_fTelegraphAdditGain,
      same_fInstrumentScaleFactor
    )) {
      out <- out_list[[1]]
      out$data <- master_data
      headers    <- lapply(out_list, `[[`, "header")
      out$header <- combine_abf_headers(headers)
    } else {
      stop("ABF metadata not uniform across all files")
    }
  }
  out
}

analyseABF2 <- function() {

  # Increase max upload size (e.g., 100MB)
  options(shiny.maxRequestSize = 100*1024^2)  # 100MB in bytes

  # 30*1024^2 = 30MB
  # 50*1024^2 = 50MB
  # 100*1024^2 = 100MB
  # 500*1024^2 = 500MB
  # 1000*1024^2 = 1GB
    

  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
        @media (prefers-color-scheme: dark) {
          body { background-color: #1e1e1e; color: #e0e0e0; }
          .well { background-color: #2d2d2d; border-color: #444; }
          .form-control { background-color: #2d2d2d; color: #c0c0c0 !important; border: 1px solid #555 !important; }
          input[type='number'], input[type='text'] { background-color: #2d2d2d !important; color: #c0c0c0 !important; }
          .selectize-input, .selectize-dropdown { background-color: #ffffff !important; color: #666666 !important; }
          .btn { background-color: #3d3d3d; color: #ffffff; border-color: #555; }
          .btn-primary, .action-button { background-color: #3c8dbc; color: #ffffff; }
          pre, code { background-color: #1a1a1a; color: #f0f0f0; }
        }
      "))
    ),
    
    titlePanel("ABF Analysis"),
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput('abfFiles', 'Upload ABF Files', 
                  multiple = TRUE, 
                  accept = '.abf'),
        
        tabsetPanel(
          id = "settingsTabs",
          
          tabPanel("Main Settings",
            numericInput('baseline', 'Baseline (ms):', 100, min = 0),
            numericInput('stimulation', 'Stimulation Time (ms):', value = 150),
            checkboxInput('autoDetectStim', 'Auto-detect from TTL', TRUE),
            
            hr(),
            
            textInput('levels', 'Levels (comma-separated):', 'control,drug'),
            textInput('drugApplication', 'Drug Application Times (comma-separated):', ''),
            
            hr(),
            
            numericInput('pscChannel', 'PSC Channel:', 1, min = 1, max = 10),
            numericInput('hpChannel', 'Holding Potential Channel:', 2, min = 1, max = 10),
            numericInput('ttlChannel', 'TTL Channel (optional):', 3, min = 1, max = 10),
            
            hr(),
            verbatimTextOutput('fileInfo')
          ),
          
          tabPanel("Trace Selection",
        helpText("Define traces to average for each level"),
        uiOutput('traceSelectionUI'),
        hr(),
        verbatimTextOutput('tracesInfo')
      ),
          
          tabPanel("Graph Settings",
            numericInput('width', 'Plot Width:', 6, min = 1, max = 20),
            numericInput('height', 'Plot Height:', 8, min = 1, max = 20),
            
            hr(),
            h5("Peak Amplitude Plot"),
            numericInput('xmajor_tick_amp', 'X Tick:', 5, min = 1),
            numericInput('ymajor_tick_amp', 'Y Tick:', 100, min = 1),
            
            hr(),
            h5("Holding Current Plot"),
            numericInput('ymajor_tick_hc', 'Y Tick:', 100, min = 1),
            
            hr(),
            h5("Holding Potential Plot"),
            numericInput('xmajor_tick_hp', 'X Tick:', 5, min = 1),
            numericInput('ymajor_tick_hp', 'Y Tick:', 10, min = 1),
            
            hr(),
            h5("Appearance"),
            numericInput('cex_points', 'Point Size:', 1.5, min = 0.1, max = 3, step = 0.1),
            numericInput('lwd_graph', 'Line Thickness:', 1.5, min = 0.5, max = 5, step = 0.1),
            numericInput('cex_labels', 'Label Size:', 1.2, min = 0.5, max = 3, step = 0.1),
            numericInput('cex_axis', 'Axis Text Size:', 1.6, min = 0.5, max = 3, step = 0.1),
            
            hr(),
            textInput('xlab', 'X Label:', 'time (minutes)'),
            textInput('ylab', 'Y Label:', '|PSC| (pA)'),
            textInput('color', 'Plot Color:', '#CD5C5C')
          )
        ),
        
      hr(),

      actionButton('loadData', 'Load ABF Data', class = 'btn-primary'),
      actionButton('runAnalysis', 'Run Analysis', class = 'btn-primary'),

      hr(),

      downloadButton('downloadExcel', 'Download Excel'),
      downloadButton('downloadRData', 'Download RData'),
      actionButton('clearAll', 'Clear All', class = 'btn-default')
      ),
      
      mainPanel(
      width = 9,
      tabsetPanel(
        id = "mainTabs",
        
        tabPanel("Summary",
          h4("Data Summary"),
          verbatimTextOutput('summaryText'),
          hr(),
          tableOutput('summaryTable')
        ),
        
        tabPanel("Time Series Plot",
          plotOutput('timeSeriesPlot', height = '800px')
        ),
        
        tabPanel("Single Examples",
          fluidRow(
            column(12,
              downloadButton('downloadSVG', 'Download SVG Plots', class = 'btn-default')
            )
          ),
          hr(),
          uiOutput('examplePlotsUI')
        ),
        
        tabPanel("Review Traces",
          fluidRow(
            column(12,
              h4("Review Individual Traces"),
              helpText("Review and accept/reject individual traces for each level")
            )
          ),
          hr(),
          fluidRow(
            column(4,
              selectInput('reviewLevelSelect', 'Select Level to Review:', choices = NULL),
              actionButton('startReview', 'Start Review', class = 'btn-primary'),
              hr(),
              verbatimTextOutput('reviewProgress')
            ),
            column(8,
              plotOutput('reviewPlot', height = '500px'),
              hr(),
              fluidRow(
                column(6, actionButton('acceptTrace', 'Accept', class = 'btn-success', style = 'width: 100%;')),
                column(6, actionButton('rejectTrace', 'Reject', class = 'btn-danger', style = 'width: 100%;'))
              )
            )
          )
        ),
        
        tabPanel("Data Export",
          h4("Available Data"),
          verbatimTextOutput('exportInfo'),
          hr(),
          h5("Raw Data Preview"),
          tableOutput('rawDataPreview')
        )
      )
    )
    )
  )

  
  server <- function(input, output, session) {
    
    # Reactive values
    state <- reactiveValues(
      abf_dataset = NULL,
      I_data = NULL,
      I_data2 = NULL,
      holding_potential = NULL,
      holding_current = NULL,
      stimulation_time = NULL,
      dt = NULL,
      summary = NULL,
      Apeak = NULL,
      charge = NULL,
      traces2average = list(),
      split_include = list(),
      single_examples = NULL,
    review_level = NULL,    
    review_index = 1,       
    review_active = FALSE     
    )

  # Clear all data when new files are browsed
  observeEvent(input$abfFiles, {
    state$abf_dataset <- NULL
    state$I_data <- NULL
    state$I_data2 <- NULL
    state$holding_potential <- NULL
    state$holding_current <- NULL
    state$stimulation_time <- NULL
    state$dt <- NULL
    state$summary <- NULL
    state$Apeak <- NULL
    state$charge <- NULL
    state$traces2average <- list()
    state$split_include <- list()
    state$single_examples <- NULL
  })

  # Clear all button
  observeEvent(input$clearAll, {
    state$abf_dataset <- NULL
    state$I_data <- NULL
    state$I_data2 <- NULL
    state$holding_potential <- NULL
    state$holding_current <- NULL
    state$stimulation_time <- NULL
    state$dt <- NULL
    state$summary <- NULL
    state$Apeak <- NULL
    state$charge <- NULL
    state$traces2average <- list()
    state$split_include <- list()
    state$single_examples <- NULL
    
    updateTextInput(session, 'levels', value = 'control,drug')
    updateTextInput(session, 'drugApplication', value = '')
    updateNumericInput(session, 'baseline', value = 100)
    updateNumericInput(session, 'stimulation', value = 150)
    updateCheckboxInput(session, 'autoDetectStim', value = TRUE)
    updateNumericInput(session, 'pscChannel', value = 1)
    updateNumericInput(session, 'hpChannel', value = 2)
    updateNumericInput(session, 'ttlChannel', value = 3)
    
    levels <- c('control', 'drug')
    for (i in 1:10) {
      input_id <- paste0('traces_level_', i)
      if (!is.null(input[[input_id]])) {
        updateTextInput(session, input_id, value = '')
      }
    }
    
    showNotification("All data cleared", type = "message", duration = 2)
  })

    
    # File info output
    output$fileInfo <- renderPrint({
      req(state$abf_dataset)
      
      cat("Channels available:", ncol(state$abf_dataset$data[[1]]), "\n")
      cat("Traces:", length(state$abf_dataset$data), "\n")
      cat("Channel names:", paste(state$abf_dataset$channelNames, collapse = ", "), "\n")
      cat("Channel units:", paste(state$abf_dataset$channelUnits, collapse = ", "), "\n")
    })
    
    # Load ABF data
    observeEvent(input$loadData, {
      req(input$abfFiles)
      
      withProgress(message = 'Loading ABF files...', value = 0, {
        
        paths <- input$abfFiles$datapath
        names(paths) <- input$abfFiles$name
        
        incProgress(0.2, detail = "Reading files...")
        
        tryCatch({
          # Read ABF file(s)
          if (length(paths) == 1) {
            state$abf_dataset <- readABF(paths[1])
          } else {
            old_wd <- getwd()
            temp_dir <- dirname(paths[1])
            setwd(temp_dir)
            
            for (i in seq_along(paths)) {
              file.copy(paths[i], input$abfFiles$name[i], overwrite = TRUE)
            }
            
            state$abf_dataset <- readABFs(input$abfFiles$name)
            setwd(old_wd)
          }
          
          incProgress(0.4, detail = "Extracting data...")
          
          # Extract metadata
          state$dt <- state$abf_dataset$samplingIntervalInSec * 1000
          
          # Validate channel numbers
          n_channels <- ncol(state$abf_dataset$data[[1]])
          psc_ch <- as.integer(input$pscChannel)
          hp_ch <- as.integer(input$hpChannel)
          ttl_ch <- as.integer(input$ttlChannel)
          
          if (is.na(psc_ch) || psc_ch < 1 || psc_ch > n_channels) {
            stop(paste("PSC channel", psc_ch, "is invalid. Available channels: 1 to", n_channels))
          }
          
          if (is.na(hp_ch) || hp_ch < 1 || hp_ch > n_channels) {
            stop(paste("HP channel", hp_ch, "is invalid. Available channels: 1 to", n_channels))
          }
          
          # Extract PSC data
          state$I_data <- sapply(1:length(state$abf_dataset$data), function(ii) {
            state$abf_dataset$data[[ii]][, psc_ch]
          })
          colnames(state$I_data) <- seq(ncol(state$I_data))
          rownames(state$I_data) <- seq(nrow(state$I_data))
          
          # Extract holding potential
          state$holding_potential <- sapply(1:length(state$abf_dataset$data), function(ii) {
            state$abf_dataset$data[[ii]][, hp_ch]
          })
          
          incProgress(0.6, detail = "Detecting stimulation...")
          
          # Auto-detect stimulation from TTL if enabled
          if (input$autoDetectStim && !is.na(ttl_ch) && ttl_ch > 0 && ttl_ch <= n_channels) {
            tryCatch({
              TTL_pulse <- sapply(1:length(state$abf_dataset$data), function(ii) {
                state$abf_dataset$data[[ii]][, ttl_ch]
              })
              
              threshold <- 0.5
              pulse_on <- which(TTL_pulse[, 1] > threshold)
              
              if (length(pulse_on) > 0) {
                idx2 <- pulse_on[1]
                state$stimulation_time <- idx2 * state$dt - state$dt
                updateNumericInput(session, 'stimulation', value = round(state$stimulation_time, 2))
                showNotification(paste("Auto-detected stimulation at", round(state$stimulation_time, 2), "ms"), 
                               type = "message")
              } else {
                state$stimulation_time <- input$stimulation
                showNotification("No TTL pulse detected, using manual stimulation time", type = "warning")
              }
            }, error = function(e) {
              state$stimulation_time <- input$stimulation
              showNotification(paste("TTL detection failed:", e$message, "- using manual value"), 
                             type = "warning")
            })
          } else {
            state$stimulation_time <- input$stimulation
          }
          
          # Validate stimulation time
          if (is.null(state$stimulation_time) || is.na(state$stimulation_time)) {
            stop("Stimulation time must be specified")
          }
          
          incProgress(0.8, detail = "Processing traces...")
          
          # Process all traces
          out <- lapply(1:ncol(state$I_data), function(ii) {
            peak.fun2(state$I_data[, ii], 
                     dt = state$dt, 
                     stimulation_time = state$stimulation_time, 
                     baseline = input$baseline, 
                     smooth = 5)
          })
          
          state$Apeak <- sapply(out, function(x) x$peak)
          state$charge <- sapply(out, function(x) x$charge)
          state$I_data2 <- sapply(out, function(x) x$response)
          rownames(state$I_data2) <- seq(nrow(state$I_data2))
          colnames(state$I_data2) <- seq(ncol(state$I_data2))
          state$holding_current <- sapply(out, function(x) x$baseline)
          
          # Process holding potential
          out1 <- lapply(1:ncol(state$holding_potential), function(ii) {
            peak.fun2(state$holding_potential[, ii], 
                     dt = state$dt, 
                     stimulation_time = state$stimulation_time, 
                     baseline = input$baseline, 
                     smooth = 5)
          })
          state$holding_potential <- sapply(out1, function(x) x$baseline)
          
          # Create summary
          state$summary <- data.frame(
            'time_minutes' = 1:length(state$Apeak),
            'holding_potential_mV' = state$holding_potential,
            'holding_current_pA' = state$holding_current,
            'peak_amplitude_pA' = -state$Apeak,
            'charge_transfer_pC' = -state$charge,
            check.names = FALSE
          )
          
          incProgress(1.0, detail = "Done!")
          
          showNotification("ABF data loaded successfully!", type = "message", duration = 5)
          
        }, error = function(e) {
          showNotification(paste("Error loading ABF:", e$message), type = "error", duration = 10)
        })
      })
    })
    
    # Dynamic trace selection UI
    output$traceSelectionUI <- renderUI({
      req(state$I_data)
      n_traces <- ncol(state$I_data)
      
      tagList(
        helpText(paste("Total traces available:", n_traces)),
        uiOutput('dynamicTraceInputs')
      )
    })
    
  output$dynamicTraceInputs <- renderUI({
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    tagList(
      h5("Plot Axis Limits (applies to all Single Examples)"),
      fluidRow(
        column(6, 
          numericInput('xlim_min_all', 'X min:', value = NULL, step = 10)
        ),
        column(6, 
          numericInput('xlim_max_all', 'X max:', value = NULL, step = 10)
        )
      ),
      fluidRow(
        column(6, 
          numericInput('ylim_min_all', 'Y min:', value = NULL, step = 10)
        ),
        column(6, 
          numericInput('ylim_max_all', 'Y max:', value = NULL, step = 10)
        )
      ),
      hr(),
      h5("Scale Bar Settings"),
      fluidRow(
        column(6,
          numericInput('xbar_length', 'X bar:', value = 50, min = 1, step = 10)
        ),
        column(6,
          numericInput('ybar_length', 'Y bar:', value = 50, min = 1, step = 10)
        )
      ),
      fluidRow(
        column(6,
          numericInput('bar_lwd', 'Bar thickness:', value = 2, min = 0.5, max = 5, step = 0.5)
        )
      ),
      hr(),
      h5("Trace Selection"),
      lapply(seq_along(levels), function(i) {
        textInput(paste0('traces_level_', i), 
                 paste0('Traces for ', levels[i], ':'),
                 placeholder = 'e.g., 1,2,3,4,5 or 1:5')
      })
    )
  })
    
    # Parse traces2average from inputs
    parseTraces2Average <- reactive({
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    traces_list <- list()
    for (i in seq_along(levels)) {
      input_id <- paste0('traces_level_', i)
      trace_str <- input[[input_id]]
      
      if (!is.null(trace_str) && nchar(trace_str) > 0) {
        tryCatch({
          if (grepl(':', trace_str)) {
            parts <- strsplit(trace_str, ':')[[1]]
            start_val <- as.integer(parts[1])
            end_val <- as.integer(parts[2])
            
            if (!is.na(start_val) && !is.na(end_val)) {
              traces_list[[i]] <- start_val:end_val
            } else {
              traces_list[[i]] <- integer(0)
            }
          } else {
            vals <- as.integer(strsplit(trace_str, ',')[[1]])
            traces_list[[i]] <- vals[!is.na(vals)]
          }
        }, error = function(e) {
          traces_list[[i]] <- integer(0)
        })
      } else {
        traces_list[[i]] <- integer(0)
      }
    }
    
    traces_list
  })

    
    # Trace selection output
    output$tracesInfo <- renderPrint({
      traces_list <- parseTraces2Average()
      levels <- strsplit(input$levels, ',')[[1]]
      levels <- trimws(levels)
      
      cat("Configured trace groups:\n\n")
      for (i in seq_along(traces_list)) {
        if (length(traces_list[[i]]) > 0) {
          cat(levels[i], ": ", paste(traces_list[[i]], collapse = ", "), "\n")
        }
      }
    })
    
    # Run analysis
    observeEvent(input$runAnalysis, {
      req(state$I_data2, state$summary)
      
      withProgress(message = 'Running analysis...', value = 0, {
        
        traces2average <- parseTraces2Average()
        
        # Validate traces
        if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
          showNotification("Please specify traces to average in the Trace Selection tab", 
                         type = "error", duration = 5)
          return()
        }
        
        levels <- strsplit(input$levels, ',')[[1]]
        levels <- trimws(levels)
        
        incProgress(0.3, detail = "Averaging traces...")
        
        # Auto-include all traces
        split_include <- lapply(traces2average, function(traces) {
          rep(1, length(traces))
        })
        
        # Calculate averages
        result <- lapply(seq_along(split_include), function(iii) {
          mask <- split_include[[iii]]
          cols <- seq_along(mask)
          state$I_data2[, traces2average[[iii]], drop = FALSE][, cols[mask == 1], drop = FALSE]
        })
        
        state$average_traces <- sapply(result, function(mat) rowMeans(mat))
        
        time <- seq(nrow(state$I_data2)) * state$dt - state$dt
        state$single_examples <- cbind(time, state$average_traces)
        colnames(state$single_examples) <- c('time', levels)
        
        state$traces2average <- traces2average
        state$split_include <- split_include
        
        incProgress(1.0, detail = "Done!")
        
        showNotification("Analysis complete!", type = "message", duration = 3)
      })
    })

  # Auto-populate axis limits after analysis (common scale for all)
  observeEvent(state$single_examples, {
    req(state$single_examples)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    time <- state$single_examples[, 'time']
    xlim_common <- c(min(time), max(time))
    
    # Calculate common y limits across all levels
    all_y <- c()
    for (level in levels) {
      if (level %in% colnames(state$single_examples)) {
        all_y <- c(all_y, state$single_examples[, level])
      }
    }
    
    max_abs_y <- max(abs(all_y), na.rm = TRUE)
    ylim_common <- c(-max_abs_y * 1.1, 5)  # Negative-going trace, just +5 for noise
    
    # Round Y to nearest 10
    ylim_common[1] <- floor(ylim_common[1] / 10) * 10
    
    # Round X max to nearest 100
    xlim_common[2] <- ceiling(xlim_common[2] / 100) * 100
    
    updateNumericInput(session, 'xlim_min_all', value = round(xlim_common[1], 1))
    updateNumericInput(session, 'xlim_max_all', value = xlim_common[2])
    updateNumericInput(session, 'ylim_min_all', value = ylim_common[1])
    updateNumericInput(session, 'ylim_max_all', value = ylim_common[2])
  })

    
  # Update review level choices when analysis completes
  observe({
    req(state$traces2average)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    updateSelectInput(session, 'reviewLevelSelect', choices = levels, selected = levels[1])
  })

  # Start review
  observeEvent(input$startReview, {
    req(state$traces2average, state$split_include)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    level_idx <- which(levels == input$reviewLevelSelect)
    
    if (length(level_idx) == 0 || length(state$traces2average[[level_idx]]) == 0) {
      showNotification("No traces selected for this level", type = "warning")
      return()
    }
    
    state$review_level <- level_idx
    state$review_index <- 1
    state$review_active <- TRUE
    
    showNotification(paste("Reviewing", input$reviewLevelSelect), type = "message")
  })

  # Accept trace
  observeEvent(input$acceptTrace, {
    req(state$review_active, state$review_level)
    
    level_idx <- state$review_level
    trace_idx <- state$review_index
    
    state$split_include[[level_idx]][trace_idx] <- 1
    
    if (trace_idx < length(state$traces2average[[level_idx]])) {
      state$review_index <- trace_idx + 1
    } else {
      state$review_active <- FALSE
      recalculateAverages()
      showNotification("Review complete! Averages updated.", type = "message", duration = 5)
    }
  })

  # Reject trace
  observeEvent(input$rejectTrace, {
    req(state$review_active, state$review_level)
    
    level_idx <- state$review_level
    trace_idx <- state$review_index
    
    state$split_include[[level_idx]][trace_idx] <- 0
    
    if (trace_idx < length(state$traces2average[[level_idx]])) {
      state$review_index <- trace_idx + 1
    } else {
      state$review_active <- FALSE
      recalculateAverages()
      showNotification("Review complete! Averages updated.", type = "message", duration = 5)
    }
  })

  # Recalculate averages after review
  recalculateAverages <- function() {
    req(state$I_data2, state$traces2average, state$split_include)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    result <- lapply(seq_along(state$split_include), function(iii) {
      mask <- state$split_include[[iii]]
      cols <- seq_along(mask)
      accepted <- cols[mask == 1]
      
      if (length(accepted) > 0) {
        state$I_data2[, state$traces2average[[iii]], drop = FALSE][, accepted, drop = FALSE]
      } else {
        NULL
      }
    })
    
    non_empty <- !sapply(result, is.null)
    
    if (!any(non_empty)) {
      state$average_traces <- NULL
      state$single_examples <- NULL
      return(invisible(NULL))
    }
    
    avg_mat <- sapply(result[non_empty], function(mat) rowMeans(mat))
    if (is.null(dim(avg_mat))) {
      avg_mat <- matrix(avg_mat, ncol = 1)
    }
    
    time <- seq(nrow(state$I_data2)) * state$dt - state$dt
    state$single_examples <- cbind(time, avg_mat)
    
    colnames(state$single_examples) <- c('time', levels[non_empty])
    state$average_traces <- avg_mat
  }

  # Review progress output
  output$reviewProgress <- renderPrint({
    if (!state$review_active) {
      cat("Click 'Start Review' to begin\n")
      return()
    }
    
    req(state$review_level, state$traces2average)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    level_name <- levels[state$review_level]
    
    total_traces <- length(state$traces2average[[state$review_level]])
    current_trace <- state$traces2average[[state$review_level]][state$review_index]
    
    accepted <- sum(state$split_include[[state$review_level]] == 1)
    rejected <- sum(state$split_include[[state$review_level]] == 0)
    
    cat("Reviewing:", level_name, "\n")
    cat("Progress:", state$review_index, "/", total_traces, "\n")
    cat("Current trace:", current_trace, "\n\n")
    cat("Accepted:", accepted, "\n")
    cat("Rejected:", rejected, "\n")
  })

  # Review plot
  output$reviewPlot <- renderPlot({
    req(state$review_active, state$review_level, state$I_data2, state$single_examples)
    
    level_idx <- state$review_level
    trace_idx <- state$review_index
    actual_trace_num <- state$traces2average[[level_idx]][trace_idx]
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    level_name <- levels[level_idx]
    
    time <- state$single_examples[, 'time']
    avg_trace <- state$single_examples[, level_name]
    individual_trace <- state$I_data2[, actual_trace_num]
    
    # Use the same limits as Single Examples plots
    xlim_min <- input$xlim_min_all
    xlim_max <- input$xlim_max_all
    ylim_min <- input$ylim_min_all
    ylim_max <- input$ylim_max_all
    
    # Calculate defaults if not set
    if (is.null(xlim_min) || is.na(xlim_min)) xlim_min <- min(time)
    if (is.null(xlim_max) || is.na(xlim_max)) xlim_max <- max(time)
    
    if (is.null(ylim_min) || is.na(ylim_min) || is.null(ylim_max) || is.na(ylim_max)) {
      ylim_range <- range(c(avg_trace, individual_trace), na.rm = TRUE)
      ylim_min <- ylim_range[1] * 1.1
      ylim_max <- 5
    }
    
    xlim <- c(xlim_min, xlim_max)
    ylim <- c(ylim_min, ylim_max)
    
    # Get scale bar settings
    xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
    ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
    bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
    
    par(mar = c(2, 2, 3, 2))
    
    # Draw individual trace FIRST
    plot(time, individual_trace, type = 'l', col = '#CD5C5C', lwd = 2,
         ylim = ylim, xlim = xlim,
         xlab = '', ylab = '',
         main = paste(level_name, "- Trace", actual_trace_num),
         bty = 'n', axes = FALSE,
         cex.main = 1.3)
    
    abline(h = 0, lty = 2, col = 'grey70')
    
    # Draw average trace ON TOP
    lines(time, avg_trace, col = 'darkgrey', lwd = 3)
    
    # Add scale bars (no labels)
    usr <- par('usr')
    x_range <- usr[1:2]
    y_range <- usr[3:4]
    
    ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
    x_start <- x_range[2] - xbar - (x_range[2] - x_range[1]) * 0.05
    y_start <- ybar_start
    x_end <- x_start + xbar
    y_end <- y_start + ybar
    
    segments(x_start, y_start, x_end, y_start, lwd = bar_lwd, col = 'black')
    segments(x_start, y_start, x_start, y_end, lwd = bar_lwd, col = 'black')
    
    legend('bottomleft', 
           legend = c('Current Average', paste('Trace', actual_trace_num)),
           col = c('darkgrey', '#CD5C5C'),
           lwd = c(3, 2),
           bty = 'n',
           cex = 1.1)
  })

    # Summary output
    output$summaryText <- renderPrint({
      req(state$abf_dataset)
      
      cat("ABF Dataset Information\n")
      cat("=======================\n\n")
      cat("Sampling interval:", state$dt, "ms\n")
      cat("Number of traces:", ncol(state$I_data), "\n")
      cat("Trace length:", nrow(state$I_data), "samples\n")
      cat("Duration:", round(nrow(state$I_data) * state$dt / 1000, 2), "seconds\n")
      
      if (!is.null(state$stimulation_time)) {
        cat("Stimulation time:", state$stimulation_time, "ms\n")
      }
      
      if (!is.null(state$summary)) {
        cat("\nPeak amplitude range:", 
            round(min(state$summary$peak_amplitude_pA, na.rm = TRUE), 2), "to",
            round(max(state$summary$peak_amplitude_pA, na.rm = TRUE), 2), "pA\n")
      }
    })
    
    output$summaryTable <- renderTable({
      req(state$summary)
      head(state$summary, 20)
    }, striped = TRUE, hover = TRUE, bordered = TRUE)
   

  output$timeSeriesPlot <- renderPlot({
    req(state$summary)
    
    layout(matrix(1:3, ncol = 1), heights = c(3, 2, 2))
    
    tmax <- max(state$summary$time_minutes)
    xlim <- c(0, 5 * ceiling(tmax / 5))
    
    lwd_axes <- input$lwd_graph
    cex_lab <- input$cex_labels
    cex_ax <- input$cex_axis
    pt_size <- input$cex_points
    
    ylim1 <- c(0, ceiling(max(state$summary$peak_amplitude_pA, na.rm = TRUE) / 100) * 100)
    if (ylim1[2] == 0) ylim1[2] <- 100
    
    par(mar = c(0.5, 6, 3, 1), mgp = c(3.5, 0.7, 0))
    
    plot(state$summary$time_minutes, state$summary$peak_amplitude_pA,
         type = 'n', xlim = xlim, ylim = ylim1,
         xlab = '', ylab = '', main = '',
         axes = FALSE, xaxs = 'r', yaxs = 'r', cex.main = cex_lab)
    
    axis(2, at = seq(0, ylim1[2], by = input$ymajor_tick_amp),
         labels = seq(0, ylim1[2], by = input$ymajor_tick_amp),
         las = 1, tcl = -0.3, cex.axis = cex_ax, lwd = lwd_axes)
    
    mtext(input$ylab, side = 2, line = 4, cex = cex_lab)    
    
    points(state$summary$time_minutes, state$summary$peak_amplitude_pA,
           pch = 16, col = input$color, cex = pt_size)
    
    if (length(state$traces2average) > 0) {
      for (iii in seq_along(state$split_include)) {
        mask <- state$split_include[[iii]]
        cols <- which(mask == 1)
        if (length(cols) > 0 && length(state$traces2average[[iii]]) > 0) {
          exY <- state$summary$peak_amplitude_pA[state$traces2average[[iii]]][cols]
          exX <- state$summary$time_minutes[state$traces2average[[iii]]][cols]
          points(exX, exY, col = 'darkgrey', pch = 16, cex = pt_size * 1.2)
        }
      }
    }
    
    if (!is.null(input$drugApplication) && nchar(input$drugApplication) > 0) {
      drug_times <- as.numeric(strsplit(input$drugApplication, ',')[[1]])
      for (i in seq_along(drug_times)) {
        usr <- par('usr')
        rect(xleft = drug_times[i], 
             ybottom = usr[4] - (usr[4] - usr[3]) * 0.05,
             xright = tmax,
             ytop = usr[4],
             col = if (i %% 2 == 1) rgb(0.5, 0.5, 0.5, 0.3) else rgb(0.7, 0.7, 0.7, 0.3),
             border = NA)
      }
    }
    
    ylim2 <- c(floor(min(state$summary$holding_current_pA, na.rm = TRUE) / 100) * 100, 0)
    
    par(mar = c(0.5, 6, 0.5, 1), mgp = c(3.5, 0.7, 0))
    
    plot(state$summary$time_minutes, state$summary$holding_current_pA,
         type = 'n', xlim = xlim, ylim = ylim2,
         xlab = '', ylab = '', main = '',
         axes = FALSE, xaxs = 'r', yaxs = 'r')
    
    axis(2, at = seq(ylim2[1], 0, by = input$ymajor_tick_hc),
         labels = seq(ylim2[1], 0, by = input$ymajor_tick_hc),
         las = 1, tcl = -0.3, cex.axis = cex_ax, lwd = lwd_axes)
    
    mtext('Holding Current (pA)', side = 2, line = 4, cex = cex_lab)    
    
    points(state$summary$time_minutes, state$summary$holding_current_pA,
           pch = 16, col = input$color, cex = pt_size)
    
    if (length(state$traces2average) > 0) {
      for (iii in seq_along(state$split_include)) {
        mask <- state$split_include[[iii]]
        cols <- which(mask == 1)
        if (length(cols) > 0 && length(state$traces2average[[iii]]) > 0) {
          exY <- state$summary$holding_current_pA[state$traces2average[[iii]]][cols]
          exX <- state$summary$time_minutes[state$traces2average[[iii]]][cols]
          points(exX, exY, col = 'darkgrey', pch = 16, cex = pt_size * 1.2)
        }
      }
    }
    
    ylim3 <- range(state$summary$holding_potential_mV, na.rm = TRUE)
    ylim3 <- c(floor(ylim3[1] / 10) * 10 - 10, ceiling(ylim3[2] / 10) * 10 + 10)
    
    par(mar = c(5, 6, 0.5, 1), mgp = c(3.5, 0.7, 0))
    
    plot(state$summary$time_minutes, state$summary$holding_potential_mV,
         type = 'n', xlim = xlim, ylim = ylim3,
         xlab = '', ylab = '', main = '',
         axes = FALSE, xaxs = 'r', yaxs = 'r')
    
    axis(2, at = seq(ylim3[1], ylim3[2], by = input$ymajor_tick_hp),
         labels = seq(ylim3[1], ylim3[2], by = input$ymajor_tick_hp),
         las = 1, tcl = -0.3, cex.axis = cex_ax, lwd = lwd_axes)
    
    mtext('Holding Potential (mV)', side = 2, line = 4, cex = cex_lab)
    
    axis(1, at = seq(xlim[1], xlim[2], by = input$xmajor_tick_hp),
         labels = seq(xlim[1], xlim[2], by = input$xmajor_tick_hp),
         tcl = -0.3, cex.axis = cex_ax, lwd = lwd_axes)
    
    mtext(input$xlab, side = 1, line = 2.5, cex = cex_lab)
    
    points(state$summary$time_minutes, state$summary$holding_potential_mV,
           pch = 16, col = input$color, cex = pt_size)
    
    if (length(state$traces2average) > 0) {
      for (iii in seq_along(state$split_include)) {
        mask <- state$split_include[[iii]]
        cols <- which(mask == 1)
        if (length(cols) > 0 && length(state$traces2average[[iii]]) > 0) {
          exY <- state$summary$holding_potential_mV[state$traces2average[[iii]]][cols]
          exX <- state$summary$time_minutes[state$traces2average[[iii]]][cols]
          points(exX, exY, col = 'darkgrey', pch = 16, cex = pt_size * 1.2)
        }
      }
    }
    
    layout(1)
  })


  output$examplePlotsUI <- renderUI({
    req(state$single_examples)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    plot_outputs <- lapply(seq_along(levels), function(i) {
      plotOutput(paste0('examplePlot_', i), height = '400px')
    })
    
    do.call(fluidRow, lapply(plot_outputs, function(p) column(6, p)))
  })

  observe({
    req(state$single_examples)
    
    levels <- strsplit(input$levels, ',')[[1]]
    levels <- trimws(levels)
    
    x <- state$single_examples[, 'time']
    
    # Get common limits (same for all plots)
    xlim_min <- input$xlim_min_all
    xlim_max <- input$xlim_max_all
    ylim_min <- input$ylim_min_all
    ylim_max <- input$ylim_max_all
    
    # Calculate defaults if not set
    if (is.null(xlim_min) || is.na(xlim_min)) xlim_min <- min(x)
    if (is.null(xlim_max) || is.na(xlim_max)) xlim_max <- max(x)
    
    if (is.null(ylim_min) || is.na(ylim_min) || is.null(ylim_max) || is.na(ylim_max)) {
      all_y <- c()
      for (level in levels) {
        if (level %in% colnames(state$single_examples)) {
          all_y <- c(all_y, state$single_examples[, level])
        }
      }
      max_abs_y <- max(abs(all_y), na.rm = TRUE)
      ylim_min <- -max_abs_y * 1.1
      ylim_max <- 5
    }
    
    xlim <- c(xlim_min, xlim_max)
    ylim <- c(ylim_min, ylim_max)
    
    # Get scale bar settings
    xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
    ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
    bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
    
    # Create a plot for each level (all with same scale)
    for (i in seq_along(levels)) {
      local({
        my_i <- i
        my_level <- levels[my_i]
        
        output[[paste0('examplePlot_', my_i)]] <- renderPlot({
          if (my_level %in% colnames(state$single_examples)) {
            y <- state$single_examples[, my_level]
            
            par(mar = c(2, 2, 3, 2))
            
            plot(x, y, type = 'l', col = 'darkgrey', lwd = 1.5,
                 xlim = xlim, ylim = ylim,
                 xlab = '', ylab = '',
                 main = my_level,
                 bty = 'n', axes = FALSE,
                 cex.main = 1.3)
            
            abline(h = 0, lty = 2, col = 'grey')
            
            # Add scale bars (no labels)
            usr <- par('usr')
            x_range <- usr[1:2]
            y_range <- usr[3:4]
            
            ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
            x_start <- x_range[2] - xbar - (x_range[2] - x_range[1]) * 0.05
            y_start <- ybar_start
            x_end <- x_start + xbar
            y_end <- y_start + ybar
            
            segments(x_start, y_start, x_end, y_start, lwd = bar_lwd, col = 'black')
            segments(x_start, y_start, x_start, y_end, lwd = bar_lwd, col = 'black')
          }
        })
      })
    }
  })
   
  # Download handlers
  output$downloadExcel <- downloadHandler(
    filename = function() {
      # Generate filename based on loaded ABF files
      if (!is.null(input$abfFiles)) {
        file_names <- tools::file_path_sans_ext(input$abfFiles$name)
        
        if (length(file_names) == 1) {
          # Single file
          paste0(file_names[1], '_analysis.xlsx')
        } else {
          # Multiple files - use first and last
          first_name <- file_names[1]
          last_name <- file_names[length(file_names)]
          paste0(first_name, '_', last_name, '_analysis.xlsx')
        }
      } else {
        # Fallback if no files
        paste0('ABF_analysis_', Sys.Date(), '.xlsx')
      }
    },
    content = function(file) {
      req(state$summary, state$I_data, state$I_data2)
      
      wb <- createWorkbook()
      
      addWorksheet(wb, "Summary")
      writeData(wb, "Summary", state$summary)
      
      addWorksheet(wb, "Raw PSC Data")
      writeData(wb, "Raw PSC Data", as.data.frame(state$I_data))
      
      addWorksheet(wb, "Baseline Corrected")
      writeData(wb, "Baseline Corrected", as.data.frame(state$I_data2))
      
      if (!is.null(state$single_examples)) {
        addWorksheet(wb, "Single Examples")
        writeData(wb, "Single Examples", as.data.frame(state$single_examples))
      }
      
      # Add trace selection worksheet
      if (length(state$traces2average) > 0) {
        levels <- strsplit(input$levels, ',')[[1]]
        levels <- trimws(levels)
        
        trace_selection <- data.frame(
          Level = character(),
          Accepted = character(),
          Rejected = character(),
          stringsAsFactors = FALSE
        )
        
        for (i in seq_along(levels)) {
          if (length(state$traces2average[[i]]) > 0) {
            accepted_str <- ""
            rejected_str <- ""
            
            # Add accepted/rejected status if available
            if (length(state$split_include) >= i) {
              accepted <- state$traces2average[[i]][state$split_include[[i]] == 1]
              rejected <- state$traces2average[[i]][state$split_include[[i]] == 0]
              
              if (length(accepted) > 0) {
                accepted_str <- paste(accepted, collapse = ", ")
              }
              if (length(rejected) > 0) {
                rejected_str <- paste(rejected, collapse = ", ")
              }
            } else {
              # If no review done, all are accepted
              accepted_str <- paste(state$traces2average[[i]], collapse = ", ")
            }
            
            trace_selection <- rbind(trace_selection, 
              data.frame(Level = levels[i], 
                        Accepted = accepted_str,
                        Rejected = rejected_str))
          }
        }
        
        addWorksheet(wb, "Trace Selection")
        writeData(wb, "Trace Selection", trace_selection)
      }
      
      metadata <- data.frame(
        Parameter = c('dt (ms)', 'Stimulation Time (ms)', 'Baseline (ms)', 'N Traces'),
        Value = c(state$dt, state$stimulation_time, input$baseline, ncol(state$I_data))
      )
      addWorksheet(wb, "Metadata")
      writeData(wb, "Metadata", metadata)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
    
    output$downloadRData <- downloadHandler(
      filename = function() {
        paste0('ABF_analysis_', Sys.Date(), '.RData')
      },
      content = function(file) {
        results <- list(
          summary = state$summary,
          raw_data = state$I_data,
          baseline_corrected = state$I_data2,
          single_examples = state$single_examples,
          traces2average = state$traces2average,
          split_include = state$split_include,
          metadata = list(
            dt = state$dt,
            stimulation_time = state$stimulation_time,
            baseline = input$baseline
          )
        )
        save(results, file = file)
      }
    )

    output$downloadSVG <- downloadHandler(
      filename = function() {
        # Generate filename based on loaded ABF files
        if (!is.null(input$abfFiles)) {
          file_names <- tools::file_path_sans_ext(input$abfFiles$name)
          
          if (length(file_names) == 1) {
            paste0(file_names[1], '_plots.zip')
          } else {
            first_name <- file_names[1]
            last_name <- file_names[length(file_names)]
            paste0(first_name, '_', last_name, '_plots.zip')
          }
        } else {
          paste0('ABF_plots_', Sys.Date(), '.zip')
        }
      },
      content = function(file) {
        req(state$single_examples)
        
        levels <- strsplit(input$levels, ',')[[1]]
        levels <- trimws(levels)
        
        x <- state$single_examples[, 'time']
        
        # Get common limits
        xlim_min <- input$xlim_min_all
        xlim_max <- input$xlim_max_all
        ylim_min <- input$ylim_min_all
        ylim_max <- input$ylim_max_all
        
        if (is.null(xlim_min) || is.na(xlim_min)) xlim_min <- min(x)
        if (is.null(xlim_max) || is.na(xlim_max)) xlim_max <- max(x)
        
        if (is.null(ylim_min) || is.na(ylim_min) || is.null(ylim_max) || is.na(ylim_max)) {
          all_y <- c()
          for (level in levels) {
            if (level %in% colnames(state$single_examples)) {
              all_y <- c(all_y, state$single_examples[, level])
            }
          }
          max_abs_y <- max(abs(all_y), na.rm = TRUE)
          ylim_min <- -max_abs_y * 1.1
          ylim_max <- 5
        }
        
        xlim <- c(xlim_min, xlim_max)
        ylim <- c(ylim_min, ylim_max)
        
        # Get scale bar settings
        xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
        ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
        bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
        
        # Generate base filename
        if (!is.null(input$abfFiles)) {
          file_names <- tools::file_path_sans_ext(input$abfFiles$name)
          
          if (length(file_names) == 1) {
            base_name <- file_names[1]
          } else {
            first_name <- file_names[1]
            last_name <- file_names[length(file_names)]
            base_name <- paste0(first_name, '_', last_name)
          }
        } else {
          base_name <- paste0('ABF_', Sys.Date())
        }
        
        # Create temp directory for SVG files
        temp_dir <- tempdir()
        svg_files <- c()
        
        # Create separate SVG for each level
        for (level in levels) {
          if (level %in% colnames(state$single_examples)) {
            y <- state$single_examples[, level]
            
            svg_filename <- file.path(temp_dir, paste0(base_name, '_', level, '.svg'))
            svg_files <- c(svg_files, svg_filename)
            
            svg(svg_filename, width = 7, height = 5)
            par(mar = c(2, 2, 3, 2))
            
            plot(x, y, type = 'l', col = 'darkgrey', lwd = 1.5,
                 xlim = xlim, ylim = ylim,
                 xlab = '', ylab = '',
                 main = level,
                 bty = 'n', axes = FALSE,
                 cex.main = 1.3)
            
            abline(h = 0, lty = 2, col = 'grey')
            
            # Add scale bars
            usr <- par('usr')
            x_range <- usr[1:2]
            y_range <- usr[3:4]
            
            ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
            x_start <- x_range[2] - xbar - (x_range[2] - x_range[1]) * 0.05
            y_start <- ybar_start
            x_end <- x_start + xbar
            y_end <- y_start + ybar
            
            segments(x_start, y_start, x_end, y_start, lwd = bar_lwd, col = 'black')
            segments(x_start, y_start, x_start, y_end, lwd = bar_lwd, col = 'black')
            
            dev.off()
          }
        }
        
        # Create zip file
        zip(file, svg_files, flags = '-j')
      }
    )

    output$exportInfo <- renderPrint({
      if (is.null(state$summary)) {
        cat("No data loaded yet.\n")
        return()
      }
      
      cat("Available data for export:\n\n")
      cat("- Summary table (", nrow(state$summary), " rows)\n", sep = "")
      cat("- Raw PSC data (", ncol(state$I_data), " traces)\n", sep = "")
      cat("- Baseline corrected data\n")
      if (!is.null(state$single_examples)) {
        cat("- Single example traces (", ncol(state$single_examples) - 1, " levels)\n", sep = "")
      }
      cat("\nUse download buttons to export data.")
    })
    
    output$rawDataPreview <- renderTable({
      req(state$summary)
      head(state$summary, 10)
    }, striped = TRUE, hover = TRUE)
    
  }
  
  shinyApp(ui, server)
}


analyseABFshiny2 <- function() {

  # Increase max upload size (e.g., 100MB)
  options(shiny.maxRequestSize = 100*1024^2)  # 100MB in bytes

  # 30*1024^2 = 30MB
  # 50*1024^2 = 50MB
  # 100*1024^2 = 100MB
  # 500*1024^2 = 500MB
  # 1000*1024^2 = 1GB

  extract_metadata <- function(abf_dataset) {
    list(
      path                  = abf_dataset$path,
      formatVersion         = abf_dataset$formatVersion,
      channelNames          = abf_dataset$channelNames,
      channelUnits          = abf_dataset$channelUnits,
      samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
      header                = abf_dataset$header
    )
  }

  choose_data_column <- function(channelUnits, experiment) {
    if (experiment == 'voltage clamp') {
      idx <- grep('A', channelUnits, ignore.case = TRUE)
    } else if (experiment == 'current clamp') {
      idx <- grep('V', channelUnits, ignore.case = TRUE)
    } else {
      idx <- integer(0)
    }
    if (length(idx) > 0) return(idx[1])
    NA_integer_
  }

  combine_abf_data <- function(result) {
    master_abf <- list(data = list(), samplingIntervalInSec = result$datasets[[1]]$samplingIntervalInSec)
    for (i in seq_along(result$datasets)) {
      ds <- result$datasets[[i]]
      master_abf$data <- c(master_abf$data, ds$data)
    }
    master_abf
  }

  egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1,
                       show_text = FALSE, xbar = 100, ybar = 50,
                       color = 'darkgray', show_bar = FALSE, cex = 0.6) {
    if (is.null(ylim))
      ylim <- if (sign==1) c(0,max(y)) else c(-max(-y),0)
    if (is.null(xlim))
      xlim <- c(min(x), max(x))
    idx1 <- which.min(abs(x - xlim[1]))
    idx2 <- which.min(abs(x - xlim[2]))
    plot(x[idx1:idx2], y[idx1:idx2], type='l', col=color, xlim=xlim, ylim=ylim,
         axes=FALSE, xlab='time (ms)', ylab='')
    if (show_bar) {
      segments(max(x)-xbar, min(ylim), max(x), min(ylim))
      segments(max(x), min(ylim), max(x), min(ylim)+ybar)
      if (show_text) {
        text(max(x)-xbar/2, min(ylim)-0.05*diff(ylim), paste(xbar,'ms'))
        text(max(x)+0.02*diff(xlim), min(ylim)+ybar/2, paste(ybar,'pA'), srt=90)
      }
    }
  }

  ui <- fluidPage(
    # Add dark mode CSS
    tags$head(
      tags$style(HTML("
        /* Light mode (default) */
        body {
          background-color: white;
          color: black;
        }
        
        /* Dark mode - auto-detect browser preference */
        @media (prefers-color-scheme: dark) {
          body {
            background-color: #1e1e1e;
            color: #e0e0e0;
          }
          
          .well {
            background-color: #2d2d2d;
            border-color: #444;
          }
          
          /* Regular input fields - keep dark */
          .form-control {
            background-color: #2d2d2d;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* Numeric and text inputs - keep dark */
          input[type='number'],
          input[type='text'] {
            background-color: #2d2d2d !important;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* DROPDOWN MENUS ONLY - WHITE boxes with dark gray text */
          .selectize-input, .selectize-dropdown {
            background-color: #ffffff !important;
            color: #666666 !important;
            border: none !important;
          }

          /* Dropdown options */
          .selectize-dropdown .option {
            background-color: #ffffff;
            color: #666666;
          }

          .selectize-dropdown .option:hover {
            background-color: #f0f0f0;
            color: #333333;
          }

          /* Selected option in dropdown */
          .selectize-dropdown .selected,
          .selectize-dropdown .active {
            background-color: #e0e0e0;
            color: #333333;
          }

          /* Selectize item (selected chips) - NO GREY BOX */
          .selectize-input .item {
            background-color: #ffffff !important;
            color: #666666 !important;
            border: none !important;
          }
          
          /* Labels - also brighter */
          .shiny-input-container label, h4, h3 {
            color: #f0f0f0;
            font-weight: 500;
          }
          
          /* Radio buttons and checkbox labels */
          .radio label, .checkbox label {
            color: #f0f0f0;
          }
          
          /* Tab navigation */
          .nav-tabs {
            border-bottom-color: #444;
          }
          
          .nav-tabs > li > a {
            background-color: #2d2d2d;
            color: #e0e0e0;
            border-color: #444;
          }
          
          .nav-tabs > li.active > a,
          .nav-tabs > li.active > a:hover,
          .nav-tabs > li.active > a:focus {
            background-color: #1e1e1e;
            color: #fff;
            border-color: #444 #444 transparent;
          }
          
          /* Buttons */
          .btn {
            background-color: #3d3d3d;
            color: #ffffff;
            border-color: #555;
          }
          
          .btn:hover {
            background-color: #4d4d4d;
            border-color: #666;
            color: #ffffff;
          }
          
          .btn-default:hover {
            color: #fff;
          }
          
          /* Download buttons */
          .btn-default {
            color: #ffffff;
          }
          
          /* Action buttons - make them stand out more */
          .btn-primary, .action-button {
            background-color: #3c8dbc;
            color: #ffffff;
            border-color: #357ca5;
          }
          
          .btn-primary:hover, .action-button:hover {
            background-color: #4a9dd1;
            border-color: #428bca;
          }
          
          /* Verbatim output */
          pre, code {
            background-color: #1a1a1a;
            color: #f0f0f0;
            border-color: #444;
          }
          
          hr {
            border-top-color: #444;
          }
          
          /* File input styling */
          .btn-file {
            background-color: #3d3d3d;
            color: #ffffff;
          }
          
          /* Progress bar in dark mode */
          .shiny-progress .progress {
            background-color: #2d2d2d;
          }
          
          .shiny-progress .progress-bar {
            background-color: #3c8dbc;
          }
          
          .shiny-progress-notification {
            background-color: #2d2d2d;
            color: #f0f0f0;
            border-color: #444;
          }
          
          /* Focus states for inputs */
          .form-control:focus, .selectize-input.focus {
            border-color: #3c8dbc;
            box-shadow: 0 0 0 0.2rem rgba(60, 141, 188, 0.25);
          }
          
          /* Tables in dark mode */
          table {
            color: #e0e0e0;
          }
          
          .table-bordered {
            border-color: #444;
          }
          
          .table-bordered th,
          .table-bordered td {
            border-color: #444;
          }
        }
      "))
    ),
    
    titlePanel("ABF Analysis"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          id = "sideTabs",
          tabPanel("Main",
            fileInput("abfFiles","Upload ABF Files", multiple=TRUE, accept=".abf"),
            checkboxInput("concatenate","Concatenate ABFs",FALSE),
            actionButton("load","Load Data"), br(), br(),
            actionButton("review","Review Recordings"), br(), br(),
            actionButton("accept","Accept"), actionButton("reject","Reject"),
            actionButton("nextReview","Next"), br(), br(),
            actionButton("addGroup","Add Selected Group"), br(), br(),
            actionButton("completeSel","Selection Complete"), br(), br(),
            actionButton("average","Average Approved Traces"), br(), br(),
            actionButton("nextAvg","Next Average"), br(), br(),
            downloadButton("downloadData","Download Averaged CSV")
          ),
          tabPanel("Settings",
            selectInput("experiment","Experiment", c("voltage clamp","current clamp")),
            uiOutput("columnUI"),
            verbatimTextOutput("unitsText"),
            numericInput("dt","dt (ms)", NA),
            numericInput("ntraces","# traces", NA),
            numericInput("baseline","Baseline (ms)", 100),
            numericInput("stimTime","Stimulation time (ms)", 150),
            numericInput("xbar","x-bar length (ms)", 100),
            numericInput("ybar","y-bar length (pA)", 50)
          )
        )
      ),
      mainPanel(
        tabsetPanel(
          id = "mainTabs",
          tabPanel("Metadata",
            verbatimTextOutput("metaText"),
            tableOutput("firstTable")
          ),
          tabPanel("Review",
            fluidRow(
              column(8,   # ~66% for the plot
                plotOutput("reviewPlot", height="600px", width="100%")
              ),
              column(4,   # ~33% for the console box
                wellPanel(
                  style="height:600px; overflow:auto; padding:10px;",
                  verbatimTextOutput("console")
                )
              )
            )
          ),
          tabPanel("Average",
            fluidRow(
              column(8,
                plotOutput("avgPlot", height="600px", width="100%")
              ),
              column(4,
                wellPanel(
                  style="height:600px; overflow:auto; padding:10px;",
                  verbatimTextOutput("avgInfo")
                )
              )
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    vals <- reactiveValues(
      datasets=NULL, metadata=NULL,
      master=NULL, traces2avg=NULL,
      mode=NULL, ct=1, total=0,
      curGroup=NULL, groups=list(),
      avg=NULL, ca=1,
      log="", avgLog=""
    )

    baseline2zero <- function(y, dt, stim, baseline) {
      idx_baseline <- round(baseline / dt)
      idx_start    <- round((stim - baseline) / dt) + 1
      y0 <- y - mean(y[1:idx_baseline])
      y0[idx_start:length(y0)]
    }

    observeEvent(input$load, {
      req(input$abfFiles)
      paths <- input$abfFiles$datapath
      names(paths) <- input$abfFiles$name
      vals$datasets <- lapply(paths, readABF)
      vals$metadata <- lapply(vals$datasets, extract_metadata)
      m1 <- vals$metadata[[1]]
      updateNumericInput(session, "dt", value = m1$samplingIntervalInSec * 1000)
      updateNumericInput(session, "ntraces", value = m1$header$lActualEpisodes)
      sel <- choose_data_column(m1$channelUnits, input$experiment)
      updateSelectInput(session, "column", choices = seq_along(m1$channelUnits), selected = sel)
      vals$log <- ""
      updateTabsetPanel(session, "sideTabs", selected = "Settings")
    })

    observe({
      req(vals$datasets)    # only run once data are loaded

      total <- if (input$concatenate) {
        # sum traces across every ABF
        sum(vapply(vals$datasets, function(ds) length(ds$data), integer(1)))
      } else {
        # original sweep count in the first ABF
        vals$metadata[[1]]$header$lActualEpisodes
      }

      updateNumericInput(session, "ntraces", value = total)
    })

    output$columnUI <- renderUI({
      req(vals$metadata)
      selectInput("column","Data Column", seq_along(vals$metadata[[1]]$channelUnits))
    })
    output$unitsText <- renderText({
      req(input$column)
      paste0("Units: ", vals$metadata[[1]]$channelUnits[as.integer(input$column)])
    })
    output$metaText <- renderText({
      req(vals$metadata)
      m <- vals$metadata[[1]]; first <- vals$datasets[[1]]$data[[1]]
      c(
        paste0("Format version: ", m$formatVersion),
        paste0("Sampling interval: ", m$samplingIntervalInSec, " s"),
        paste0("Channel names: ", paste(m$channelNames, collapse=" ")),
        paste0("Channel units: ", paste(m$channelUnits, collapse=" ")),
        paste0("Number of sweeps: ", m$header$lActualEpisodes),
        paste0("Length of first sweep: ", nrow(first)),
        paste0("Path: ", m$path)
      )
    })
    output$firstTable <- renderTable({
      req(vals$datasets, input$column)
      out <- vals$datasets[[1]]$data[[1]][1:10, ]
      colnames(out) <- vals$metadata[[1]]$channelUnits
      head(out, 10)
    })

    # Review
    observeEvent(input$review, {
      req(vals$datasets, input$column)
      vals$mode <- if (input$concatenate) "concat" else "sep"
      if (vals$mode == "concat") {
        vals$master <- combine_abf_data(list(datasets = vals$datasets))
        vals$total  <- length(vals$master$data)
        vals$ct     <- 1
        vals$curGroup <- integer(0)
        vals$groups   <- list()
      } else {
        vals$curFile   <- 1; vals$ct <- 1
        vals$traces2avg <- vector("list", length(vals$datasets))
        for (i in seq_along(vals$traces2avg)) vals$traces2avg[[i]] <- integer(0)
      }
      vals$log <- ""
      updateTabsetPanel(session, "mainTabs", selected = "Review")
    })

    output$reviewPlot

      observeEvent(input$accept, {
        if (is.null(vals$mode)) {
          showNotification(
            "Error: No review in progress. Click 'Review Recordings' first.",
            type = "error",
            duration = 5
          )
          return()
        }

        if (vals$mode == "concat") {
          # concatenated mode
          vals$curGroup <- union(vals$curGroup, vals$ct)
          vals$log <- paste0(vals$log, "Accepted trace ", vals$ct, "\n")

          isolate({ input$nextReview })
          if (vals$ct < vals$total) {
            vals$ct <- vals$ct + 1
          } else {
            vals$log <- paste0(vals$log, "Review complete for all traces.\n")
          }

        } else {
          # separate mode
          fidx <- vals$curFile
          ds   <- vals$datasets[[fidx]]
          fname <- names(vals$datasets)[fidx]

          # log the accept
          vals$traces2avg[[fidx]] <- union(vals$traces2avg[[fidx]], vals$ct)
          vals$log <- paste0(vals$log, "Accepted ", fname, " trace ", vals$ct, "\n")

          isolate({ input$nextReview })
          if (vals$ct < length(ds$data)) {
            vals$ct <- vals$ct + 1
          } else if (vals$curFile < length(vals$datasets)) {
            # finished this file (but not the last) → log and advance
            vals$log <- paste0(vals$log, fname, " complete\n")
            vals$curFile <- vals$curFile + 1
            vals$ct      <- 1
          } else {
            # last file → log complete and final message
            vals$log <- paste0(vals$log, fname, " complete\n")
            vals$log <- paste0(vals$log, "Review complete: Approved recordings stored.\n")
          }
        }
      })

    observeEvent(input$reject, {
      if (is.null(vals$mode)) {
        showNotification(
          "Error: No review in progress. Click 'Review Recordings' first.",
          type = "error",
          duration = 5
        )
        return()
      }

      if (vals$mode == "concat") {
        # concatenated mode
        vals$log <- paste0(vals$log, "Rejected trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < vals$total) {
          vals$ct <- vals$ct + 1
        } else {
          vals$log <- paste0(vals$log, "Review complete for all traces.\n")
        }

      } else {
        # separate mode
        fidx <- vals$curFile
        ds   <- vals$datasets[[fidx]]
        fname <- names(vals$datasets)[fidx]

        # log the reject
        vals$log <- paste0(vals$log, "Rejected ", fname, " trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < length(ds$data)) {
          vals$ct <- vals$ct + 1
        } else if (vals$curFile < length(vals$datasets)) {
          # finished this file → log and advance
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$curFile <- vals$curFile + 1
          vals$ct      <- 1
        } else {
          # last file → log complete and final message
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$log <- paste0(vals$log, "Review complete: Approved recordings stored.\n")
        }
      }
    })

    output$console <- renderText(vals$log)

    observeEvent(input$addGroup, {
      if (length(vals$curGroup) == 0) {
        vals$log <- paste0(vals$log, "No traces selected in current group.\n")
      } else {
        vals$groups[[length(vals$groups) + 1]] <- vals$curGroup
        vals$log <- paste0(vals$log,
                           "Group ", length(vals$groups),
                           " selected: ",
                           paste(vals$curGroup, collapse = ","),
                           "\n")
        vals$curGroup <- integer(0)
      }
    })
    observeEvent(input$completeSel, {
      vals$log <- paste0(vals$log, "Review complete: Approved traces stored.\n")
    })

    # average traces
    observeEvent(input$average, {
      req(vals$mode)
      updateTabsetPanel(session, "mainTabs", selected = "Average")

      # common dt in ms
      dt <- if (vals$mode == "concat") {
        vals$master$samplingIntervalInSec * 1000
      } else {
        vals$datasets[[1]]$samplingIntervalInSec * 1000
      }

      if (vals$mode == "concat") {
        if (length(vals$groups) == 0) {
          vals$avgLog <- "No groups to average.\n"
          return()
        }
        vals$avg <- lapply(vals$groups, function(gr) {
          # combine accepted traces, compute mean
          y_full <- rowMeans(
            do.call(cbind, lapply(gr, function(i)
              vals$master$data[[i]][, as.integer(input$column)]))
          )
          baseline2zero(y_full, dt, input$stimTime, input$baseline)
        })

      } else {
        # separate files
        bc2 <- mapply(function(ds, idxs) {
          if (length(idxs) == 0) return(NULL)
          y_full <- rowMeans(
            do.call(cbind, lapply(idxs, function(i)
              ds$data[[i]][, as.integer(input$column)]))
          )
          baseline2zero(y_full, dt, input$stimTime, input$baseline)
        }, vals$datasets, vals$traces2avg,
        SIMPLIFY = FALSE)

        vals$avg <- bc2[!sapply(bc2, is.null)]
      }

      vals$ca     <- 1
      vals$avgLog <- "Averaging complete.\n"
    })

    observeEvent(input$nextAvg, {
      req(vals$avg)
      n <- length(vals$avg)
      vals$ca <- if (vals$ca < n) vals$ca + 1 else 1
    })

    output$avgPlot <- renderPlot({
      req(vals$avg)
      par(mar = c(2, 2, 1, 1))

      # data
      y   <- vals$avg[[vals$ca]]
      dt  <- vals$metadata[[1]]$samplingIntervalInSec * 1000
      time <- seq(0, by = dt, length.out = length(y))

      # plot trace
      egs_plot(time, y, show_bar  = FALSE, show_text = FALSE, color = 'darkgray')

      usr    <- par("usr")
      x_min  <- usr[1]; x_max <- usr[2]
      y_min  <- usr[3]; y_max <- usr[4]
      x_span <- x_max - x_min
      y_span <- y_max - y_min

      margin_x <- 0.05 * x_span
      margin_y <- 0.05 * y_span

      # bottom-right origin for bars
      x0 <- x_max - input$xbar - margin_x
      y0 <- y_min + margin_y

      # draw bars
      segments(x0, y0, x0 + input$xbar, y0, lwd = 1)           # xbar
      segments(x0, y0, x0, y0 + input$ybar, lwd = 1)           # ybar

      text(x0 + input$xbar/2,
           y0 - 0.03 * y_span,
           paste0(input$xbar, " ms"),
           adj = c(0.5, 1),
           cex = 1.2)

      text(
        x = x0 - margin_x/2,
        y = y0 + input$ybar/2,
        labels = paste0(input$ybar, " pA"),
        srt    = 90,
        adj    = c(0.5, 0.5),  # center in both directions
        cex    = 1.2)

      # stimulation marker at baseline
      text(input$baseline, 0, "*", col = "black", cex = 2.5)
      text(input$baseline, 0, labels = "stim", pos = 3, cex = 1)
    })

    output$avgInfo <- renderText({
      req(vals$avg)
      paste0("Average ", vals$ca, " of ", length(vals$avg), "\n", vals$avgLog)
    })

    # download as csv
    output$downloadData <- downloadHandler(
      filename = function() "averaged_data.csv",
      content = function(file) {
        # combine only the averaged traces (no time column)
        df <- as.data.frame(do.call(cbind, lapply(vals$avg, as.vector)))
        colnames(df) <- as.character(seq_along(vals$avg))
        write.csv(df, file, row.names = FALSE)
      }
    )
  }

  shinyApp(ui, server)

}

analysePSCtk <- function() {

  PSC_analysis_tk <- function() {
    tt <- tktoplevel()
    tkwm.title(tt, 'PSC Analysis')

    if (.Platform$OS.type == "windows") {
      hscale <- 2  
      vscale <- 2
    } else {
      # keep your old 3″×3″ DPI math on non-Windows if you like:
      dpi    <- as.numeric(tclvalue(tcl('winfo','pixels', tt, '1i')))
      w_in   <- 7; h_in <- 7
      hscale <- (w_in * dpi) / 480
      vscale <- (h_in * dpi) / 480
    }
    
    # divide window into sidebar and main panels
    sidebarFrame <- tkframe(tt)
    mainFrame <- tkframe(tt)
    plotWidget <<- NULL
    tkgrid(sidebarFrame, row=0, column=0, sticky='ns')
    tkgrid(mainFrame, row=0, column=1, sticky='nsew')
    tkgrid.rowconfigure(tt, 0, weight=0)
    tkgrid.columnconfigure(tt, 1, weight=1)
    
    # sidebar controls
    fileLabel <- tklabel(sidebarFrame, text='Upload csv or xlsx:')
    tkgrid(fileLabel, row=0, column=0, sticky='w')
    filePathVar <- tclVar('')
    fileEntry <- tkentry(sidebarFrame, textvariable=filePathVar, width=30)
    tkgrid(fileEntry, row=0, column=1, sticky='w')
    browseButton <- tkbutton(sidebarFrame, text='Browse', command=function() {
      filePath <- tclvalue(tkgetOpenFile(filetypes='{{CSV Files} {.csv}} {{Excel Files} {.xlsx .xls}}'))
      if (nchar(filePath) > 0) {
        tclvalue(filePathVar) <- filePath
        ext <- tools::file_ext(filePath)
        if (tolower(ext) == 'csv') {
          uploaded_data <<- read.csv(filePath)
        } else {
          uploaded_data <<- readxl::read_excel(filePath)
        }
        columns <<- names(uploaded_data)
        tkconfigure(columnCombo, values=columns)
      }
    })
    tkgrid(browseButton, row=0, column=2, padx=5)
    
    colLabel <- tklabel(sidebarFrame, text='Select column:')
    tkgrid(colLabel, row=1, column=0, sticky='w')
    columnVar <- tclVar('')
    columnCombo <- ttkcombobox(sidebarFrame, textvariable=columnVar, values='', width=20)
    tkgrid(columnCombo, row=1, column=1, columnspan=2, sticky='w')
    
    # Notebook for option tabs
    nb <- ttknotebook(sidebarFrame)
    tkgrid(nb, row=2, column=0, columnspan=3, pady=5, sticky='nsew')
    
    mainOptionsFrame   <- tkframe(nb)
    fitOptionsFrame    <- tkframe(nb)
    mleSettingsFrame   <- tkframe(nb)
    advancedFrame      <- tkframe(nb)
    graphSettingsFrame <- tkframe(nb)
    
    tkadd(nb, mainOptionsFrame, text='Main Options')
    tkadd(nb, fitOptionsFrame, text='Fit Options')
    tkadd(nb, mleSettingsFrame, text='MLE Settings')
    tkadd(nb, advancedFrame, text='Advanced')
    tkadd(nb, graphSettingsFrame, text='Plot Settings')
    
    # Main Options Tab
    dtVar <- tclVar('0.1')
    stimTimeVar <- tclVar('100')
    baselineVar <- tclVar('50')
    nVar <- tclVar('30')
    yAblineVar <- tclVar('0.1')
    funcVar <- tclVar('product1N')
    tkgrid(tklabel(mainOptionsFrame, text='dt (ms):'), row=0, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=dtVar, width=10), row=0, column=1)
    tkgrid(tklabel(mainOptionsFrame, text='Stimulation Time:'), row=1, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=stimTimeVar, width=10), row=1, column=1)
    tkgrid(tklabel(mainOptionsFrame, text='Baseline:'), row=2, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=baselineVar, width=10), row=2, column=1)
    tkgrid(tklabel(mainOptionsFrame, text='n:'), row=3, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=nVar, width=10), row=3, column=1)
    tkgrid(tklabel(mainOptionsFrame, text='Fit cutoff:'), row=4, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=yAblineVar, width=10), row=4, column=1)
    tkgrid(tklabel(mainOptionsFrame, text='Function:'), row=5, column=0, sticky='w')
    funcChoices <- c('product1N', 'product2N', 'product3N')
    funcCombo <- ttkcombobox(mainOptionsFrame, textvariable=funcVar, values=funcChoices, width=10)
    tkgrid(funcCombo, row=5, column=1)
    
    dsVar <- tclVar('1')
    tkgrid(tklabel(mainOptionsFrame, text='Downsample Factor:'), row=6, column=0, sticky='w')
    tkgrid(tkentry(mainOptionsFrame, textvariable=dsVar, width=10), row=6, column=1)
    
    # Fit Options Tab
    NVar <- tclVar('1')
    IEIVar <- tclVar('50')
    smoothVar <- tclVar('5')
    methodVar <- tclVar('BF.LM')
    weightMethodVar <- tclVar('none')
    sequentialFitVar <- tclVar('0')
    intervalMinVar <- tclVar('0.1')
    intervalMaxVar <- tclVar('0.9')
    lowerVar <- tclVar('')
    upperVar <- tclVar('')
    latencyLimitVar <- tclVar('')
    tkgrid(tklabel(fitOptionsFrame, text='N:'), row=0, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=NVar, width=10), row=0, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='IEI:'), row=1, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=IEIVar, width=10), row=1, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Smooth:'), row=2, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=smoothVar, width=10), row=2, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Method:'), row=3, column=0, sticky='w')
    methodChoices <- c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')
    methodCombo <- ttkcombobox(fitOptionsFrame, textvariable=methodVar, values=methodChoices, width=10)
    tkgrid(methodCombo, row=3, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Weighting:'), row=4, column=0, sticky='w')
    weightChoices <- c('none', '~y_sqrt', '~y')
    weightCombo <- ttkcombobox(fitOptionsFrame, textvariable=weightMethodVar, values=weightChoices, width=10)
    tkgrid(weightCombo, row=4, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Sequential Fit:'), row=5, column=0, sticky='w')
    sequentialFitCheck <- tkcheckbutton(fitOptionsFrame, variable=sequentialFitVar)
    tkgrid(sequentialFitCheck, row=5, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Min interval:'), row=6, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=intervalMinVar, width=10), row=6, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Max interval:'), row=7, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=intervalMaxVar, width=10), row=7, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Lower bounds:'), row=8, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=lowerVar, width=10), row=8, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Upper bounds:'), row=9, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=upperVar, width=10), row=9, column=1)
    tkgrid(tklabel(fitOptionsFrame, text='Latency limit:'), row=10, column=0, sticky='w')
    tkgrid(tkentry(fitOptionsFrame, textvariable=latencyLimitVar, width=10), row=10, column=1)
    
    # MLE Settings Tab
    iterVar <- tclVar('1000')
    metropolisScaleVar <- tclVar('1.5')
    fitAttemptsVar <- tclVar('10')
    RWmVar <- tclVar('0')
    tkgrid(tklabel(mleSettingsFrame, text='MLE Iterations:'), row=0, column=0, sticky='w')
    tkgrid(tkentry(mleSettingsFrame, textvariable=iterVar, width=10), row=0, column=1)
    tkgrid(tklabel(mleSettingsFrame, text='Metropolis Scale:'), row=1, column=0, sticky='w')
    tkgrid(tkentry(mleSettingsFrame, textvariable=metropolisScaleVar, width=10), row=1, column=1)
    tkgrid(tklabel(mleSettingsFrame, text='Fit Attempts:'), row=2, column=0, sticky='w')
    tkgrid(tkentry(mleSettingsFrame, textvariable=fitAttemptsVar, width=10), row=2, column=1)
    tkgrid(tklabel(mleSettingsFrame, text='Random Walk Metropolis:'), row=3, column=0, sticky='w')
    RWmCheck <- tkcheckbutton(mleSettingsFrame, variable=RWmVar)
    tkgrid(RWmCheck, row=3, column=1)
    
    # Advanced Tab
    filterVar <- tclVar('0')
    fcVar <- tclVar('1000')
    # relDecayFitLimitVar <- tclVar('0.1')
    halfWidthFitLimitVar <- tclVar('500')
    seedVar <- tclVar('42')
    dpVar <- tclVar('3')
    fastConstraintVar <- tclVar('0')
    fastConstraintMethodVar <- tclVar('rise')
    fastDecayLimitVar <- tclVar('')
    firstDelayConstraintVar <- tclVar('0')
    tkgrid(tklabel(advancedFrame, text='Filter:'), row=0, column=0, sticky='w')
    filterCheck <- tkcheckbutton(advancedFrame, variable=filterVar)
    tkgrid(filterCheck, row=0, column=1)
    tkgrid(tklabel(advancedFrame, text='Filter cutoff (Hz):'), row=1, column=0, sticky='w')
    tkgrid(tkentry(advancedFrame, textvariable=fcVar, width=10), row=1, column=1)
    tkgrid(tklabel(advancedFrame, text='Half-width fit limit:'), row=3, column=0, sticky='w')
    tkgrid(tkentry(advancedFrame, textvariable=halfWidthFitLimitVar, width=10), row=3, column=1)
    tkgrid(tklabel(advancedFrame, text='Seed:'), row=4, column=0, sticky='w')
    tkgrid(tkentry(advancedFrame, textvariable=seedVar, width=10), row=4, column=1)
    tkgrid(tklabel(advancedFrame, text='Decimal points:'), row=5, column=0, sticky='w')
    tkgrid(tkentry(advancedFrame, textvariable=dpVar, width=10), row=5, column=1)
    tkgrid(tklabel(advancedFrame, text='Fast constraint:'), row=6, column=0, sticky='w')
    fastConstraintCheck <- tkcheckbutton(advancedFrame, variable=fastConstraintVar)
    tkgrid(fastConstraintCheck, row=6, column=1)
    tkgrid(tklabel(advancedFrame, text='Fast constraint method:'), row=7, column=0, sticky='w')
    fastConstraintChoices <- c('rise', 'peak')
    fastConstraintCombo <- ttkcombobox(advancedFrame, textvariable=fastConstraintMethodVar, values=fastConstraintChoices, width=10)
    tkgrid(fastConstraintCombo, row=7, column=1)
    tkgrid(tklabel(advancedFrame, text='Fast decay limit(s):'), row=8, column=0, sticky='w')
    tkgrid(tkentry(advancedFrame, textvariable=fastDecayLimitVar, width=10), row=8, column=1)
    tkgrid(tklabel(advancedFrame, text='First delay constraint:'), row=9, column=0, sticky='w')
    firstDelayCheck <- tkcheckbutton(advancedFrame, variable=firstDelayConstraintVar)
    tkgrid(firstDelayCheck, row=9, column=1)
    
    # Graph Settings Tab
    lwdVar <- tclVar('1.2')
    xbarVar <- tclVar('50')
    ybarVar <- tclVar('50')
    xbarLabVar <- tclVar('ms')
    ybarLabVar <- tclVar('pA')

    tkgrid(tklabel(graphSettingsFrame, text='Line width:'), row=0, column=0, sticky='w')
    tkgrid(tkentry(graphSettingsFrame, textvariable=lwdVar, width=10), row=0, column=1)
    tkgrid(tklabel(graphSettingsFrame, text='x-bar length:'), row=1, column=0, sticky='w')
    tkgrid(tkentry(graphSettingsFrame, textvariable=xbarVar, width=10), row=1, column=1)
    tkgrid(tklabel(graphSettingsFrame, text='x-bar units:'), row=2, column=0, sticky='w')
    tkgrid(tkentry(graphSettingsFrame, textvariable=xbarLabVar, width=10), row=2, column=1)
    tkgrid(tklabel(graphSettingsFrame, text='y-bar length:'), row=3, column=0, sticky='w')
    tkgrid(tkentry(graphSettingsFrame, textvariable=ybarVar, width=10), row=3, column=1)
    tkgrid(tklabel(graphSettingsFrame, text='y-bar units:'), row=4, column=0, sticky='w')
    tkgrid(tkentry(graphSettingsFrame, textvariable=ybarLabVar, width=10), row=4, column=1)
    
    xlimVar <- tclVar('')
    tkgrid(tklabel(graphSettingsFrame, text='x limits (e.g., 25,400):'), row=7, column=0, sticky='w')
    xlimEntry <- tkentry(graphSettingsFrame, textvariable=xlimVar, width=20)
    tkgrid(xlimEntry, row=7, column=1)

    tkbind(xlimEntry, "<Return>", function() {
      if (!is.null(plotWidget)) tkrreplot(plotWidget, fun=drawPlotXlim, silent=TRUE)
    })

    cexVar <- if (Sys.info()["sysname"] == "Darwin") tclVar('0.6') else tclVar('1.4')
    tkgrid(
      tklabel(graphSettingsFrame, text='Text scale (cex):'),
      row=8, column=0, sticky='w'
    )
    cexEntry <- tkentry(graphSettingsFrame, textvariable=cexVar, width=20)
    tkgrid(cexEntry, row=8, column=1)
    tkbind(cexEntry, "<Return>", function(widget, ...) {
      if (!is.null(plotWidget)) {
        # replot initial view
        tkrreplot(plotWidget, fun=drawPlot1, silent=TRUE)
      }
    })


    # Additional sidebar controls
    userTmaxVar <- tclVar('')
    tkgrid(tklabel(sidebarFrame, text='User maximum time for fit:'), row=3, column=0, sticky='w', pady=5)
    tkgrid(tkentry(sidebarFrame, textvariable=userTmaxVar, width=10), row=3, column=1, pady=5)
    
    repeatConstraintVar <- tclVar('0')
    tkgrid(tklabel(sidebarFrame, text='Add fast constraint:'), row=4, column=0, sticky='w')
    repeatConstraintCheck <- tkcheckbutton(sidebarFrame, variable=repeatConstraintVar)
    tkgrid(repeatConstraintCheck, row=4, column=1)
    
    buttonFrame <- tkframe(sidebarFrame)
    tkgrid(buttonFrame, row=5, column=0, columnspan=3, pady=10, sticky='ew')

    tkgrid.columnconfigure(sidebarFrame, 0, weight=1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight=0)
    tkgrid.columnconfigure(sidebarFrame, 2, weight=1)

    tkgrid.columnconfigure(buttonFrame, 0, weight=1)
    tkgrid.columnconfigure(buttonFrame, 1, weight=0)
    tkgrid.columnconfigure(buttonFrame, 2, weight=0)
    tkgrid.columnconfigure(buttonFrame, 3, weight=1)

    # Analysis action buttons
    runAnalysisButton <- tkbutton(buttonFrame, text='Run Initial Analysis', command=function() {
      filePath <- tclvalue(filePathVar)
      if (nchar(filePath) == 0) {
        tkinsert(consoleText, 'end', 'Please select a file first\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      if (nchar(tclvalue(columnVar)) == 0) {
        tkinsert(consoleText, 'end', 'Please select a column\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      ext <- tools::file_ext(filePath)
      if (tolower(ext) == 'csv') {
        uploaded_data <<- read.csv(filePath)
      } else {
        uploaded_data <<- readxl::read_excel(filePath)
      }
      response_data <<- uploaded_data[[tclvalue(columnVar)]]
      ds <- as.numeric(tclvalue(dsVar))
      if (ds > 1) {
        response_data <<- response_data[seq(1, length(response_data), by=ds)]
      }
      # tkrreplot(plotWidget, fun=drawPlot1)

      if (is.null(plotWidget)) {
        plotWidget <<- tkrplot(
          mainFrame,
          fun    = drawPlot1,
          hscale = hscale,
          vscale = vscale
        )
        tkgrid(plotWidget, row=0, column=0, sticky='nsew')
      } else {
        tkrreplot(plotWidget, fun=drawPlot1)
      }

    })
      
    runMainAnalysisButton <- tkbutton(buttonFrame, text='Run Main Analysis', command=function() {
      fast.constraint        <- as.logical(as.numeric(tclvalue(repeatConstraintVar)))
      ds                     <- as.numeric(tclvalue(dsVar))
      dt                     <- as.numeric(tclvalue(dtVar)) * ds
      stimulation_time       <- as.numeric(tclvalue(stimTimeVar))
      baseline               <- as.numeric(tclvalue(baselineVar))
      smooth                 <- as.numeric(tclvalue(smoothVar))
      n                      <- as.numeric(tclvalue(nVar))
      N                      <- as.numeric(tclvalue(NVar))
      IEI                    <- as.numeric(tclvalue(IEIVar))
      func                   <- get(tclvalue(funcVar))
      method                 <- tclvalue(methodVar)
      weight_method          <- tclvalue(weightMethodVar)
      sequential.fit         <- as.logical(as.numeric(tclvalue(sequentialFitVar)))
      fit.limits             <- as.numeric(tclvalue(userTmaxVar))
      rel.decay.fit.limit    <- as.numeric(tclvalue(yAblineVar))
      lwd                    <- as.numeric(tclvalue(lwdVar))
      fc                     <- as.numeric(tclvalue(fcVar))
      interval               <- c(as.numeric(tclvalue(intervalMinVar)), as.numeric(tclvalue(intervalMaxVar)))
      lower                  <- if (nchar(tclvalue(lowerVar)) > 0) as.numeric(unlist(strsplit(tclvalue(lowerVar), ','))) else NULL
      upper                  <- if (nchar(tclvalue(upperVar)) > 0) as.numeric(unlist(strsplit(tclvalue(upperVar), ','))) else NULL
      iter                   <- as.numeric(tclvalue(iterVar))
      metropolis.scale       <- as.numeric(tclvalue(metropolisScaleVar))
      fit.attempts           <- as.numeric(tclvalue(fitAttemptsVar))
      RWm                    <- as.logical(as.numeric(tclvalue(RWmVar)))
      fast.decay.limit       <- if (nchar(tclvalue(fastDecayLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ','))) else NULL
      fast.constraint.method <- tclvalue(fastConstraintMethodVar)
      first.delay.constraint <- as.logical(as.numeric(tclvalue(firstDelayConstraintVar)))
      dp                     <- as.numeric(tclvalue(dpVar))
      seed                   <- as.numeric(tclvalue(seedVar))
      filter                 <- as.logical(as.numeric(tclvalue(filterVar)))
      
      y <- response_data
      if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
        y <- y[!is.na(y)]
      }
      x <- seq(0, (length(y) - 1) * dt, by=dt)
      
      if (!sequential.fit) {
        tmax <- fit.limits
        x_limit <- determine_tmax2(y=y, N=N, stimulation_time=stimulation_time, baseline=baseline, lwd=lwd, 
                                  smooth=smooth, tmax=tmax, y_abline=rel.decay.fit.limit, xbar=as.numeric(tclvalue(xbarVar)),
                                  ybar=as.numeric(tclvalue(ybarVar)), xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
        adjusted_response <- y[x < x_limit]
        
        out <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func, method=method,
                    weight_method=weight_method, MLEsettings=list(iter=iter, metropolis.scale=metropolis.scale, 
                    fit.attempts=fit.attempts, RWm=RWm), stimulation_time=stimulation_time, baseline=baseline,
                    filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, 
                    fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, 
                    first.delay.constraint=first.delay.constraint, lower=lower, upper=upper,
                    latency.limit=if (nchar(tclvalue(latencyLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ','))) else NULL,
                    return.output=TRUE, show.plot=FALSE, half_width_fit_limit=as.numeric(tclvalue(halfWidthFitLimitVar)),
                    dp=dp, height=5, width=5,seed=seed)
        
        out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time,
                                  baseline=baseline, func=func, filter=filter, fc=fc)
        xlim_input <- tclvalue(xlimVar)
        if (nchar(xlim_input) > 0) {
          xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ',')))
          if (length(xlim_vals) == 2) {
            out$traces <- out$traces[out$traces$x >= xlim_vals[1] & out$traces$x <= xlim_vals[2], ]
          }
        }

        tkrreplot(plotWidget, fun=function() {
          drawPlot2(traces=out$traces, func=func, lwd=lwd, cex=as.numeric(tclvalue(cexVar)), filter=filter,
                    xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                    xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
        })

      } else {
        out <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method,
                    stimulation_time=stimulation_time, baseline=baseline, fit.limits=fit.limits, fast.decay.limit=fast.decay.limit,
                    fast.constraint=as.logical(as.numeric(tclvalue(fastConstraintVar))), fast.constraint.method=fast.constraint.method,
                    first.delay.constraint=first.delay.constraint, latency.limit=if (nchar(tclvalue(latencyLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ','))) else NULL,
                    lower=lower, upper=upper, filter=filter, fc=fc, interval=interval,
                    MLEsettings=list(iter=iter, metropolis.scale=metropolis.scale, fit.attempts=fit.attempts, RWm=RWm),
                    MLE.method=method, half_width_fit_limit=as.numeric(tclvalue(halfWidthFitLimitVar)),
                    dp=dp, lwd=lwd, xlab='', ylab='', width=5, height=5, return.output=TRUE, show.output=TRUE,
                    show.plot=TRUE, seed=seed)
                }
      
      analysis_output <<- out
      df_out <- out$output
      if(sum(grepl('^A\\d+$', names(df_out))) == 1) names(df_out)[which(grepl('^A\\d+$', names(df_out)))] <- 'A'
      if(sum(grepl('^area\\d+$', names(df_out)))==1) names(df_out)[which(grepl('^area\\d+$', names(df_out)))] <- 'area'
      names(df_out) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df_out))
      names(df_out) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df_out))
      names(df_out)[names(df_out) == 'half_width'] <- 'half width'
      # tkdelete(consoleText, '1.0', 'end')
      # tkinsert(consoleText, 'end', 'Analysis complete.')
      tkdelete(fitOutputText, '1.0', 'end')
      tkinsert(fitOutputText, 'end', paste(capture.output(print(df_out)), collapse='\n'))
    })
    
    downloadOutputButton <- tkbutton(buttonFrame, text='Download RData', 
      command=function() {
        if (!exists('analysis_output') || is.null(analysis_output)) {
          tkinsert(consoleText, 'end', 'No analysis output available!\n')
          tkyview.moveto(consoleText, 1.0)
          return()
        }
        saveFile <- tclvalue(tkgetSaveFile(filetypes='{{Rdata Files} {.Rdata}} {{All Files} *}'))
        if (nchar(saveFile) > 0) {
          # Build metadata list using the tk variable values.
          metadata <- list(
            dt=as.numeric(tclvalue(dtVar)),
            stimTime=as.numeric(tclvalue(stimTimeVar)),
            baseline=as.numeric(tclvalue(baselineVar)),
            n=as.numeric(tclvalue(nVar)),
            y_abline=as.numeric(tclvalue(yAblineVar)),
            func=tclvalue(funcVar),
            ds=as.numeric(tclvalue(dsVar)),
            N=as.numeric(tclvalue(NVar)),
            IEI=as.numeric(tclvalue(IEIVar)),
            smooth=as.numeric(tclvalue(smoothVar)),
            method=tclvalue(methodVar),
            weight_method=tclvalue(weightMethodVar),
            sequential_fit=as.logical(as.numeric(tclvalue(sequentialFitVar))),
            interval=c(as.numeric(tclvalue(intervalMinVar)), as.numeric(tclvalue(intervalMaxVar))),
            lower=if (nchar(tclvalue(lowerVar)) > 0)
                      as.numeric(unlist(strsplit(tclvalue(lowerVar), ",")))
                    else NULL,
            upper=if (nchar(tclvalue(upperVar)) > 0)
                      as.numeric(unlist(strsplit(tclvalue(upperVar), ",")))
                    else NULL,
            latency_limit=if (nchar(tclvalue(latencyLimitVar)) > 0)
                              as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
                            else NULL,
            iter=as.numeric(tclvalue(iterVar)),
            metropolis_scale=as.numeric(tclvalue(metropolisScaleVar)),
            fit_attempts=as.numeric(tclvalue(fitAttemptsVar)),
            RWm=as.logical(as.numeric(tclvalue(RWmVar))),
            fast_decay_limit=if (nchar(tclvalue(fastDecayLimitVar)) > 0)
                                 as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ",")))
                               else NULL,
            fast_constraint=as.logical(as.numeric(tclvalue(fastConstraintVar))),
            fast_constraint_method=tclvalue(fastConstraintMethodVar),
            first_delay_constraint=as.logical(as.numeric(tclvalue(firstDelayConstraintVar))),
            dp=as.numeric(tclvalue(dpVar)),
            seed=as.numeric(tclvalue(seedVar)),
            filter=as.logical(as.numeric(tclvalue(filterVar))),
            fc=as.numeric(tclvalue(fcVar)),
            userTmax=as.numeric(tclvalue(userTmaxVar)),
            data_col=tclvalue(columnVar)
          )
          # Combine analysis output and metadata into one list.
          results <- list(
            analysis=analysis_output,
            metadata=metadata
          )
          save(results, file=saveFile)
          tkinsert(consoleText, 'end', 'Output saved successfully.\n')
          tkyview.moveto(consoleText, 1.0)
        }
      }
    )
    
      downloadResultsButton <- tkbutton(buttonFrame, text='Download Output (csv/xlsx)', command=function() {
        if (!exists('analysis_output') || is.null(analysis_output)) {
          tkinsert(consoleText, 'end', 'No analysis output available!\n')
          tkyview.moveto(consoleText, 1.0)
          return()
        }
        filePath <- tclvalue(tkgetSaveFile(filetypes='{{Excel File} {.xlsx}} {{CSV File} {.csv}}'))
        if (nchar(filePath) == 0) return()
        ext <- tolower(tools::file_ext(filePath))

        data_list <- list(
          output          = analysis_output$output,
          traces          = analysis_output$traces,
          `fit criterion` = data.frame(AIC = analysis_output$AIC, BIC = analysis_output$BIC),
          `model message` = data.frame(message = analysis_output$model.message)
        )

        metadata_labels <- c(
          'Data column:','dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:','Function:',
          'Downsample Factor:','User maximum time for fit:','Add fast constraint:','N:','IEI:','Smooth:',
          'Method:','Weighting:','Sequential Fit:','Min interval:','Max interval:',
          'Lower bounds (comma-separated):','Upper bounds (comma-separated):','Latency limit:',
          'MLE Iterations:','Metropolis Scale:','Fit Attempts:','Random Walk Metropolis:',
          'Filter:','Filter cutoff (Hz):','Half-width fit limit:','Seed:','Decimal points:',
          'Fast constraint method:','Fast decay limit(s):','First delay constraint:'
        )

        metadata_values <- c(
          tclvalue(columnVar),tclvalue(dtVar),tclvalue(stimTimeVar),tclvalue(baselineVar),
          tclvalue(nVar),tclvalue(yAblineVar),tclvalue(funcVar),tclvalue(dsVar),
          tclvalue(userTmaxVar),as.character(as.logical(as.numeric(tclvalue(repeatConstraintVar)))),
          tclvalue(NVar),tclvalue(IEIVar),tclvalue(smoothVar),tclvalue(methodVar),
          tclvalue(weightMethodVar),as.character(as.logical(as.numeric(tclvalue(sequentialFitVar)))),
          tclvalue(intervalMinVar),tclvalue(intervalMaxVar),tclvalue(lowerVar),tclvalue(upperVar),
          tclvalue(latencyLimitVar),tclvalue(iterVar),tclvalue(metropolisScaleVar),
          tclvalue(fitAttemptsVar),as.character(as.logical(as.numeric(tclvalue(RWmVar)))),
          as.character(as.logical(as.numeric(tclvalue(filterVar)))),tclvalue(fcVar),
          tclvalue(halfWidthFitLimitVar),tclvalue(seedVar),tclvalue(dpVar),
          tclvalue(fastConstraintMethodVar),tclvalue(fastDecayLimitVar),
          as.character(as.logical(as.numeric(tclvalue(firstDelayConstraintVar))))
        )

        coerce_val <- function(v) {
          if (nzchar(v)) {
            n <- suppressWarnings(as.numeric(v))
            if (!is.na(n)) return(n)
          }
          v
        }

        if (ext == 'csv') {
          for (nm in names(data_list)) {
            write.csv(data_list[[nm]], file = sub("\\.csv$", paste0("_", nm, ".csv"), filePath), row.names = FALSE)
          }
          meta_df <- data.frame(name = metadata_labels, value = sapply(metadata_values, coerce_val), stringsAsFactors = FALSE)
          write.csv(meta_df, file = sub("\\.csv$", "_metadata.csv", filePath), row.names = FALSE)

        } else if (ext == 'xlsx') {
          wb <- createWorkbook()
          for (nm in names(data_list)) {
            addWorksheet(wb, nm)
            writeData(wb, nm, data_list[[nm]])
          }
          addWorksheet(wb, "metadata")
          writeData(wb, "metadata", c("parameter", "value"), startRow = 1, startCol = 1, colNames = FALSE)
          for (i in seq_along(metadata_labels)) {
            writeData(wb, "metadata", metadata_labels[i], startRow = i + 1, startCol = 1, colNames = FALSE)
            writeData(wb, "metadata", coerce_val(metadata_values[i]), startRow = i + 1, startCol = 2, colNames = FALSE)
          }
          saveWorkbook(wb, filePath, overwrite = TRUE)

        } else {
          tkinsert(consoleText, 'end', 'Unsupported file type. Use .csv or .xlsx\n')
          tkyview.moveto(consoleText, 1.0)
        }
      })

    exportSVGButton <- tkbutton(buttonFrame, text='Export Plot to SVG', command=function() {
      if (!exists('analysis_output') || is.null(analysis_output)) {
        tkinsert(consoleText, 'end', 'No analysis available to export!\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      saveFile <- tclvalue(tkgetSaveFile(filetypes='{{SVG Files} {.svg}}'))
      if (nchar(saveFile) > 0) {
        # Parse and apply xlim if provided
        traces <- analysis_output$traces
        xlim_input <- tclvalue(xlimVar)
        if (nchar(xlim_input) > 0) {
          xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ',')))
          if (length(xlim_vals) == 2) {
            traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
          }
        }
        
        svg(filename=saveFile, width=7, height=5)
        drawPlot2(traces=traces, func=get(tclvalue(funcVar)), lwd=as.numeric(tclvalue(lwdVar)), cex=0.6,
                  filter=as.logical(as.numeric(tclvalue(filterVar))),
                  xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                  xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
        dev.off()
        tkinsert(consoleText, 'end', 'SVG plot saved successfully.\n')
        tkyview.moveto(consoleText, 1.0)
      }
    })

    
    clearOutputButton <- tkbutton(buttonFrame, text='Clear Output', command=function() {
      analysis_output <<- NULL
      # tkdelete(consoleText, '1.0', 'end')
      tkdelete(fitOutputText, '1.0', 'end')
      # tkrreplot(plotWidget, fun=drawPlot1)
      # plotWidget <<- tkrplot(mainFrame, fun=drawPlot1)
      # tkgrid(    plotWidget,    row=0, column=0, sticky='nsew')
      if (!is.null(plotWidget)) {
        tkrreplot(plotWidget, fun=drawPlot1)
      }

    })
    
    tkgrid(runAnalysisButton,      row=0, column=0, padx=5, pady=2)
    tkgrid(runMainAnalysisButton,  row=0, column=1, padx=5, pady=2)
    tkgrid(downloadResultsButton,  row=1, column=1, padx=5, pady=2)
    tkgrid(downloadOutputButton,   row=1, column=0, padx=5, pady=2)
    tkgrid(exportSVGButton,        row=2, column=0, padx=5, pady=2)
    tkgrid(clearOutputButton,      row=2, column=1, padx=5, pady=2)

    # Main Panel plot
    drawPlot1 <- function() {
      
      ds <- as.numeric(tclvalue(dsVar))
      dt <- as.numeric(tclvalue(dtVar)) * ds
      lwd <- as.numeric(tclvalue(lwdVar))
      stimTime <- as.numeric(tclvalue(stimTimeVar))
      baseline <- as.numeric(tclvalue(baselineVar))
      smooth <- as.numeric(tclvalue(smoothVar))
      y_abline <- as.numeric(tclvalue(yAblineVar))
      y_val <- if (exists('response_data') && !is.null(response_data)) response_data else rnorm(10000, 0.1)
      cex <- as.numeric(tclvalue(cexVar))
      
      determine_tmax2(y=y_val, N=1, dt=dt, stimulation_time=stimTime, baseline=baseline, smooth=smooth, lwd=lwd,
        tmax=NULL, y_abline=y_abline, xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
        xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar), cex=cex)
    }
    
    drawPlotXlim <- function() {
      xlim_input <- tclvalue(xlimVar)
      cex <- as.numeric(tclvalue(cexVar))
      traces <- if (exists("analysis_output") && !is.null(analysis_output)) analysis_output$traces else NULL
      if (is.null(traces)) return()

      if (nchar(xlim_input) > 0) {
        xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ",")))
        if (length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }
      }

      drawPlot2(traces=traces, func=get(tclvalue(funcVar)), lwd=as.numeric(tclvalue(lwdVar)), cex=cex,
                filter=as.logical(as.numeric(tclvalue(filterVar))),
                xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
    }
    # plotWidget <- tkrplot(tt, fun=drawPlot1)
    # tkgrid(plotWidget, row=0, column=1, sticky='nsew')
    
    # consoleText <- tktext(mainFrame, width=80, height=4)
    # tkgrid(consoleText, row=1, column=1, sticky='nsew')
    
    fitOutputLabel <- tklabel(sidebarFrame, text='Fit Output:')
    tkgrid(fitOutputLabel, row=11, column=0, columnspan=3, sticky='w', pady=c(10,2), padx=20)
    
    fitOutputText <- tktext(sidebarFrame, width=90, height=5)
    tkgrid(fitOutputText, row=12, column=0, columnspan=3, sticky='w', padx=20)
    
    tkfocus(tt)
    tkwait.window(tt)
  }

  PSC_analysis_tk()

}


# MODULAR VERSION OF analysePSC()



# UI HELPER FUNCTIONS


create_dark_mode_css_psc <- function() {
  tags$head(
    tags$style(HTML("
      @media (prefers-color-scheme: dark) {
        body { background-color: #1e1e1e; color: #e0e0e0; }
        .well { background-color: #2d2d2d; border-color: #444; }
        .form-control { background-color: #2d2d2d; color: #c0c0c0 !important; border: 1px solid #555 !important; }
        input[type='number'], input[type='text'] { background-color: #2d2d2d !important; color: #c0c0c0 !important; }
        .selectize-input, .selectize-dropdown { background-color: #ffffff !important; color: #666666 !important; }
        .btn { background-color: #3d3d3d; color: #ffffff; border-color: #555; }
        .btn-primary, .action-button { background-color: #3c8dbc; color: #ffffff; }
        pre, code { background-color: #1a1a1a; color: #f0f0f0; }
      }
    "))
  )
}

create_main_options_ui <- function() {
  tabPanel('Main Options',
    numericInput('dt', 'dt (ms):', 0.1),
    numericInput('stimulation_time', 'Stimulation Time:', 100),
    numericInput('baseline', 'Baseline:', 50),
    numericInput('n', 'n:', 30),
    numericInput('y_abline', 'Fit Cutoff:', 0.1),
    selectInput('func', 'Function:', choices=c('product1N', 'product2N', 'product3N')),
    numericInput('ds', 'Downsample Factor:', 1, min=1)
  )
}

create_fit_options_ui <- function() {
  tabPanel('Fit Options',
    numericInput('N', 'N:', 1),
    numericInput('IEI', 'IEI:', 50),
    numericInput('smooth', 'Smooth:', 5),
    selectInput('method', 'Method:', choices=c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')),
    selectInput('weight_method', 'Weighting:', choices=c('none', '~y_sqrt', '~y')),
    checkboxInput('sequential_fit', 'Sequential Fit', FALSE),
    numericInput('interval_min', 'Min Interval:', 0.1),
    numericInput('interval_max', 'Max Interval:', 0.9),
    textInput('lower', 'Lower Bounds (comma-separated):', ''),
    textInput('upper', 'Upper Bounds (comma-separated):', ''),
    textInput('latency_limit', 'Latency Limit:', '')
  )
}

create_mle_settings_ui <- function() {
  tabPanel('MLE Settings',
    numericInput('iter', 'MLE Iterations:', 1000),
    numericInput('metropolis_scale', 'Metropolis Scale:', 1.5),
    numericInput('fit_attempts', 'Fit Attempts:', 10),
    checkboxInput('RWm', 'Random Walk Metropolis', FALSE)
  )
}

create_advanced_ui <- function() {
  tabPanel('Advanced',
    checkboxInput('filter', 'Filter', FALSE),
    numericInput('fc', 'Filter Cutoff (Hz):', 1000),
    numericInput('half_width_fit_limit', 'Half-width Fit Limit:', 500),
    numericInput('seed', 'Seed:', 42),
    numericInput('dp', 'Decimal Points:', 3),
    checkboxInput('fast_constraint', 'Fast Constraint', FALSE),
    selectInput('fast_constraint_method', 'Fast Constraint Method:', choices=c('rise', 'peak')),
    textInput('fast_decay_limit', 'Fast Decay Limit(s) (comma-separated):', ''),
    checkboxInput('first_delay_constraint', 'First Delay Constraint', FALSE)
  )
}

create_plot_settings_ui <- function() {
  tabPanel('Plot Settings',
    numericInput('lwd', 'Line Width:', 1.2),
    numericInput('xbar', 'x-bar Length:', 50),
    numericInput('ybar', 'y-bar Length:', 50),
    textInput('xbar_lab', 'x-axis Units:', 'ms'),
    textInput('ybar_lab', 'y-axis Units:', 'pA'),
    textInput('xlim', 'x limits (e.g., 0,400):', '')
  )
}

create_sidebar_panel_psc <- function() {
  sidebarPanel(
    fileInput('file', 'Upload csv or xlsx', accept=c('.csv', '.xlsx')),
    uiOutput('column_selector'),
    tabsetPanel(
      create_main_options_ui(),
      create_fit_options_ui(),
      create_mle_settings_ui(),
      create_advanced_ui(),
      create_plot_settings_ui()
    ),
    numericInput('userTmax', 'User Maximum Time for Fit:', NA),
    actionButton('run_initial', 'Run Initial Analysis'),
    actionButton('run_main', 'Run Main Analysis'),
    actionButton('add_result', 'Add to Results'),
    actionButton('clear_results', 'Clear Results'), 
    downloadButton('download_xlsx',  'Download Output (*.xlsx)'),
    downloadButton('download_output', 'Download RData'),
    downloadButton('download_svg', 'Download SVG Plot'),
    actionButton('clear_output', 'Clear Output')
  )
}

create_main_panel_psc <- function() {
  mainPanel(
    plotOutput('plot', height='500px'),
    verbatimTextOutput('console'),
    hr(),
    h4("Summary"),
    verbatimTextOutput('accumulated_summary')
  )
}


# DATA PROCESSING HELPER FUNCTIONS


load_uploaded_data <- function(file_path, file_ext) {
  if (tolower(file_ext) == 'csv') {
    read.csv(file_path)
  } else {
    readxl::read_excel(file_path)
  }
}

downsample_data <- function(data_col, ds) {
  if (ds > 1) {
    data_col[seq(1, length(data_col), by=ds)]
  } else {
    data_col
  }
}

clean_column_names <- function(df) {
  names(df) <- gsub('^A\\d+$', 'A', names(df))
  names(df) <- gsub('^area\\d+$', 'area', names(df))
  names(df) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df))
  names(df) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df))
  names(df) <- gsub('half_width', 'half width', names(df))
  df
}

parse_comma_separated <- function(text_input) {
  if (nchar(text_input) > 0) {
    as.numeric(unlist(strsplit(text_input, ',')))
  } else {
    NULL
  }
}

# ANALYSIS HELPER FUNCTIONS
calculate_tmax <- function(y, N, dt, stim_time, baseline, smooth, user_tmax, y_abline, xbar, ybar, xbar_lab, ybar_lab) {
  png(tempfile())
  tmax_value <- determine_tmax2(
    y = y, N = N, dt = dt, 
    stimulation_time = stim_time, 
    baseline = baseline, smooth = smooth,
    # tmax = if (!is.na(user_tmax)) user_tmax else NULL,
    tmax = if (length(user_tmax) > 0 && !is.na(user_tmax)) user_tmax else NULL,
    y_abline = y_abline, xbar = xbar, ybar = ybar,
    xbar_lab = xbar_lab, ybar_lab = ybar_lab
  )
  dev.off()
  tmax_value
}

run_psc_analysis <- function(y, input, tmax_value) {
  dt <- as.numeric(input$dt) * as.numeric(input$ds)
  x <- seq(0, (length(y) - 1) * dt, by=dt)
  x_limit <- tmax_value
  adjusted_response <- y[x < x_limit]
  
  func <- switch(input$func,
    'product1N'=product1N, 'product2N'=product2N, 
    'product3N'=product3N, product1N)
  
  if (!input$sequential_fit) {
    result <- nFIT(
      response=adjusted_response, n=as.numeric(input$n), 
      N=as.numeric(input$N), IEI=as.numeric(input$IEI), 
      dt=dt, func=func, method=input$method, 
      weight_method=input$weight_method,
      MLEsettings=list(
        iter=as.numeric(input$iter), 
        metropolis.scale=as.numeric(input$metropolis_scale), 
        fit.attempts=as.numeric(input$fit_attempts), 
        RWm=input$RWm
      ),
      stimulation_time=as.numeric(input$stimulation_time), 
      baseline=as.numeric(input$baseline), 
      filter=input$filter, fc=as.numeric(input$fc),
      interval=c(as.numeric(input$interval_min), as.numeric(input$interval_max)),
      fast.decay.limit=parse_comma_separated(input$fast_decay_limit), 
      fast.constraint=input$fast_constraint,
      fast.constraint.method=input$fast_constraint_method, 
      first.delay.constraint=input$first_delay_constraint,
      lower=parse_comma_separated(input$lower), 
      upper=parse_comma_separated(input$upper), 
      latency.limit=parse_comma_separated(input$latency_limit),
      return.output=TRUE, show.plot=FALSE, 
      half_width_fit_limit=as.numeric(input$half_width_fit_limit),
      dp=as.numeric(input$dp), height=5, width=5, 
      seed=as.numeric(input$seed)
    )
    result$traces <- traces_fun2(
      y=y, fits=result$fits, dt=dt, 
      N=as.numeric(input$N), IEI=as.numeric(input$IEI),
      stimulation_time=as.numeric(input$stimulation_time), 
      baseline=as.numeric(input$baseline), func=func,
      filter=input$filter, fc=as.numeric(input$fc)
    )
  } else {
    result <- nFIT_sequential(
      response=y, n=as.numeric(input$n), dt=dt, func=func, 
      method=input$method, weight_method=input$weight_method,
      stimulation_time=as.numeric(input$stimulation_time), 
      baseline=as.numeric(input$baseline), 
      fit.limits=as.numeric(input$userTmax),
      fast.decay.limit=parse_comma_separated(input$fast_decay_limit), 
      fast.constraint=input$fast_constraint,
      fast.constraint.method=input$fast_constraint_method, 
      first.delay.constraint=input$first_delay_constraint,
      latency.limit=parse_comma_separated(input$latency_limit), 
      lower=parse_comma_separated(input$lower), 
      upper=parse_comma_separated(input$upper), 
      filter=input$filter, fc=as.numeric(input$fc), 
      interval=c(as.numeric(input$interval_min), as.numeric(input$interval_max)),
      MLEsettings=list(
        iter=as.numeric(input$iter), 
        metropolis.scale=as.numeric(input$metropolis_scale), 
        fit.attempts=as.numeric(input$fit_attempts), 
        RWm=input$RWm
      ),
      MLE.method=input$method, 
      half_width_fit_limit=as.numeric(input$half_width_fit_limit),
      dp=as.numeric(input$dp), lwd=as.numeric(input$lwd), 
      xlab='', ylab='', width=5, height=5,
      return.output=TRUE, show.output=TRUE, show.plot=TRUE, 
      seed=as.numeric(input$seed)
    )
  }
  result
}


# DOWNLOAD HELPER FUNCTIONS
create_summary_dataframe <- function(accumulated_results) {
  all_cols <- list()
  
  for (result in accumulated_results) {
    col_name <- result$column
    df <- result$output
    
    for (row_idx in 1:nrow(df)) {
      row_data <- df[row_idx, , drop = FALSE]
      component <- if (!is.null(rownames(row_data)) && rownames(row_data)[1] != as.character(row_idx)) {
        rownames(row_data)[1]
      } else {
        row_idx
      }
      
      for (col_idx in 1:ncol(row_data)) {
        col_label <- names(row_data)[col_idx]
        value <- row_data[1, col_idx]
        new_col_name <- paste0(col_label, component)
        
        if (is.null(all_cols[[col_name]])) {
          all_cols[[col_name]] <- list(Experiment = col_name)
        }
        all_cols[[col_name]][[new_col_name]] <- value
      }
    }
  }
  
  summary_df <- do.call(rbind, lapply(all_cols, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  rownames(summary_df) <- NULL
  clean_column_names(summary_df)
}

create_metadata_list <- function(input, col_name=NULL) {
  list(
    data_col=if (!is.null(col_name)) col_name else input$data_col,
    dt=input$dt, stimulation_time=input$stimulation_time,
    baseline=input$baseline, n=input$n, y_abline=input$y_abline,
    func=input$func, ds=input$ds, userTmax=input$userTmax,
    fast_constraint=input$fast_constraint, N=input$N, IEI=input$IEI,
    smooth=input$smooth, method=input$method, 
    weight_method=input$weight_method, sequential_fit=input$sequential_fit,
    interval_min=input$interval_min, interval_max=input$interval_max,
    lower=input$lower, upper=input$upper, latency_limit=input$latency_limit,
    iter=input$iter, metropolis_scale=input$metropolis_scale,
    fit_attempts=input$fit_attempts, RWm=input$RWm,
    filter=input$filter, fc=input$fc, 
    half_width_fit_limit=input$half_width_fit_limit,
    seed=input$seed, dp=input$dp,
    fast_constraint_method=input$fast_constraint_method,
    fast_decay_limit=input$fast_decay_limit,
    first_delay_constraint=input$first_delay_constraint
  )
}

write_metadata_sheet <- function(wb, sheet_name, metadata, col_name) {
  metadata_labels <- c(
    'Data column:', 'dt (ms):', 'Stimulation Time:', 'Baseline:', 'n:', 
    'Fit cutoff:', 'Function:', 'Downsample Factor:', 
    'User maximum time for fit:', 'Add fast constraint:', 'N:', 'IEI:', 
    'Smooth:', 'Method:', 'Weighting:', 'Sequential Fit:', 'Min interval:', 
    'Max interval:', 'Lower bounds (comma-separated):', 
    'Upper bounds (comma-separated):', 'Latency limit:', 'MLE Iterations:', 
    'Metropolis Scale:', 'Fit Attempts:', 'Random Walk Metropolis:', 
    'Filter:', 'Filter cutoff (Hz):', 'Half-width fit limit:', 'Seed:', 
    'Decimal points:', 'Fast constraint method:', 'Fast decay limit(s):', 
    'First delay constraint:'
  )
  
  metadata_values <- list(
    col_name, metadata$dt, metadata$stimulation_time, metadata$baseline,
    metadata$n, metadata$y_abline, metadata$func, metadata$ds, 
    metadata$userTmax, metadata$fast_constraint, metadata$N, metadata$IEI,
    metadata$smooth, metadata$method, metadata$weight_method, 
    metadata$sequential_fit, metadata$interval_min, metadata$interval_max,
    metadata$lower, metadata$upper, metadata$latency_limit, metadata$iter,
    metadata$metropolis_scale, metadata$fit_attempts, metadata$RWm,
    metadata$filter, metadata$fc, NA, metadata$seed, metadata$dp,
    metadata$fast_constraint_method, metadata$fast_decay_limit,
    metadata$first_delay_constraint
  )
  
  numeric_labels <- c(
    'dt (ms):', 'Stimulation Time:', 'Baseline:', 'n:', 'Fit cutoff:',
    'Downsample Factor:', 'User maximum time for fit:', 'N:', 'IEI:', 
    'Smooth:', 'Min interval:', 'Max interval:', 'Latency limit:', 
    'MLE Iterations:', 'Metropolis Scale:', 'Fit Attempts:', 
    'Filter cutoff (Hz):', 'Half-width fit limit:', 'Seed:', 'Decimal points:'
  )
  
  logical_labels <- c(
    'Add fast constraint:', 'Sequential Fit:', 'Random Walk Metropolis:', 
    'Filter:', 'First delay constraint:'
  )
  
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, c("Parameter", "Value"), 
                     startRow = 1, startCol = 1, colNames = FALSE)
  
  for (i in seq_along(metadata_labels)) {
    lbl <- metadata_labels[i]
    val <- metadata_values[[i]]
    openxlsx::writeData(wb, sheet_name, lbl, 
                       startRow = i+1, startCol = 1, colNames = FALSE)
    if (lbl %in% numeric_labels) {
      openxlsx::writeData(wb, sheet_name, as.numeric(val), 
                         startRow = i+1, startCol = 2, colNames = FALSE)
    } else if (lbl %in% logical_labels) {
      openxlsx::writeData(wb, sheet_name, as.logical(val), 
                         startRow = i+1, startCol = 2, colNames = FALSE)
    } else {
      openxlsx::writeData(wb, sheet_name, val, 
                         startRow = i+1, startCol = 2, colNames = FALSE)
    }
  }
}

# STATE MANAGEMENT FUNCTIONS
clear_psc_state <- function(state) {
  state$response <- NULL
  state$analysis <- NULL
}

# MAIN FUNCTION - analysePSC (Modular Version)

analysePSC <- function() {
  ui <- fluidPage(
    create_dark_mode_css_psc(),
    add_busy_spinner(spin = "fading-circle", position = "top-right", 
                    color = "#3c8dbc", height = "60px", width = "60px"),
    titlePanel('PSC Analysis'),
    sidebarLayout(create_sidebar_panel_psc(), create_main_panel_psc())
  )

  server <- function(input, output, session) {
    state <- reactiveValues(
      response=NULL, analysis=NULL, accumulated_results=list()
    )

    baseline_debounced <- debounce(reactive(input$baseline), 800)
    stimulation_time_debounced <- debounce(reactive(input$stimulation_time), 800)
    dt_debounced <- debounce(reactive(input$dt), 800)
    smooth_debounced <- debounce(reactive(input$smooth), 800)

    uploaded_data <- reactive({
      req(input$file)
      ext <- tools::file_ext(input$file$name)
      load_uploaded_data(input$file$datapath, ext)
    })
    
    output$column_selector <- renderUI({
      req(uploaded_data())
      selectInput('data_col', 'Select Column to Analyse', choices=names(uploaded_data()))
    })
    
    observeEvent(input$run_initial, {
      req(uploaded_data(), input$data_col)
      clear_psc_state(state)
      
      data_col <- uploaded_data()[[input$data_col]]
      ds <- as.numeric(input$ds)
      state$response <- downsample_data(data_col, ds)
      
      dt <- as.numeric(input$dt) * ds
      adjusted_tmax <- calculate_tmax(
        state$response, as.numeric(input$N), dt, 
        as.numeric(input$stimulation_time), as.numeric(input$baseline),
        as.numeric(input$smooth), NULL, as.numeric(input$y_abline),
        as.numeric(input$xbar), as.numeric(input$ybar),
        input$xbar_lab, input$ybar_lab
      )
      
      displayed_tmax <- adjusted_tmax - as.numeric(input$stimulation_time) + as.numeric(input$baseline)
      updateNumericInput(session, "userTmax", value = displayed_tmax)
    })
    
    output$accumulated_summary <- renderPrint({
      if (length(state$accumulated_results) == 0) {
        cat("No accumulated results yet. Analyze columns and click 'Add to Results'.\n")
      } else {
        cat(paste("Analysed:", length(state$accumulated_results), "\n\n"))
        for (result in state$accumulated_results) {
          cat("Experiment:", result$column, "\n")
          print(clean_column_names(result$output))
          cat("\n")
        }
      }
    })

    observeEvent(input$add_result, {
      req(state$analysis)
      col_name <- req(input$data_col)
      
      column_result <- list(
        column = col_name,
        output = state$analysis$output,
        traces = state$analysis$traces,
        AIC = state$analysis$AIC,
        BIC = state$analysis$BIC,
        model_message = state$analysis$model.message,
        metadata = create_metadata_list(input, col_name)
      )
      
      existing_cols <- sapply(state$accumulated_results, function(x) x$column)
      existing_idx <- which(existing_cols == col_name)
      
      if (length(existing_idx) > 0) {
        state$accumulated_results[[existing_idx[1]]] <- column_result
        showNotification(paste0("Updated: ", col_name), type = "message", duration = 3)
      } else {
        state$accumulated_results[[length(state$accumulated_results) + 1]] <- column_result
        showNotification(paste0("Added: ", col_name, " (Total: ", length(state$accumulated_results), ")"), 
                        type = "message", duration = 3)
      }
    })

    observeEvent(input$clear_results, {
      state$accumulated_results <- list()
      showNotification("All accumulated results cleared", type = "warning", duration = 3)
    })

    observeEvent(input$ds, {
      req(uploaded_data(), input$data_col)
      if (!is.null(state$response)) {
        data_col <- uploaded_data()[[input$data_col]]
        state$response <- downsample_data(data_col, as.numeric(input$ds))
        state$analysis <- NULL
      }
    }, ignoreInit=TRUE)
    
    observeEvent(list(input$func, input$N, input$IEI), {
      if (!is.null(state$analysis) && !is.null(state$response)) {
        state$analysis <- NULL
        showNotification("Model changed. Please re-run analysis.", type = "warning", duration = 3)
      }
    }, ignoreInit = TRUE)

    output$plot <- renderPlot({
      req(state$response)
      req(baseline_debounced(), stimulation_time_debounced())
      req(baseline_debounced() > 0, stimulation_time_debounced() > 0)
      
      dt <- as.numeric(dt_debounced()) * as.numeric(input$ds)
      
      if (is.null(state$analysis)) {
        determine_tmax2(
          y=state$response, N=as.numeric(input$N), dt=dt, 
          stimulation_time=as.numeric(stimulation_time_debounced()),
          baseline=as.numeric(baseline_debounced()), 
          smooth=as.numeric(smooth_debounced()),
          lwd=as.numeric(input$lwd), cex=1, tmax=NULL, 
          y_abline=as.numeric(input$y_abline),
          xbar=as.numeric(input$xbar), ybar=as.numeric(input$ybar),
          xbar_lab=input$xbar_lab, ybar_lab=input$ybar_lab
        )
      } else {
        req(state$analysis$traces)
        func <- switch(input$func, 'product1N'=product1N, 
                      'product2N'=product2N, 'product3N'=product3N, product1N)
        xlim_vals <- parse_comma_separated(input$xlim)
        traces <- state$analysis$traces
        
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }
        
        drawPlot2(traces=traces, func=func, lwd=as.numeric(input$lwd),
                 filter=input$filter, xbar=as.numeric(input$xbar), 
                 ybar=as.numeric(input$ybar),
                 xbar_lab=input$xbar_lab, ybar_lab=input$ybar_lab)
      }
    })
    
    observeEvent(input$run_main, {
      req(state$response)
      dt <- as.numeric(input$dt) * as.numeric(input$ds)
      
      y <- state$response
      if (any(is.na(y))) y <- y[!is.na(y)]
      
      tmax_value <- calculate_tmax(
        y, as.numeric(input$N), dt, 
        as.numeric(input$stimulation_time), as.numeric(input$baseline),
        as.numeric(input$smooth), 
        if (!is.na(as.numeric(input$userTmax))) as.numeric(input$userTmax) else NULL,
        as.numeric(input$y_abline), as.numeric(input$xbar), 
        as.numeric(input$ybar), input$xbar_lab, input$ybar_lab
      )
      
      # tryCatch to handle fit failures
      tryCatch({
        state$analysis <- run_psc_analysis(y, input, tmax_value)
        showNotification("Analysis complete!", type = "message", duration = 3)
      }, error = function(e) {
        # Don't update state$analysis if error to preserve previous results
        showNotification(
          paste("Analysis failed:", e$message, "- Try adjusting fit parameters or constraints"),
          type = "error",
          duration = 10
        )
      })
    })
    
    observeEvent(input$clear_output, { state$analysis <- NULL })
    
    output$console <- renderPrint({
      if (!is.null(state$analysis)) {
        print(clean_column_names(state$analysis$output))
      } else {
        cat('No analysis output performed')
      }
    })
    
    output$download_output <- downloadHandler(
      filename=function() {
        req(input$file)
        paste0(tools::file_path_sans_ext(basename(input$file$name)),
              "_", input$data_col, "_PSC_analysis.RData")
      },
      content=function(file) {
        results <- list(
          analysis=state$analysis,
          metadata=create_metadata_list(input)
        )
        save(results, file=file)
      }
    )

    output$download_svg <- downloadHandler(
      filename = function() paste0('PSC_plot_', Sys.Date(), '.svg'),
      content = function(file) {
        req(state$analysis)
        func <- switch(input$func, 'product1N'=product1N, 
                      'product2N'=product2N, 'product3N'=product3N, product1N)
        traces <- state$analysis$traces
        xlim_vals <- parse_comma_separated(input$xlim)
        
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }
        
        svg(filename = file, width = 7, height = 5)
        drawPlot2(traces=traces, func=func, lwd=as.numeric(input$lwd),
                 filter=input$filter, xbar=as.numeric(input$xbar), 
                 ybar=as.numeric(input$ybar),
                 xbar_lab=input$xbar_lab, ybar_lab=input$ybar_lab)
        dev.off()
      }
    )

    output$download_xlsx <- downloadHandler(
      filename = function() {
        if (length(state$accumulated_results) > 0) {
          paste0(tools::file_path_sans_ext(basename(input$file$name)), "_PSC_analyses.zip")
        } else {
          paste0(tools::file_path_sans_ext(basename(input$file$name)), 
                "_", input$data_col, "_PSC_analysis.xlsx")
        }
      },
      content = function(file) {
        if (length(state$accumulated_results) > 0) {
          temp_dir <- tempdir()
          files_to_zip <- c()
          
          # Summary file
          summary_file <- file.path(temp_dir, "summary.xlsx")
          wb_summary <- openxlsx::createWorkbook()
          summary_df <- create_summary_dataframe(state$accumulated_results)
          openxlsx::addWorksheet(wb_summary, "summary")
          openxlsx::writeData(wb_summary, "summary", summary_df)
          openxlsx::saveWorkbook(wb_summary, summary_file, overwrite = TRUE)
          files_to_zip <- c(files_to_zip, summary_file)
          
          # Individual files
          for (result in state$accumulated_results) {
            col_name <- result$column
            safe_name <- gsub("[:\\\\/?*\\[\\]]", "_", col_name)
            individual_file <- file.path(temp_dir, paste0(safe_name, "_PSC_analysis.xlsx"))
            wb_individual <- openxlsx::createWorkbook()
            
            openxlsx::addWorksheet(wb_individual, "output")
            openxlsx::writeData(wb_individual, "output", result$output)
            
            openxlsx::addWorksheet(wb_individual, "traces")
            openxlsx::writeData(wb_individual, "traces", result$traces)
            
            openxlsx::addWorksheet(wb_individual, "fit criterion")
            openxlsx::writeData(wb_individual, "fit criterion", 
                               data.frame(AIC = result$AIC, BIC = result$BIC))
            
            openxlsx::addWorksheet(wb_individual, "model message")
            openxlsx::writeData(wb_individual, "model message", 
                               data.frame(message = result$model_message))
            
            write_metadata_sheet(wb_individual, "metadata", result$metadata, col_name)
            openxlsx::saveWorkbook(wb_individual, individual_file, overwrite = TRUE)
            files_to_zip <- c(files_to_zip, individual_file)
          }
          
          oldwd <- getwd()
          setwd(temp_dir)
          zip::zip(zipfile = file, files = basename(files_to_zip))
          setwd(oldwd)
          
        } else {
          req(state$analysis)
          wb <- openxlsx::createWorkbook()
          
          openxlsx::addWorksheet(wb, "output")
          openxlsx::writeData(wb, "output", state$analysis$output)
          
          openxlsx::addWorksheet(wb, "traces")
          openxlsx::writeData(wb, "traces", state$analysis$traces)
          
          openxlsx::addWorksheet(wb, "fit criterion")
          openxlsx::writeData(wb, "fit criterion", 
                             data.frame(AIC = state$analysis$AIC, BIC = state$analysis$BIC))
          
          openxlsx::addWorksheet(wb, "model message")
          openxlsx::writeData(wb, "model message", 
                             data.frame(message = state$analysis$model.message))
          
          write_metadata_sheet(wb, "metadata", create_metadata_list(input), input$data_col)
          openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        }
      }
    )
  }

  shinyApp(ui=ui, server=server)
}

analysePSC2 <- function() {

  ui <- fluidPage(
    # Add dark mode CSS
    tags$head(
      tags$style(HTML("
        /* Light mode (default) */
        body {
          background-color: white;
          color: black;
        }
        
        /* Dark mode - auto-detect browser preference */
        @media (prefers-color-scheme: dark) {
          body {
            background-color: #1e1e1e;
            color: #e0e0e0;
          }
          
          .well {
            background-color: #2d2d2d;
            border-color: #444;
          }
          
          /* Regular input fields - keep dark */
          .form-control {
            background-color: #2d2d2d;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* Numeric and text inputs - keep dark */
          input[type='number'],
          input[type='text'] {
            background-color: #2d2d2d !important;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* DROPDOWN MENUS ONLY - WHITE boxes with dark gray text */
          .selectize-input, .selectize-dropdown {
            background-color: #ffffff !important;
            color: #666666 !important;
            border: none !important;  /* NO BORDER AT ALL */
          }

          /* Dropdown options */
          .selectize-dropdown .option {
            background-color: #ffffff;
            color: #666666;
          }

          .selectize-dropdown .option:hover {
            background-color: #f0f0f0;
            color: #333333;
          }

          /* Selected option in dropdown */
          .selectize-dropdown .selected,
          .selectize-dropdown .active {
            background-color: #e0e0e0;
            color: #333333;
          }

          /* Selectize item (selected chips) - NO GREY BOX */
          .selectize-input .item {
            background-color: #ffffff !important;  /* WHITE, not grey */
            color: #666666 !important;
            border: none !important;  /* NO border */
          }
          
          /* Labels - also brighter */
          .shiny-input-container label, h4, h3 {
            color: #f0f0f0;
            font-weight: 500;
          }
          
          /* Radio buttons and checkbox labels */
          .radio label, .checkbox label {
            color: #f0f0f0;
          }
          
          /* Tab navigation */
          .nav-tabs {
            border-bottom-color: #444;
          }
          
          .nav-tabs > li > a {
            background-color: #2d2d2d;
            color: #e0e0e0;
            border-color: #444;
          }
          
          .nav-tabs > li.active > a,
          .nav-tabs > li.active > a:hover,
          .nav-tabs > li.active > a:focus {
            background-color: #1e1e1e;
            color: #fff;
            border-color: #444 #444 transparent;
          }
          
          /* Buttons */
          .btn {
            background-color: #3d3d3d;
            color: #ffffff;
            border-color: #555;
          }
          
          .btn:hover {
            background-color: #4d4d4d;
            border-color: #666;
            color: #ffffff;
          }
          
          .btn-default:hover {
            color: #fff;
          }
          
          /* Download buttons */
          .btn-default {
            color: #ffffff;
          }
          
          /* Action buttons - make them stand out more */
          .btn-primary, .action-button {
            background-color: #3c8dbc;
            color: #ffffff;
            border-color: #357ca5;
          }
          
          .btn-primary:hover, .action-button:hover {
            background-color: #4a9dd1;
            border-color: #428bca;
          }
          
          /* Verbatim output */
          pre, code {
            background-color: #1a1a1a;
            color: #f0f0f0;
            border-color: #444;
          }
          
          hr {
            border-top-color: #444;
          }
          
          /* File input styling */
          .btn-file {
            background-color: #3d3d3d;
            color: #ffffff;
          }
          
          /* Progress bar in dark mode */
          .shiny-progress .progress {
            background-color: #2d2d2d;
          }
          
          .shiny-progress .progress-bar {
            background-color: #3c8dbc;
          }
          
          .shiny-progress-notification {
            background-color: #2d2d2d;
            color: #f0f0f0;
            border-color: #444;
          }
          
          /* Focus states for inputs */
          .form-control:focus, .selectize-input.focus {
            border-color: #3c8dbc;
            box-shadow: 0 0 0 0.2rem rgba(60, 141, 188, 0.25);
          }
        }
      "))
    ),

    # Add automatic busy indicator (ADD THIS)
    add_busy_spinner(
      spin = "fading-circle",
      position = "top-right",
      color = "#3c8dbc",
      height = "60px",
      width = "60px"
    ),

    titlePanel('PSC Analysis'),
    sidebarLayout(
      sidebarPanel(
        fileInput('file', 'Upload csv or xlsx', accept=c('.csv', '.xlsx')),
        uiOutput('column_selector'),
        
        tabsetPanel(
          tabPanel('Main Options',
                   numericInput('dt', 'dt (ms):', 0.1),
                   numericInput('stimulation_time', 'Stimulation Time:', 100),
                   numericInput('baseline', 'Baseline:', 50),
                   numericInput('n', 'n:', 30),
                   numericInput('y_abline', 'Fit Cutoff:', 0.1),
                   selectInput('func', 'Function:', choices=c('product1N', 'product2N', 'product3N')),
                   numericInput('ds', 'Downsample Factor:', 1, min=1)
          ),
          tabPanel('Fit Options',
                   numericInput('N', 'N:', 1),
                   numericInput('IEI', 'IEI:', 50),
                   numericInput('smooth', 'Smooth:', 5),
                   selectInput('method', 'Method:', choices=c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')),
                   selectInput('weight_method', 'Weighting:', choices=c('none', '~y_sqrt', '~y')),
                   checkboxInput('sequential_fit', 'Sequential Fit', FALSE),
                   numericInput('interval_min', 'Min Interval:', 0.1),
                   numericInput('interval_max', 'Max Interval:', 0.9),
                   textInput('lower', 'Lower Bounds (comma-separated):', ''),
                   textInput('upper', 'Upper Bounds (comma-separated):', ''),
                   textInput('latency_limit', 'Latency Limit:', '')
          ),
          tabPanel('MLE Settings',
                   numericInput('iter', 'MLE Iterations:', 1000),
                   numericInput('metropolis_scale', 'Metropolis Scale:', 1.5),
                   numericInput('fit_attempts', 'Fit Attempts:', 10),
                   checkboxInput('RWm', 'Random Walk Metropolis', FALSE)
          ),
          tabPanel('Advanced',
                   checkboxInput('filter', 'Filter', FALSE),
                   numericInput('fc', 'Filter Cutoff (Hz):', 1000),
                   numericInput('half_width_fit_limit', 'Half-width Fit Limit:', 500),
                   numericInput('seed', 'Seed:', 42),
                   numericInput('dp', 'Decimal Points:', 3),
                   checkboxInput('fast_constraint', 'Fast Constraint', FALSE),
                   selectInput('fast_constraint_method', 'Fast Constraint Method:', choices=c('rise', 'peak')),
                   textInput('fast_decay_limit', 'Fast Decay Limit(s) (comma-separated):', ''),
                   checkboxInput('first_delay_constraint', 'First Delay Constraint', FALSE)
          ),
          tabPanel('Plot Settings',
                   numericInput('lwd', 'Line Width:', 1.2),
                   numericInput('xbar', 'x-bar Length:', 50),
                   numericInput('ybar', 'y-bar Length:', 50),
                   textInput('xbar_lab', 'x-axis Units:', 'ms'),
                   textInput('ybar_lab', 'y-axis Units:', 'pA'),
                   textInput('xlim', 'x limits (e.g., 0,400):', '')
          )
        ),
        
        numericInput('userTmax', 'User Maximum Time for Fit:', NA),
        actionButton('run_initial', 'Run Initial Analysis'),
        actionButton('run_main', 'Run Main Analysis'),
        actionButton('add_result', 'Add to Results'),
        actionButton('clear_results', 'Clear Results'), 
        downloadButton('download_xlsx',  'Download Output (*.xlsx)'),
        downloadButton('download_output', 'Download RData'),
        downloadButton('download_svg', 'Download SVG Plot'),
        actionButton('clear_output', 'Clear Output')
      ),
      mainPanel(
        plotOutput('plot', height='500px'),
        verbatimTextOutput('console'),
        
        hr(),
        h4("Summary"),
        verbatimTextOutput('accumulated_summary')
      )
    )
  )

  server <- function(input, output, session) {
    
    # Reactive values to store data and analysis results.
    state <- reactiveValues(
      response=NULL,
      analysis=NULL,
      accumulated_results=list()
    )
    
    baseline_debounced <- debounce(reactive(input$baseline), 800)
    stimulation_time_debounced <- debounce(reactive(input$stimulation_time), 800)
    dt_debounced <- debounce(reactive(input$dt), 800)
    smooth_debounced <- debounce(reactive(input$smooth), 800)



    # upload file
    uploaded_data <- reactive({
      req(input$file)
      ext <- tools::file_ext(input$file$name)
      if (tolower(ext) == 'csv') {
        read.csv(input$file$datapath)
      } else {
        readxl::read_excel(input$file$datapath)
      }
    })
    
    # update the column selector
    output$column_selector <- renderUI({
      req(uploaded_data())
      selectInput('data_col', 'Select Column to Analyse', choices=names(uploaded_data()))
    })
    
    # # run initial analysis
    # observeEvent(input$run_initial, {
    #   req(uploaded_data(), input$data_col)
    #   # Clear any previous response and analysis.
    #   state$response <- NULL
    #   state$analysis <- NULL
    #   # Extract the column from the uploaded data.
    #   data_col <- uploaded_data()[[input$data_col]]
    #   ds <- as.numeric(input$ds)
    #   if (ds > 1) {
    #     data_col <- data_col[seq(1, length(data_col), by=ds)]
    #   }
    #   state$response <- data_col
    # })
    
    # run initial analysis
    observeEvent(input$run_initial, {
      req(uploaded_data(), input$data_col)
      # Clear any previous response and analysis.
      state$response <- NULL
      state$analysis <- NULL
      # Extract the column from the uploaded data.
      data_col <- uploaded_data()[[input$data_col]]
      ds <- as.numeric(input$ds)
      if (ds > 1) {
        data_col <- data_col[seq(1, length(data_col), by=ds)]
      }
      state$response <- data_col
      
    # Calculate and auto-populate the displayed tmax value (suppress graphics)
    dt <- as.numeric(input$dt) * ds
    png(tempfile())  # Create temporary graphics device (won't display)
    adjusted_tmax <- determine_tmax2(
      y = data_col, 
      N = as.numeric(input$N), 
      dt = dt, 
      stimulation_time = as.numeric(input$stimulation_time), 
      baseline = as.numeric(input$baseline), 
      smooth = as.numeric(input$smooth),
      tmax = NULL, 
      y_abline = as.numeric(input$y_abline),
      xbar = as.numeric(input$xbar), 
      ybar = as.numeric(input$ybar),
      xbar_lab = input$xbar_lab, 
      ybar_lab = input$ybar_lab
    )
    dev.off()  # Close the temporary device
      
      # Convert from adjusted value to displayed value
      displayed_tmax <- adjusted_tmax - as.numeric(input$stimulation_time) + as.numeric(input$baseline)
      updateNumericInput(session, "userTmax", value = displayed_tmax)
    })
  
    # Display accumulated results - show each column's table separately
    output$accumulated_summary <- renderPrint({
      if (length(state$accumulated_results) == 0) {
        cat("No accumulated results yet. Analyze columns and click 'Add to Results'.\n")
      } else {
        cat(paste("Analysed:", length(state$accumulated_results), "\n\n"))
        
        for (result in state$accumulated_results) {
          cat("Experiment:", result$column, "\n")
          
          df_out <- result$output
          # Clean up column names exactly as console does
          if(sum(grepl('^A\\d+$', names(df_out))) == 1) names(df_out)[which(grepl('^A\\d+$', names(df_out)))] <- 'A'
          if(sum(grepl('^area\\d+$', names(df_out)))==1) names(df_out)[which(grepl('^area\\d+$', names(df_out)))] <- 'area'
          names(df_out) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df_out))
          names(df_out) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df_out))
          names(df_out)[names(df_out) == 'half_width'] <- 'half width'
          
          print(df_out)
          cat("\n")
        }
      }
    })

    observeEvent(input$add_result, {
      req(state$analysis)
      
      col_name <- req(input$data_col)
      if (is.null(col_name) || col_name == "") {
        showNotification("Error: No column selected", type = "error", duration = 3)
        return()
      }
      
      # Store COMPLETE analysis info for this column (not just output)
      column_result <- list(
        column = col_name,
        output = state$analysis$output,
        traces = state$analysis$traces,
        AIC = state$analysis$AIC,
        BIC = state$analysis$BIC,
        model_message = state$analysis$model.message,
        metadata = list(
          dt = input$dt,
          stimulation_time = input$stimulation_time,
          baseline = input$baseline,
          n = input$n,
          y_abline = input$y_abline,
          func = input$func,
          ds = input$ds,
          userTmax = input$userTmax,
          fast_constraint = input$fast_constraint,
          N = input$N,
          IEI = input$IEI,
          smooth = input$smooth,
          method = input$method,
          weight_method = input$weight_method,
          sequential_fit = input$sequential_fit,
          interval_min = input$interval_min,
          interval_max = input$interval_max,
          lower = input$lower,
          upper = input$upper,
          latency_limit = input$latency_limit,
          iter = input$iter,
          metropolis_scale = input$metropolis_scale,
          fit_attempts = input$fit_attempts,
          RWm = input$RWm,
          filter = input$filter,
          fc = input$fc,
          half_width_fit_limit = input$half_width_fit_limit,
          seed = input$seed,
          dp = input$dp,
          fast_constraint_method = input$fast_constraint_method,
          fast_decay_limit = input$fast_decay_limit,
          first_delay_constraint = input$first_delay_constraint
        )
      )
      
      # Check if column already exists
      existing_cols <- sapply(state$accumulated_results, function(x) x$column)
      existing_idx <- which(existing_cols == col_name)
      
      if (length(existing_idx) > 0) {
        state$accumulated_results[[existing_idx[1]]] <- column_result
        showNotification(paste0("Updated: ", col_name), type = "message", duration = 3)
      } else {
        state$accumulated_results[[length(state$accumulated_results) + 1]] <- column_result
        showNotification(paste0("Added: ", col_name, " (Total: ", length(state$accumulated_results), ")"), type = "message", duration = 3)
      }
    })

    # Clear accumulated results
    observeEvent(input$clear_results, {
      state$accumulated_results <- list()
      showNotification(
        "All accumulated results cleared",
        type = "warning",
        duration = 3
      )
    })

    # update when downsampled
    observeEvent(input$ds, {
      req(uploaded_data(), input$data_col)
      # Only proceed if a response is already loaded.
      if (!is.null(state$response)) {
        data_col <- uploaded_data()[[input$data_col]]
        ds <- as.numeric(input$ds)
        if (ds > 1) {
          data_col <- data_col[seq(1, length(data_col), by=ds)]
        }
        state$response <- data_col
        # Also clear any analysis result to force re-running the main analysis.
        state$analysis <- NULL
        cat('Downsample factor changed: Updated response with length =', length(data_col), '\n')
      }
    }, ignoreInit=TRUE)
    
    # Clear analysis when function or critical parameters change
    observeEvent(list(input$func, input$N, input$IEI), {
      # Only clear if analysis exists and response is loaded
      if (!is.null(state$analysis) && !is.null(state$response)) {
        state$analysis <- NULL
        showNotification(
          "Model changed. Please re-run analysis.",
          type = "warning",
          duration = 3
        )
      }
    }, ignoreInit = TRUE)

    # plot output
    output$plot <- renderPlot({
      req(state$response)
      req(baseline_debounced(), stimulation_time_debounced())
      req(baseline_debounced() > 0, stimulation_time_debounced() > 0)
      
      # Compute effective dt using the current ds.
      dt <- as.numeric(dt_debounced()) * as.numeric(input$ds)
      lwd <- as.numeric(input$lwd)
      stim_time <- as.numeric(stimulation_time_debounced())
      baseline <- as.numeric(baseline_debounced())
      smooth <- as.numeric(smooth_debounced())
      y_abline <- as.numeric(input$y_abline)
      xbar <- as.numeric(input$xbar)
      ybar <- as.numeric(input$ybar)
      xbar_lab <- input$xbar_lab
      ybar_lab <- input$ybar_lab
      
      if (is.null(state$analysis)) {
        determine_tmax2(y=state$response, N=as.numeric(input$N), dt=dt, 
                        stimulation_time=stim_time, baseline=baseline, smooth=smooth,
                        lwd=lwd, cex=1, tmax=NULL, y_abline=y_abline, 
                        xbar=xbar, ybar=ybar, xbar_lab=xbar_lab, ybar_lab=ybar_lab)
      } else {
        req(state$analysis$traces)
        func <- switch(input$func,
                       'product1N'=product1N,
                       'product2N'=product2N,
                       'product3N'=product3N,
                       product1N)

        xlim_vals <- if (nchar(input$xlim) > 0) as.numeric(unlist(strsplit(input$xlim, ","))) else NULL
        traces <- state$analysis$traces
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }

        drawPlot2(traces=traces, func=func, lwd=lwd,
                  filter=input$filter, xbar=xbar, ybar=ybar,
                  xbar_lab=xbar_lab, ybar_lab=ybar_lab)
      }
    })
    
    # run main analysis
    observeEvent(input$run_main, {
      req(state$response)
      
      dt <- as.numeric(input$dt) * as.numeric(input$ds)
      stim_time <- as.numeric(input$stimulation_time)
      baseline <- as.numeric(input$baseline)
      smooth <- as.numeric(input$smooth)
      n <- as.numeric(input$n)
      N <- as.numeric(input$N)
      IEI <- as.numeric(input$IEI)
      method <- input$method
      weight_method <- input$weight_method
      sequential_fit <- input$sequential_fit
      interval <- c(as.numeric(input$interval_min), as.numeric(input$interval_max))
      lower <- if (nchar(input$lower) > 0) as.numeric(unlist(strsplit(input$lower, ','))) else NULL
      upper <- if (nchar(input$upper) > 0) as.numeric(unlist(strsplit(input$upper, ','))) else NULL
      latency_limit <- if (nchar(input$latency_limit) > 0) as.numeric(unlist(strsplit(input$latency_limit, ','))) else NULL
      iter <- as.numeric(input$iter)
      metropolis_scale <- as.numeric(input$metropolis_scale)
      fit_attempts <- as.numeric(input$fit_attempts)
      RWm <- input$RWm
      fast_decay_limit <- if (nchar(input$fast_decay_limit) > 0) as.numeric(unlist(strsplit(input$fast_decay_limit, ','))) else NULL
      fast_constraint <- input$fast_constraint
      fast_constraint_method <- input$fast_constraint_method
      first_delay_constraint <- input$first_delay_constraint
      dp <- as.numeric(input$dp)
      seed <- as.numeric(input$seed)
      filter_flag <- input$filter
      fc <- as.numeric(input$fc)
      
      y <- state$response
      if (any(is.na(y))) y <- y[!is.na(y)]
      x <- seq(0, (length(y) - 1) * dt, by=dt)
      
      # Suppress graphics and always use determine_tmax2
      png(tempfile())
      tmax_value <- determine_tmax2(
        y = y, 
        N = N, 
        dt = dt, 
        stimulation_time = stim_time, 
        baseline = baseline,
        smooth = smooth, 
        tmax = if (!is.na(as.numeric(input$userTmax))) as.numeric(input$userTmax) else NULL,
        y_abline = as.numeric(input$y_abline),
        xbar = as.numeric(input$xbar), 
        ybar = as.numeric(input$ybar),
        xbar_lab = input$xbar_lab, 
        ybar_lab = input$ybar_lab
      )
      dev.off()
      
      x_limit <- tmax_value

      
      adjusted_response <- y[x < x_limit]

          func <- switch(input$func,
                     'product1N'=product1N,
                     'product2N'=product2N,
                     'product3N'=product3N,
                     product1N)
      
      if (!sequential_fit) {
        result <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func,
                       method=method, weight_method=weight_method,
                       MLEsettings=list(iter=iter, metropolis.scale=metropolis_scale, fit.attempts=fit_attempts, RWm=RWm),
                       stimulation_time=stim_time, baseline=baseline, filter=filter_flag, fc=fc,
                       interval=interval, fast.decay.limit=fast_decay_limit, fast.constraint=fast_constraint,
                       fast.constraint.method=fast_constraint_method, first.delay.constraint=first_delay_constraint,
                       lower=lower, upper=upper, latency.limit=latency_limit,
                       return.output=TRUE, show.plot=FALSE, half_width_fit_limit=as.numeric(input$half_width_fit_limit),
                       dp=dp, height=5, width=5, seed=seed)
        result$traces <- traces_fun2(y=y, fits=result$fits, dt=dt, N=N, IEI=IEI,
                                     stimulation_time=stim_time, baseline=baseline, func=func,
                                     filter=filter_flag, fc=fc)
      } else {
        result <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method,
                                  stimulation_time=stim_time, baseline=baseline, fit.limits=as.numeric(input$userTmax),
                                  fast.decay.limit=fast_decay_limit, fast.constraint=fast_constraint,
                                  fast.constraint.method=fast_constraint_method, first.delay.constraint=first_delay_constraint,
                                  latency.limit=latency_limit, lower=lower, upper=upper, filter=filter_flag, fc=fc, interval=interval,
                                  MLEsettings=list(iter=iter, metropolis.scale=metropolis_scale, fit.attempts=fit_attempts, RWm=RWm),
                                  MLE.method=method, half_width_fit_limit=as.numeric(input$half_width_fit_limit),
                                  dp=dp, lwd=as.numeric(input$lwd), xlab='', ylab='', width=5, height=5,
                                  return.output=TRUE, show.output=TRUE, show.plot=TRUE, seed=seed)
      }
      
      state$analysis <- result
    })
    
    # clear output
    observeEvent(input$clear_output, {
      state$analysis <- NULL
    })
    
    # output to console
    output$console <- renderPrint({
      if (!is.null(state$analysis)) {

        df_out <- state$analysis$output
        if(sum(grepl('^A\\d+$', names(df_out))) == 1) names(df_out)[which(grepl('^A\\d+$', names(df_out)))] <- 'A'
        if(sum(grepl('^area\\d+$', names(df_out)))==1) names(df_out)[which(grepl('^area\\d+$', names(df_out)))] <- 'area'
        names(df_out) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df_out))
        names(df_out) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df_out))
        names(df_out)[names(df_out) == 'half_width'] <- 'half width'

        print(df_out)
      } else {
        cat('No analysis output performed')
      }
    })
    
    # download output
    output$download_output <- downloadHandler(
      filename=function() {
        req(input$file)
        paste0(
          tools::file_path_sans_ext(basename(input$file$name)),
          "_", input$data_col,
          "_PSC_analysis.RData"
        )
      },
      content=function(file) {
        # all settings
        metadata <- list(
          dt=as.numeric(input$dt),
          ds=as.numeric(input$ds),
          stimulation_time=as.numeric(input$stimulation_time),
          baseline=as.numeric(input$baseline),
          n=as.numeric(input$n),
          y_abline=as.numeric(input$y_abline),
          func=input$func,
          N=as.numeric(input$N),
          IEI=as.numeric(input$IEI),
          smooth=as.numeric(input$smooth),
          method=input$method,
          weight_method=input$weight_method,
          sequential_fit=input$sequential_fit,
          interval=c(as.numeric(input$interval_min),
                       as.numeric(input$interval_max)),
          lower=if(nchar(input$lower) > 0)
                    as.numeric(unlist(strsplit(input$lower, ",")))
                  else NULL,
          upper=if(nchar(input$upper) > 0)
                    as.numeric(unlist(strsplit(input$upper, ",")))
                  else NULL,
          latency_limit=if(nchar(input$latency_limit) > 0)
                            as.numeric(unlist(strsplit(input$latency_limit, ",")))
                          else NULL,
          iter=as.numeric(input$iter),
          metropolis_scale=as.numeric(input$metropolis_scale),
          fit_attempts=as.numeric(input$fit_attempts),
          RWm=input$RWm,
          fast_decay_limit=if(nchar(input$fast_decay_limit) > 0) as.numeric(unlist(strsplit(input$fast_decay_limit, ","))) else NULL,
          fast_constraint=input$fast_constraint,
          fast_constraint_method=input$fast_constraint_method,
          first_delay_constraint=input$first_delay_constraint,
          dp=as.numeric(input$dp),
          seed=as.numeric(input$seed),
          filter=input$filter,
          fc=as.numeric(input$fc),
          userTmax=as.numeric(input$userTmax),
          data_col=input$data_col
        )
        
        # save analysis and metadata
        results <- list(
          analysis=state$analysis,
          metadata=metadata
        )
        
        save(results, file=file)
      }
    )

    output$download_svg <- downloadHandler(
      filename = function() {
        paste0('PSC_plot_', Sys.Date(), '.svg')
      },
      content = function(file) {
        req(state$analysis)
        func <- switch(input$func,
                       'product1N' = product1N,
                       'product2N' = product2N,
                       'product3N' = product3N,
                       product1N)

        traces <- state$analysis$traces
        xlim_vals <- if (nchar(input$xlim) > 0) as.numeric(unlist(strsplit(input$xlim, ","))) else NULL
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }

        svg(filename = file, width = 7, height = 5)
        drawPlot2(
          traces = traces,
          func = func,
          lwd = as.numeric(input$lwd),
          filter = input$filter,
          xbar = as.numeric(input$xbar),
          ybar = as.numeric(input$ybar),
          xbar_lab = input$xbar_lab,
          ybar_lab = input$ybar_lab
        )
        dev.off()
      }
    )

    output$download_xlsx <- downloadHandler(
      filename = function() {
        if (length(state$accumulated_results) > 0) {
          paste0(
            tools::file_path_sans_ext(basename(input$file$name)),
            "_PSC_analyses.zip"
          )
        } else {
          paste0(
            tools::file_path_sans_ext(basename(input$file$name)),
            "_", input$data_col,
            "_PSC_analysis.xlsx"
          )
        }
      },
      content = function(file) {
        
        if (length(state$accumulated_results) > 0) {
          # Multi-column download: Create ZIP with summary + individual Excel files
          
          temp_dir <- tempdir()
          files_to_zip <- c()
          
          # 1. Create summary.xlsx
          summary_file <- file.path(temp_dir, "summary.xlsx")
          wb_summary <- openxlsx::createWorkbook()
          
          # Build wide format summary
          all_cols <- list()
          
          for (result in state$accumulated_results) {
            col_name <- result$column
            df <- result$output
            
            for (row_idx in 1:nrow(df)) {
              row_data <- df[row_idx, , drop = FALSE]
              
              if (!is.null(rownames(row_data)) && rownames(row_data)[1] != as.character(row_idx)) {
                component <- rownames(row_data)[1]
              } else {
                component <- row_idx
              }
              
              for (col_idx in 1:ncol(row_data)) {
                col_label <- names(row_data)[col_idx]
                value <- row_data[1, col_idx]
                new_col_name <- paste0(col_label, component)
                
                if (is.null(all_cols[[col_name]])) {
                  all_cols[[col_name]] <- list(Experiment = col_name)
                }
                all_cols[[col_name]][[new_col_name]] <- value
              }
            }
          }
          
          summary_df <- do.call(rbind, lapply(all_cols, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
          rownames(summary_df) <- NULL
          
          names(summary_df) <- gsub("^A\\d+$", "A", names(summary_df))
          names(summary_df) <- gsub("^area\\d+$", "area", names(summary_df))
          names(summary_df) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(summary_df))
          names(summary_df) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(summary_df))
          names(summary_df) <- gsub("half_width", "half width", names(summary_df))
          
          openxlsx::addWorksheet(wb_summary, "summary")
          openxlsx::writeData(wb_summary, "summary", summary_df)
          openxlsx::saveWorkbook(wb_summary, summary_file, overwrite = TRUE)
          files_to_zip <- c(files_to_zip, summary_file)
          
          # 2. Create individual Excel files for EACH experiment
          for (result in state$accumulated_results) {
            col_name <- result$column
            safe_name <- gsub("[:\\\\/?*\\[\\]]", "_", col_name)
            
            individual_file <- file.path(temp_dir, paste0(safe_name, "_PSC_analysis.xlsx"))
            wb_individual <- openxlsx::createWorkbook()
            
            # Output sheet
            openxlsx::addWorksheet(wb_individual, "output")
            openxlsx::writeData(wb_individual, "output", result$output)
            
            # Traces sheet
            openxlsx::addWorksheet(wb_individual, "traces")
            openxlsx::writeData(wb_individual, "traces", result$traces)
            
            # Fit criterion sheet
            openxlsx::addWorksheet(wb_individual, "fit criterion")
            openxlsx::writeData(wb_individual, "fit criterion", 
                               data.frame(AIC = result$AIC, BIC = result$BIC))
            
            # Model message sheet
            openxlsx::addWorksheet(wb_individual, "model message")
            openxlsx::writeData(wb_individual, "model message", 
                               data.frame(message = result$model_message))
            
            # Metadata sheet
            metadata_labels <- c(
              'Data column:','dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:','Function:',
              'Downsample Factor:','User maximum time for fit:','Add fast constraint:','N:','IEI:','Smooth:',
              'Method:','Weighting:','Sequential Fit:','Min interval:','Max interval:',
              'Lower bounds (comma-separated):','Upper bounds (comma-separated):','Latency limit:',
              'MLE Iterations:','Metropolis Scale:','Fit Attempts:','Random Walk Metropolis:',
              'Filter:','Filter cutoff (Hz):','Half-width fit limit:','Seed:','Decimal points:',
              'Fast constraint method:','Fast decay limit(s):','First delay constraint:'
            )
            metadata_values <- list(
              col_name, result$metadata$dt, result$metadata$stimulation_time,
              result$metadata$baseline, result$metadata$n, result$metadata$y_abline,
              result$metadata$func, result$metadata$ds, result$metadata$userTmax,
              result$metadata$fast_constraint, result$metadata$N, result$metadata$IEI,
              result$metadata$smooth, result$metadata$method, result$metadata$weight_method,
              result$metadata$sequential_fit, result$metadata$interval_min, result$metadata$interval_max,
              result$metadata$lower, result$metadata$upper, result$metadata$latency_limit,
              result$metadata$iter, result$metadata$metropolis_scale, result$metadata$fit_attempts,
              result$metadata$RWm, result$metadata$filter, result$metadata$fc, NA,
              result$metadata$seed, result$metadata$dp, result$metadata$fast_constraint_method,
              result$metadata$fast_decay_limit, result$metadata$first_delay_constraint
            )
            
            numeric_labels <- c(
              'dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:',
              'Downsample Factor:','User maximum time for fit:','N:','IEI:','Smooth:',
              'Min interval:','Max interval:','Latency limit:','MLE Iterations:',
              'Metropolis Scale:','Fit Attempts:','Filter cutoff (Hz):',
              'Half-width fit limit:','Seed:','Decimal points:'
            )
            logical_labels <- c(
              'Add fast constraint:','Sequential Fit:','Random Walk Metropolis:','Filter:',
              'First delay constraint:'
            )
            
            openxlsx::addWorksheet(wb_individual, "metadata")
            openxlsx::writeData(wb_individual, "metadata", c("Parameter","Value"), 
                               startRow = 1, startCol = 1, colNames = FALSE)
            for (i in seq_along(metadata_labels)) {
              lbl <- metadata_labels[i]
              val <- metadata_values[[i]]
              openxlsx::writeData(wb_individual, "metadata", lbl, 
                                 startRow = i+1, startCol = 1, colNames = FALSE)
              if (lbl %in% numeric_labels) {
                openxlsx::writeData(wb_individual, "metadata", as.numeric(val), 
                                   startRow = i+1, startCol = 2, colNames = FALSE)
              } else if (lbl %in% logical_labels) {
                openxlsx::writeData(wb_individual, "metadata", as.logical(val), 
                                   startRow = i+1, startCol = 2, colNames = FALSE)
              } else {
                openxlsx::writeData(wb_individual, "metadata", val, 
                                   startRow = i+1, startCol = 2, colNames = FALSE)
              }
            }
            
            openxlsx::saveWorkbook(wb_individual, individual_file, overwrite = TRUE)
            files_to_zip <- c(files_to_zip, individual_file)
          }
          
          # 3. Create ZIP file with proper working directory
          oldwd <- getwd()
          setwd(temp_dir)
          zip::zip(zipfile = file, files = basename(files_to_zip))
          setwd(oldwd)
          
        } else {
          # Single column download
          req(state$analysis)
          
          wb <- openxlsx::createWorkbook()
          
          data_list <- list(
            output = state$analysis$output,
            traces = state$analysis$traces,
            `fit criterion` = data.frame(AIC = state$analysis$AIC, BIC = state$analysis$BIC),
            `model message` = data.frame(message = state$analysis$model.message)
          )
          for (nm in names(data_list)) {
            openxlsx::addWorksheet(wb, nm)
            openxlsx::writeData(wb, nm, data_list[[nm]])
          }

          metadata_labels <- c(
            'Data column:','dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:','Function:',
            'Downsample Factor:','User maximum time for fit:','Add fast constraint:','N:','IEI:','Smooth:',
            'Method:','Weighting:','Sequential Fit:','Min interval:','Max interval:',
            'Lower bounds (comma-separated):','Upper bounds (comma-separated):','Latency limit:',
            'MLE Iterations:','Metropolis Scale:','Fit Attempts:','Random Walk Metropolis:',
            'Filter:','Filter cutoff (Hz):','Half-width fit limit:','Seed:','Decimal points:',
            'Fast constraint method:','Fast decay limit(s):','First delay constraint:'
          )
          metadata_values <- list(
            input$data_col, input$dt, input$stimulation_time, input$baseline,
            input$n, input$y_abline, input$func, input$ds, input$userTmax,
            input$fast_constraint, input$N, input$IEI, input$smooth,
            input$method, input$weight_method, input$sequential_fit,
            input$interval_min, input$interval_max, input$lower, input$upper,
            input$latency_limit, input$iter, input$metropolis_scale,
            input$fit_attempts, input$RWm, input$filter, input$fc,
            input$half_width_fit_limit, input$seed, input$dp,
            input$fast_constraint_method, input$fast_decay_limit,
            input$first_delay_constraint
          )

          numeric_labels <- c(
            'dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:',
            'Downsample Factor:','User maximum time for fit:','N:','IEI:','Smooth:',
            'Min interval:','Max interval:','Latency limit:','MLE Iterations:',
            'Metropolis Scale:','Fit Attempts:','Filter cutoff (Hz):',
            'Half-width fit limit:','Seed:','Decimal points:'
          )
          logical_labels <- c(
            'Add fast constraint:','Sequential Fit:','Random Walk Metropolis:','Filter:',
            'First delay constraint:'
          )

          openxlsx::addWorksheet(wb, "metadata")
          openxlsx::writeData(wb, "metadata", c("Parameter","Value"), 
                             startRow = 1, startCol = 1, colNames = FALSE)
          for (i in seq_along(metadata_labels)) {
            lbl <- metadata_labels[i]
            val <- metadata_values[[i]]
            openxlsx::writeData(wb, "metadata", lbl, 
                               startRow = i+1, startCol = 1, colNames = FALSE)
            if (lbl %in% numeric_labels) {
              openxlsx::writeData(wb, "metadata", as.numeric(val), 
                                 startRow = i+1, startCol = 2, colNames = FALSE)
            } else if (lbl %in% logical_labels) {
              openxlsx::writeData(wb, "metadata", as.logical(val), 
                                 startRow = i+1, startCol = 2, colNames = FALSE)
            } else {
              openxlsx::writeData(wb, "metadata", val, 
                                 startRow = i+1, startCol = 2, colNames = FALSE)
            }
          }
          
          openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        }
      }
    )

    
  }

  shinyApp(ui=ui, server=server)
}


widgetPSCtk <- function() {

  widget_PSC_tk <- function() {
    ## only one tktoplevel() here:
    tt <- tktoplevel()
    tkwm.title(tt, 'I (pA)')

    Tr            <- tclVar('5')
    Td            <- tclVar('25')
    decayUpperVar <- tclVar('0.9')
    decayLowerVar <- tclVar('0.1')

    updatePlot <- function() {
      currentTr   <- as.numeric(tclvalue(Tr))
      currentTd   <- as.numeric(tclvalue(Td))
      decay_range <- c(
        as.numeric(tclvalue(decayUpperVar)),
        as.numeric(tclvalue(decayLowerVar))
      )
      t_max  <- 250
      t_seq  <- seq(0, t_max, length.out = 2000)
      y_comb <- -(exp(-t_seq / currentTd) - exp(-t_seq / currentTr))
      y_Td   <- -exp(-t_seq / currentTd)
      y_Tr   <- -exp(-t_seq / currentTr)
      T1     <- currentTd * currentTr / (currentTd - currentTr)
      y_1mTr <- -(1 - exp(-t_seq / T1))

      peak_i <- which.min(y_comb)
      y_fall <- y_comb[peak_i:length(y_comb)]
      t_fall <- t_seq[peak_i:length(t_seq)]
      min_v  <- min(y_fall)
      up_v   <- decay_range[1] * min_v
      lo_v   <- decay_range[2] * min_v
      idx_up <- which(y_fall >= up_v)[1]
      idx_lo <- which(y_fall >= lo_v)[1]

      if (is.na(idx_lo)) {
        ext_t  <- seq(0, t_max*2, length.out = length(t_seq)*2)
        ext_y  <- -(exp(-ext_t/currentTd) - exp(-ext_t/currentTr))
        y_fall <- ext_y[peak_i:length(ext_y)]
        t_fall <- ext_t[peak_i:length(ext_t)]
        idx_lo <- which(y_fall >= lo_v)[1]
      }

      decay_time <- if (!is.na(idx_up) && !is.na(idx_lo))
        round(t_fall[idx_lo] - t_fall[idx_up], 2) else NA

      up_pct <- round(decay_range[1]*100)
      lo_pct <- round(decay_range[2]*100)
      title2 <- paste('normalised with', up_pct, '-', lo_pct, '% decay time')

      par(mfrow = c(2,1), mar = c(4,4,3,1))
      plot(t_seq, y_comb, type='l', col='grey', lwd=2, axes=FALSE,
           main=expression(e^{-t/tau[decay]} - e^{-t/tau[rise]}),
           xlab='', ylab='F(x)', xlim=c(0,t_max), ylim=c(-1,0))
      lines(t_seq, y_Td, col='indianred', lty=3, lwd=2)
      lines(t_seq, y_Tr, col='slateblue', lty=3, lwd=2)
      legend('bottomright', legend = c(
        expression(e^{-t/tau[decay]} - e^{-t/tau[rise]}),
        expression(e^{-t/tau[decay]}),
        expression(e^{-t/tau[rise]})
      ),
      col = c('grey', 'indianred', 'slateblue'), lty = c(1,3,3), lwd = 2,
      bty = 'n', inset = c(0.02, 0.15))
      axis(2, las=1, tcl=-0.2)

      plot(t_seq, y_comb/abs(min(y_comb)), type='l', col='grey', lwd=2, axes=FALSE,
           main=title2, xlab='time (ms)', ylab='normalised F(x)',
           xlim=c(0,t_max), ylim=c(-1,0))
      lines(t_seq, y_1mTr, col='slateblue', lty=3, lwd=2)
      lines(t_seq, y_Td, col='indianred', lty=3, lwd=2)
      legend('bottomright', legend = c(
        expression((e^{-t/tau[decay]} - e^{-t/tau[rise]}) / abs),
        expression(-(1 - e^{-t/tau[1]})),
        expression(-e^{-t/tau[decay]})
      ),
      col = c('gray','slateblue','indianred'), lty = c(1,3,3), lwd = 2,
      bty = 'n', inset = c(0.02, 0.15))
      if (!is.na(decay_time)) {
        abline(h=-up_v/min_v, col='darkgrey', lty=3)
        abline(h=-lo_v/min_v, col='darkgrey', lty=3)
        text(t_max*0.8, -0.3,
             paste(up_pct, '-', lo_pct, 'decay =', decay_time, 'ms'),
             col='darkgrey', cex=0.8)
      } else {
        text(t_max*0.8, -0.3, 'Decay not fully reached', col='indianred', cex=0.8)
      }
      axis(1, las=1, tcl=-0.2)
      axis(2, las=1, tcl=-0.2)
    }

    # build the UI
    img <- tkrplot(tt, fun=updatePlot, hscale=1.3, vscale=1.3)
    tkpack(img, side='top', expand=TRUE, fill='both')

    decayFrame <- tkframe(tt)
    tkpack(decayFrame, side='top', fill='x', pady=5)
    tkpack(tklabel(decayFrame, text='Decay range (high, low):'), side='left', padx=5)
    highEntry <- tkentry(decayFrame, textvariable=decayUpperVar, width=5)
    lowEntry  <- tkentry(decayFrame, textvariable=decayLowerVar, width=5)
    tkpack(highEntry, side='left', padx=2)
    tkpack(lowEntry,  side='left', padx=2)
    tkbind(highEntry, '<KeyRelease>', function() tkrplot::tkrreplot(img))
    tkbind(lowEntry,  '<KeyRelease>', function() tkrplot::tkrreplot(img))

    tkpack(tklabel(tt, text='\u03C4 rise'),  side='top')
    Tr_slider <- tkscale(tt, from=0.1, to=50, resolution=0.1,
                         showvalue=TRUE, variable=Tr, orient='horizontal',
                         command=function(...) tkrplot::tkrreplot(img))
    tkpack(Tr_slider, fill='x', padx=10, pady=5)

    tkpack(tklabel(tt, text='\u03C4 decay'), side='top')
    Td_slider <- tkscale(tt, from=0.1, to=200, resolution=1,
                         showvalue=TRUE, variable=Td, orient='horizontal',
                         command=function(...) tkrplot::tkrreplot(img))
    tkpack(Td_slider, fill='x', padx=10, pady=5)

    tkfocus(tt)
    tkwait.window(tt)
  }

  widget_PSC_tk()
}


widgetPSC <- function() {

  ui <- fluidPage(
    # Add dark mode CSS
    tags$head(
      tags$style(HTML("
        /* Light mode (default) */
        body {
          background-color: white;
          color: black;
        }
        
        /* Dark mode - auto-detect browser preference */
        @media (prefers-color-scheme: dark) {
          body {
            background-color: #1e1e1e;
            color: #e0e0e0;
          }
          
          .well {
            background-color: #2d2d2d;
            border-color: #444;
          }
          
          /* Regular input fields - keep dark */
          .form-control {
            background-color: #2d2d2d;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* Numeric and text inputs - keep dark */
          input[type='number'],
          input[type='text'] {
            background-color: #2d2d2d !important;
            color: #c0c0c0 !important;
            border: 1px solid #555 !important;
          }

          /* DROPDOWN MENUS ONLY - WHITE boxes with dark gray text, NO BORDER */
          .selectize-input, .selectize-dropdown {
            background-color: #ffffff !important;
            color: #666666 !important;
            border: none !important;
          }

          /* Dropdown options */
          .selectize-dropdown .option {
            background-color: #ffffff;
            color: #666666;
          }

          .selectize-dropdown .option:hover {
            background-color: #f0f0f0;
            color: #333333;
          }

          /* Selected option in dropdown */
          .selectize-dropdown .selected,
          .selectize-dropdown .active {
            background-color: #e0e0e0;
            color: #333333;
          }

          /* Selectize item (selected chips) - NO GREY BOX */
          .selectize-input .item {
            background-color: #ffffff !important;
            color: #666666 !important;
            border: none !important;
          }
          
          /* Labels - also brighter */
          .shiny-input-container label, h4, h3 {
            color: #f0f0f0;
            font-weight: 500;
          }
          
          /* Radio buttons and checkbox labels */
          .radio label, .checkbox label {
            color: #f0f0f0;
          }
          
          /* Tab navigation */
          .nav-tabs {
            border-bottom-color: #444;
          }
          
          .nav-tabs > li > a {
            background-color: #2d2d2d;
            color: #e0e0e0;
            border-color: #444;
          }
          
          .nav-tabs > li.active > a,
          .nav-tabs > li.active > a:hover,
          .nav-tabs > li.active > a:focus {
            background-color: #1e1e1e;
            color: #fff;
            border-color: #444 #444 transparent;
          }
          
          /* Buttons */
          .btn {
            background-color: #3d3d3d;
            color: #ffffff;
            border-color: #555;
          }
          
          .btn:hover {
            background-color: #4d4d4d;
            border-color: #666;
            color: #ffffff;
          }
          
          .btn-default:hover {
            color: #fff;
          }
          
          /* Download buttons */
          .btn-default {
            color: #ffffff;
          }
          
          /* Action buttons - make them stand out more */
          .btn-primary, .action-button {
            background-color: #3c8dbc;
            color: #ffffff;
            border-color: #357ca5;
          }
          
          .btn-primary:hover, .action-button:hover {
            background-color: #4a9dd1;
            border-color: #428bca;
          }
          
          /* Verbatim output */
          pre, code {
            background-color: #1a1a1a;
            color: #f0f0f0;
            border-color: #444;
          }
          
          hr {
            border-top-color: #444;
          }
          
          /* File input styling */
          .btn-file {
            background-color: #3d3d3d;
            color: #ffffff;
          }
          
          /* Progress bar in dark mode */
          .shiny-progress .progress {
            background-color: #2d2d2d;
          }
          
          .shiny-progress .progress-bar {
            background-color: #3c8dbc;
          }
          
          .shiny-progress-notification {
            background-color: #2d2d2d;
            color: #f0f0f0;
            border-color: #444;
          }
          
          /* Focus states for inputs */
          .form-control:focus, .selectize-input.focus {
            border-color: #3c8dbc;
            box-shadow: 0 0 0 0.2rem rgba(60, 141, 188, 0.25);
          }
          
          /* Slider styling */
          .irs--shiny .irs-bar {
            background: #3c8dbc;
            border-color: #357ca5;
          }
          
          .irs--shiny .irs-handle {
            background: #3d3d3d;
            border-color: #555;
          }
          
          .irs--shiny .irs-single,
          .irs--shiny .irs-from,
          .irs--shiny .irs-to {
            background: #3d3d3d;
            color: #ffffff;
          }
          
          .irs--shiny .irs-grid-text {
            color: #999;
          }
          
          .irs--shiny .irs-line {
            background: #2d2d2d;
            border-color: #555;
          }
        }
      "))
    ),
    
    titlePanel(
      HTML("Interactive Graphs:<br>exp(-t/&tau;<sub>decay</sub>) - exp(-t/&tau;<sub>rise</sub>)")
    ),
    sidebarLayout(
      sidebarPanel(
        numericInput(
          "decayUpper", 
          "Decay upper (fraction)", 
          value = 0.9, min = 0, max = 1, step = 0.01
        ),
        numericInput(
          "decayLower", 
          "Decay lower (fraction)", 
          value = 0.1, min = 0, max = 1, step = 0.001
        ),
        sliderInput(
          "Tr", 
          HTML("&tau;<sub>rise</sub>"), 
          min = 0.1, max = 50, value = 5, step = 0.1
        ),
        sliderInput(
          "Td", 
          HTML("&tau;<sub>decay</sub>"), 
          min = 0.101, max = 200, value = 25, step = 1
        )
      ),
      mainPanel(
        plotlyOutput("upperPlot"),
        plotlyOutput("lowerPlot")
      )
    )
  )

  server <- function(input, output, session) {
    
    # Throttle slider inputs - updates every 150ms while moving
    Tr_throttled <- throttle(reactive(input$Tr), 150)
    Td_throttled <- throttle(reactive(input$Td), 150)
    
    # ensure Td > Tr
    observeEvent(input$Tr, {
      updateSliderInput(
        session, "Td",
        min   = input$Tr + 0.001,
        value = max(input$Td, input$Tr + 0.001)
      )
    })
    
    output$upperPlot <- renderPlotly({
      t_max  <- 250
      t_seq  <- seq(0, t_max, length.out = 300)
      y_combined <- -(exp(-t_seq / Td_throttled()) - exp(-t_seq / Tr_throttled()))
      y_Td <- -exp(-t_seq / Td_throttled())
      y_Tr <- -exp(-t_seq / Tr_throttled())
      
      plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = t_seq, y = y_combined, name = 'exp(-t/τ_decay) - exp(-t/τ_rise)',
                  line = list(color = 'grey', width = 2)) %>%
        add_trace(x = t_seq, y = y_Td, name = 'exp(-t/τ_decay)',
                  line = list(color = 'indianred', dash = 'dot', width = 2)) %>%
        add_trace(x = t_seq, y = y_Tr, name = 'exp(-t/τ_rise)',
                  line = list(color = 'slateblue', dash = 'dot', width = 2)) %>%
        layout(
          title = 'exp(-t/τ_decay) - exp(-t/τ_rise)',
          xaxis = list(title = '', showgrid = FALSE, showline = TRUE, ticks = 'outside', zeroline = FALSE),
          yaxis = list(title = 'F(x)', showgrid = FALSE, showline = TRUE, ticks = 'outside', zeroline = FALSE),
          legend = list(x = 0.98, y = 0.02, xanchor = 'right', yanchor = 'bottom')
        )
    })
    
    output$lowerPlot <- renderPlotly({
      t_max    <- 250
      t_seq    <- seq(0, t_max, length.out = 300)
      y_comb   <- -(exp(-t_seq / Td_throttled()) - exp(-t_seq / Tr_throttled()))
      norm     <- y_comb / abs(min(y_comb))
      T1       <- Td_throttled() * Tr_throttled() / (Td_throttled() - Tr_throttled())
      y_1mTr   <- -(1 - exp(-t_seq / T1))
      y_Td     <- -exp(-t_seq / Td_throttled())
      
      peak_i   <- which.min(y_comb)
      y_fall   <- y_comb[peak_i:length(y_comb)]
      t_fall   <- t_seq[peak_i:length(t_seq)]
      min_val  <- min(y_fall)
      
      du <- input$decayUpper * min_val
      dl <- input$decayLower * min_val
      idx_up   <- which(y_fall >= du)[1]
      idx_lo   <- which(y_fall >= dl)[1]
      
      if (is.na(idx_lo)) {
        ext_t <- seq(0, t_max*2, length.out = length(t_seq)*2)
        ext_y <- -(exp(-ext_t/Td_throttled()) - exp(-ext_t/Tr_throttled()))
        y_fall <- ext_y[peak_i:length(ext_y)]
        t_fall <- ext_t[peak_i:length(ext_t)]
        idx_lo <- which(y_fall >= dl)[1]
      }
      
      decay_time <- if (!is.na(idx_up) && !is.na(idx_lo))
        round(t_fall[idx_lo] - t_fall[idx_up], 2) else NA
      
      up_pct <- round(input$decayUpper * 100)
      lo_pct <- round(input$decayLower * 100)
      title2 <- paste("normalised with", up_pct, "-", lo_pct, "% decay time")
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = t_seq, y = norm, name = '(exp(-t/τ_decay) - exp(-t/τ_rise)) / abs',
                  line = list(color = 'grey', width = 2)) %>%
        add_trace(x = t_seq, y = y_1mTr, name = '-(1 - exp(-t/τ1))',
                  line = list(color = 'slateblue', dash = 'dot', width = 2)) %>%
        add_trace(x = t_seq, y = y_Td, name = '-exp(-t/τ_decay)',
                  line = list(color = 'indianred', dash = 'dot', width = 2)) %>%
        layout(title = title2,
               xaxis = list(title = 'time (ms)', showgrid = FALSE, showline = TRUE, ticks = 'outside', zeroline = FALSE),
               yaxis = list(title = 'normalised F(x)', showgrid = FALSE, showline = TRUE, ticks = 'outside', zeroline = FALSE),
               legend = list(x = 0.98, y = 0.02, xanchor = 'right', yanchor = 'bottom'))
      
      if (!is.na(decay_time)) {
        p <- p %>%
          add_trace(x = c(0, t_max), y = rep(-du/min_val, 2),
                    showlegend = FALSE,
                    line = list(color = 'darkgrey', dash = 'dot')) %>%
          add_trace(x = c(0, t_max), y = rep(-dl/min_val, 2),
                    showlegend = FALSE,
                    line = list(color = 'darkgrey', dash = 'dot')) %>%
          add_annotations(x = t_max * 0.8, y = -0.3,
                          text = paste(up_pct, '-', lo_pct, 'decay =', decay_time, 'ms'),
                          showarrow = FALSE,
                          font = list(color = 'darkgrey', size = 11))
      } else {
        p <- p %>%
          add_annotations(x = t_max * 0.8, y = -0.3,
                          text = 'Decay not fully reached',
                          showarrow = FALSE,
                          font = list(color = 'indianred', size = 11))
      }
      
      p
    })
    
  }

  shinyApp(ui, server)

}

# analyseFIT - Shiny App for Fitting Any Custom Function

analyseFIT <- function() {
  
  ui <- fluidPage(
    create_dark_mode_css_fit(),
    add_busy_spinner(spin = "fading-circle", position = "top-right", 
                    color = "#3c8dbc", height = "60px", width = "60px"),
    titlePanel('Custom Function Fitting'),
    sidebarLayout(create_sidebar_panel_fit(), create_main_panel_fit())
  )
  
  server <- function(input, output, session) {
    state <- reactiveValues(
      response = NULL,
      analysis = NULL,
      accumulated_results = list(),
      custom_func = NULL,
      func_params = NULL
    )
    
    # Debounced inputs
    baseline_debounced <- debounce(reactive(input$baseline), 800)
    stimulation_time_debounced <- debounce(reactive(input$stimulation_time), 800)
    dt_debounced <- debounce(reactive(input$dt), 800)
    
    # Load data
    uploaded_data <- reactive({
      req(input$file)
      ext <- tools::file_ext(input$file$name)
      load_uploaded_data(input$file$datapath, ext)
    })
    
    output$column_selector <- renderUI({
      req(uploaded_data())
      selectInput('data_col', 'Select Column to Analyse', 
                  choices = names(uploaded_data()))
    })
    
    # Parse custom function
    observeEvent(input$parse_function, {
      tryCatch({
        # Parse the function string
        func_text <- input$function_text
        state$custom_func <- eval(parse(text = func_text))
        
        # Extract parameter names
        state$func_params <- generate_param_names(state$custom_func)
        
        showNotification(
          paste0("Function parsed! Parameters: ", 
                 paste(state$func_params, collapse = ", ")),
          type = "message", duration = 5
        )
        
        # Update bounds inputs
        output$bounds_ui <- renderUI({
          req(state$func_params)
          n_params <- length(state$func_params)
          
          tagList(
            h4("Parameter Bounds (optional)"),
            fluidRow(
              column(6, 
                     textInput("lower_bounds", "Lower Bounds (comma-separated)",
                               placeholder = paste(rep("-Inf", n_params), collapse = ", "))
              ),
              column(6,
                     textInput("upper_bounds", "Upper Bounds (comma-separated)",
                               placeholder = paste(rep("Inf", n_params), collapse = ", "))
              )
            )
          )
        })
        
      }, error = function(e) {
        showNotification(
          paste("Function parsing failed:", e$message),
          type = "error", duration = 10
        )
      })
    })
    
    # Run initial analysis
    observeEvent(input$run_initial, {
      req(uploaded_data(), input$data_col)
      clear_fit_state(state)
      
      data_col <- uploaded_data()[[input$data_col]]
      ds <- as.numeric(input$ds)
      state$response <- downsample_data(data_col, ds)
      
      showNotification("Data loaded successfully!", type = "message", duration = 3)
    })
    
    # Run main fitting
    observeEvent(input$run_main, {
      req(state$response)
      
      tryCatch({
        # AUTO-PARSE FUNCTION if not already parsed
        if (is.null(state$custom_func)) {
          func_text <- input$function_text
          state$custom_func <- eval(parse(text = func_text))
          state$func_params <- generate_param_names(state$custom_func)
          
          showNotification(
            paste0("Auto-parsed function. Parameters: ", 
                   paste(state$func_params, collapse = ", ")),
            type = "message", duration = 3
          )
        }
        
        dt <- as.numeric(input$dt) * as.numeric(input$ds)
        y <- state$response
        if (any(is.na(y))) y <- y[!is.na(y)]
        
        # Parse bounds
        lower <- parse_bounds(input$lower_bounds, length(state$func_params))
        upper <- parse_bounds(input$upper_bounds, length(state$func_params))
        
        # Run FITN
        state$analysis <- FITN(
          response = y,
          dt = dt,
          func = state$custom_func,
          N = as.numeric(input$N),
          IEI = as.numeric(input$IEI),
          method = input$method,
          weight_method = input$weight_method,
          stimulation_time = as.numeric(input$stimulation_time),
          baseline = as.numeric(input$baseline),
          lower = lower,
          upper = upper,
          filter = input$filter,
          fc = as.numeric(input$fc),
          return.output = TRUE,
          show.output = FALSE,
          show.plot = FALSE
        )
        
        showNotification("Fit complete!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(
          paste("Fitting failed:", e$message),
          type = "error", duration = 10
        )
      })
    })
    
    # Add to results
    observeEvent(input$add_result, {
      req(state$analysis)
      col_name <- req(input$data_col)
      
      column_result <- list(
        column = col_name,
        fits = state$analysis$fits,
        fits.se = state$analysis$fits.se,
        gof = state$analysis$gof,
        AIC = state$analysis$AIC,
        BIC = state$analysis$BIC,
        model_message = state$analysis$model.message,
        traces = state$analysis$traces,
        func_text = input$function_text,
        metadata = create_fit_metadata(input, col_name)
      )
      
      existing_cols <- sapply(state$accumulated_results, function(x) x$column)
      existing_idx <- which(existing_cols == col_name)
      
      if (length(existing_idx) > 0) {
        state$accumulated_results[[existing_idx[1]]] <- column_result
        showNotification(paste0("Updated: ", col_name), type = "message", duration = 3)
      } else {
        state$accumulated_results[[length(state$accumulated_results) + 1]] <- column_result
        showNotification(
          paste0("Added: ", col_name, " (Total: ", length(state$accumulated_results), ")"),
          type = "message", duration = 3
        )
      }
    })
    
    # Clear results
    observeEvent(input$clear_results, {
      state$accumulated_results <- list()
      showNotification("All accumulated results cleared", type = "warning", duration = 3)
    })
    
    observeEvent(input$clear_output, { state$analysis <- NULL })
    
    # Plot
    output$plot <- renderPlot({
      req(state$response)
      req(state$analysis$traces)
      
      traces <- state$analysis$traces
      plot(traces$x, traces$y, type = 'l', col = 'black', lwd = 2,
           xlab = input$xlab, ylab = input$ylab, main = "Fit Results")
      lines(traces$x, traces$yfit, col = 'indianred', lwd = 2)
      legend("topright", legend = c("Data", "Fit"), 
             col = c("black", "red"), lwd = 2, bty = "n")
    })
    
    # Console output
    output$console <- renderPrint({
      if (!is.null(state$analysis)) {
        cat("=== Fit Results ===\n\n")
        cat("Parameters:\n")
        print(data.frame(
          Parameter = state$func_params,
          Estimate = state$analysis$fits,
          Std.Error = state$analysis$fits.se
        ))
        cat("\n")
        cat("Goodness of Fit:\n")
        cat("  Residual Std Error:", round(state$analysis$gof, 4), "\n")
        cat("  AIC:", round(state$analysis$AIC, 2), "\n")
        cat("  BIC:", round(state$analysis$BIC, 2), "\n")
        cat("\nModel Message:", state$analysis$model.message, "\n")
      } else {
        cat('No analysis output performed')
      }
    })
    
    # Summary
    output$accumulated_summary <- renderPrint({
      if (length(state$accumulated_results) == 0) {
        cat("No accumulated results yet. Analyze columns and click 'Add to Results'.\n")
      } else {
        cat(paste("Analysed:", length(state$accumulated_results), "\n\n"))
        for (result in state$accumulated_results) {
          cat("Column:", result$column, "\n")
          cat("Fits:", paste(round(result$fits, 4), collapse = ", "), "\n")
          cat("GoF:", round(result$gof, 4), "\n\n")
        }
      }
    })
    
    # Download handlers
    output$download_rdata <- downloadHandler(
      filename = function() {
        paste0("fit_analysis_", Sys.Date(), ".RData")
      },
      content = function(file) {
        results <- list(
          accumulated_results = state$accumulated_results,
          metadata = create_fit_metadata(input)
        )
        save(results, file = file)
      }
    )
    
    output$download_xlsx <- downloadHandler(
      filename = function() {
        paste0("fit_analysis_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        req(length(state$accumulated_results) > 0)
        
        wb <- openxlsx::createWorkbook()
        
        # Summary sheet
        summary_df <- do.call(rbind, lapply(state$accumulated_results, function(r) {
          data.frame(
            Column = r$column,
            matrix(r$fits, nrow = 1),
            GoF = r$gof,
            AIC = r$AIC,
            BIC = r$BIC
          )
        }))
        
        openxlsx::addWorksheet(wb, "Summary")
        openxlsx::writeData(wb, "Summary", summary_df)
        
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
  }
  
  shinyApp(ui, server)
}

# UI Components

create_sidebar_panel_fit <- function() {
  sidebarPanel(
    fileInput('file', 'Upload csv or xlsx', accept = c('.csv', '.xlsx')),
    uiOutput('column_selector'),
    
    tabsetPanel(
      id = "tabs",
      
      # Main Options
      tabPanel("Main Options",
        br(),
        numericInput('dt', 'dt (ms):', 0.1, min = 0.001, step = 0.001),
        numericInput('stimulation_time', 'Stimulation Time:', 100, min = 0),
        numericInput('baseline', 'Baseline:', 50, min = 0),
        numericInput('N', 'N:', 1, min = 1, step = 1),
        numericInput('IEI', 'IEI (ms):', 50, min = 1),
        
        selectInput('method', 'Fit Method:', 
                    choices = c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'),
                    selected = 'LM'),
        
        selectInput('weight_method', 'Weight Method:',
                    choices = c('none', '~y_sqrt', '~y'),
                    selected = 'none'),
        
        numericInput('ds', 'Downsample Factor:', 1, min = 1, step = 1),
        checkboxInput('filter', 'Filter', FALSE),
        numericInput('fc', 'Filter Cutoff (Hz):', 1000, min = 1)
      ),
      
      # Function Input
      tabPanel("Function",
        br(),
        h4("Define Custom Function"),
        p("Enter your function as R code. Must accept (params, x, N, IEI) arguments."),
        p("Example: function(params, x, N=1, IEI=50) { params[1]*exp(-params[2]*x) + params[3] }"),
        
        textAreaInput('function_text', 'Function Definition:', 
                      value = 'function(params, x, N=1, IEI=50) {\n  params[1]*exp(-params[2]*x) + params[3]\n}',
                      rows = 8, width = '100%'),
        
        actionButton('parse_function', 'Parse Function', 
                     class = 'btn-primary', icon = icon('code')),
        br(), br(),
        uiOutput('bounds_ui')
      ),
      
      # Plot Settings
      tabPanel("Plot Settings",
        br(),
        textInput('xlab', 'X-axis Label:', 'time (ms)'),
        textInput('ylab', 'Y-axis Label:', 'Response (pA)'),
        numericInput('lwd', 'Line Width:', 1.2, min = 0.1, step = 0.1)
      )
    ),
    
    br(),
    actionButton('run_initial', 'Load Data', class = 'btn-info', 
                 icon = icon('upload')),
    actionButton('run_main', 'Run Fit', class = 'btn-info', 
                 icon = icon('play')),
    actionButton('add_result', 'Add to Results', class = 'btn-success',
                 icon = icon('plus')),
    actionButton('clear_results', 'Clear Results', class = 'btn-warning',
                 icon = icon('trash')),
    br(), br(),
    downloadButton('download_rdata', 'Download RData', class = 'btn-primary'),
    downloadButton('download_xlsx', 'Download XLSX', class = 'btn-primary'),
    actionButton('clear_output', 'Clear Output', class = 'btn-danger')
  )
}

create_main_panel_fit <- function() {
  mainPanel(
    plotOutput('plot', height = '400px'),
    h3('Fit Output'),
    verbatimTextOutput('console'),
    h3('Summary'),
    verbatimTextOutput('accumulated_summary')
  )
}

create_dark_mode_css_fit <- function() {
  tags$head(
    tags$style(HTML("
      body {
        background-color: #2c3e50;
        color: #ecf0f1;
      }
      .well {
        background-color: #34495e;
        border-color: #2c3e50;
      }
      .btn-info {
        background-color: #3498db;
        border-color: #2980b9;
      }
      .btn-success {
        background-color: #27ae60;
        border-color: #229954;
      }
      .btn-warning {
        background-color: #f39c12;
        border-color: #e67e22;
      }
      .btn-danger {
        background-color: #e74c3c;
        border-color: #c0392b;
      }
      .btn-primary {
        background-color: #3c8dbc;
        border-color: #367fa9;
      }
      .nav-tabs > li.active > a {
        background-color: #34495e !important;
        color: #ecf0f1 !important;
      }
      .tab-content {
        background-color: #34495e;
        padding: 15px;
        border: 1px solid #2c3e50;
        border-top: none;
      }
      pre {
        background-color: #1a252f;
        color: #ecf0f1;
        border: 1px solid #0d1117;
      }
      h3, h4 {
        color: #3498db;
      }
    "))
  )
}

# Helper Functions

clear_fit_state <- function(state) {
  state$analysis <- NULL
}

create_fit_metadata <- function(input, col_name = NULL) {
  list(
    column = col_name,
    dt = input$dt,
    stimulation_time = input$stimulation_time,
    baseline = input$baseline,
    N = input$N,
    IEI = input$IEI,
    method = input$method,
    weight_method = input$weight_method,
    downsample = input$ds,
    filter = input$filter,
    fc = input$fc,
    function_text = input$function_text,
    timestamp = Sys.time()
  )
}

parse_bounds <- function(bounds_text, n_params) {
  if (is.null(bounds_text) || bounds_text == "") return(NULL)
  
  tryCatch({
    vals <- as.numeric(strsplit(bounds_text, ",")[[1]])
    if (length(vals) != n_params) {
      stop(paste("Expected", n_params, "values, got", length(vals)))
    }
    vals
  }, error = function(e) {
    NULL
  })
}

# analysePSP files and functions

# ==============================================
# ABF and extract stimulus-triggered traces
# ==============================================

extract_stim_triggered_traces <- function(
    abf_filename,
    signal_channel = 1,
    trigger_channel = 4,
    threshold = 0.3,
    baseline = 10.0,
    poststimulation = 300.0,
    minstimulation_interval = 5.0
) {
  abf <- readABF::readABF(abf_filename)
  
  sample_rate <- 1 / abf$samplingIntervalInSec
  samples_per_ms <- sample_rate / 1000
  
  pre_samples <- as.integer(baseline * samples_per_ms)
  post_samples <- as.integer(poststimulation * samples_per_ms)
  min_interval_samples <- as.integer(minstimulation_interval * samples_per_ms)
  
  # Data is matrix: rows = samples, columns = channels
  data_matrix <- abf$data[[1]]
  signal_data <- data_matrix[, signal_channel]
  trigger_data <- data_matrix[, trigger_channel]
  
  cat(sprintf('trigger channel range: %.4f to %.4f\n', min(trigger_data), max(trigger_data)))
  
  # rising edges
  above_threshold <- trigger_data > threshold
  edges <- diff(as.integer(above_threshold))
  rising_edges <- which(edges == 1) + 1
# ==========================================================================  
  cat(sprintf('raw rising edges found: %d\n', length(rising_edges)))
  
  # Filter by minimum interval
  if (length(rising_edges) > 1) {
    keep <- c(TRUE, diff(rising_edges) >= min_interval_samples)
    stim_onsets <- rising_edges[keep]
  } else {
    stim_onsets <- rising_edges
  }
  
  cat(sprintf('after interval filter: %d stimulation events\n', length(stim_onsets)))
  
  # Extract traces
  traces <- list()
  valid_onsets <- c()
  
  for (onset in stim_onsets) {
    start_idx <- onset - pre_samples
    end_idx <- onset + post_samples - 1
    
    if (start_idx >= 1 && end_idx <= length(signal_data)) {
      traces[[length(traces) + 1]] <- signal_data[start_idx:end_idx]
      valid_onsets <- c(valid_onsets, onset)
    }
  }
  
  traces_matrix <- do.call(cbind, traces)
  stim_onset_times <- valid_onsets / sample_rate
  
  n_samples <- pre_samples + post_samples
  time <- seq(-baseline, poststimulation, length.out = n_samples)
  
  cat(sprintf('extracted %d valid traces\n', ncol(traces_matrix)))
  cat(sprintf('trace dimensions: %d samples x %d traces\n', nrow(traces_matrix), ncol(traces_matrix)))
  
  return(list(
    traces = traces_matrix,
    stim_onset_times = stim_onset_times,
    time = time,
    sample_rate = sample_rate
  ))
}


# ==============================================
# filtering functions
# ==============================================

hann_filter <- function(data, window_size) {
  window <- signal::hanning(window_size)
  window <- window / sum(window)
  signal::filtfilt(filt = window, a = 1, x = data)
}

butter_lowpass_filter <- function(data, cutoff, fs, order = 4) {
  nyq <- 0.5 * fs
  normal_cutoff <- cutoff / nyq
  bf <- signal::butter(order, normal_cutoff, type = 'low')
  signal::filtfilt(bf, data)
}

# ==============================================
# kinetics calculations
# ==============================================

crossings_fractional_idx <- function(trace, thr, i_from, i_to) {
  n <- length(trace)
  if (!is.numeric(trace) || n < 2L || !is.finite(thr)) return(numeric(0))

  i_from <- as.integer(max(1L, i_from))
  i_to   <- as.integer(min(n - 1L, i_to))
  if (!is.finite(i_from) || !is.finite(i_to) || i_to < i_from) return(numeric(0))

  idx <- i_from:i_to
  y1 <- trace[idx]
  y2 <- trace[idx + 1L]

  ok <- is.finite(y1) & is.finite(y2)
  if (!any(ok)) return(numeric(0))

  s1 <- y1 - thr
  s2 <- y2 - thr

  # crossing if sign changes or touches threshold
  hit <- ok & ((s1 == 0) | (s2 == 0) | ((s1 > 0) & (s2 < 0)) | ((s1 < 0) & (s2 > 0)))
  if (!any(hit)) return(numeric(0))

  ii <- idx[hit]
  x_cross <- vapply(ii, function(k) {
    a <- trace[k]
    b <- trace[k + 1L]
    if (!is.finite(a) || !is.finite(b)) return(NA_real_)
    if (abs(b - a) < .Machine$double.eps) return(as.numeric(k))
    k + (thr - a) / (b - a)
  }, numeric(1))

  x_cross[is.finite(x_cross)]
}

calc_rise_time <- function(trace, peak_idx, baseline_idx, limits, sample_rate) {
  if (length(limits) < 2 || any(!is.finite(limits[1:2])) || !is.finite(sample_rate) || sample_rate <= 0) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  lower_pct <- limits[1] / 100
  upper_pct <- limits[2] / 100

  n <- length(trace)
  if (!is.numeric(trace) || n < 2L) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  peak_idx <- as.integer(peak_idx)
  baseline_idx <- as.integer(baseline_idx)
  if (!is.finite(peak_idx) || !is.finite(baseline_idx)) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  peak_idx <- max(1L, min(n, peak_idx))
  baseline_idx <- max(1L, min(n, baseline_idx))
  if (peak_idx <= baseline_idx) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  peak_val <- trace[peak_idx]
  if (!is.finite(peak_val)) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  lower_val <- peak_val * lower_pct
  upper_val <- peak_val * upper_pct

  i_from <- baseline_idx
  i_to <- peak_idx - 1L

  lower_cross <- crossings_fractional_idx(trace, lower_val, i_from, i_to)
  upper_cross <- crossings_fractional_idx(trace, upper_val, i_from, i_to)

  if (length(lower_cross) == 0 || length(upper_cross) == 0) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  idx_lower <- as.numeric(stats::median(lower_cross))
  idx_upper <- as.numeric(stats::median(upper_cross))
  if (!is.finite(idx_lower) || !is.finite(idx_upper) || idx_upper < idx_lower) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  rise_time_ms <- (idx_upper - idx_lower) / sample_rate * 1000
  if (!is.finite(rise_time_ms) || rise_time_ms < 0) {
    return(list(rise_time_ms = NA, idx_lower = NULL, idx_upper = NULL))
  }

  list(rise_time_ms = rise_time_ms, idx_lower = idx_lower, idx_upper = idx_upper)
}

calc_decay_time <- function(trace, peak_idx, end_idx, limits, sample_rate) {
  if (length(limits) < 2 || any(!is.finite(limits[1:2])) || !is.finite(sample_rate) || sample_rate <= 0) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  upper_pct <- limits[1] / 100
  lower_pct <- limits[2] / 100

  n <- length(trace)
  if (!is.numeric(trace) || n < 2L) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  peak_idx <- as.integer(peak_idx)
  end_idx  <- as.integer(end_idx)
  if (!is.finite(peak_idx) || !is.finite(end_idx)) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  peak_idx <- max(1L, min(n, peak_idx))
  end_idx  <- max(1L, min(n, end_idx))
  if (end_idx <= peak_idx) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  peak_val <- trace[peak_idx]
  if (!is.finite(peak_val)) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  upper_val <- peak_val * upper_pct
  lower_val <- peak_val * lower_pct

  i_from <- peak_idx
  i_to <- end_idx - 1L

  upper_cross <- crossings_fractional_idx(trace, upper_val, i_from, i_to)
  lower_cross <- crossings_fractional_idx(trace, lower_val, i_from, i_to)

  if (length(upper_cross) == 0 || length(lower_cross) == 0) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  idx_upper <- as.numeric(stats::median(upper_cross))
  idx_lower <- as.numeric(stats::median(lower_cross))
  if (!is.finite(idx_upper) || !is.finite(idx_lower) || idx_lower < idx_upper) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  decay_time_ms <- (idx_lower - idx_upper) / sample_rate * 1000
  if (!is.finite(decay_time_ms) || decay_time_ms < 0) {
    return(list(decay_time_ms = NA, idx_upper = NULL, idx_lower = NULL))
  }

  list(decay_time_ms = decay_time_ms, idx_upper = idx_upper, idx_lower = idx_lower)
}

trapz <- function(x, y) {
  # Trapezoidal integration (base R replacement for pracma::trapz)
  n <- length(y)
  sum(diff(x) * (y[-n] + y[-1]) / 2)
}

calc_area <- function(trace, start_idx, end_idx, sample_rate) {
  dt_ms <- 1000 / sample_rate
  segment <- trace[start_idx:end_idx]
  x <- seq(0, by = dt_ms, length.out = length(segment))
  trapz(x, segment)
}

find_peak_with_smoothing <- function(trace, smooth_win = 1) {
  if (smooth_win > 1) {
    smoothed <- stats::filter(trace, rep(1/smooth_win, smooth_win), sides = 2)
    smoothed[is.na(smoothed)] <- trace[is.na(smoothed)]
  } else {
    smoothed <- trace
  }
  peak_idx <- which.max(smoothed)
  peak_val <- smoothed[peak_idx]
  return(list(peak_idx = peak_idx, peak_val = peak_val))
}

# ==============================================
# main analysis functions
# ==============================================

create_psp_result <- function(trace_idx, baseline_mean, baseline_sd, is_failure, peak_amplitude, peak_time_ms, 
  rise_time_ms, decay_time_ms, area_mVms, raw_subtracted, filtered_subtracted) {
  
  list(
    trace_idx = trace_idx,
    baseline_mean = baseline_mean,
    baseline_sd = baseline_sd,
    is_failure = is_failure,
    peak_amplitude = peak_amplitude,
    peak_time_ms = peak_time_ms,
    rise_time_ms = rise_time_ms,
    decay_time_ms = decay_time_ms,
    area_mVms = area_mVms,
    raw_subtracted = raw_subtracted,
    filtered_subtracted = filtered_subtracted
  )
}

PSPanalysis <- function(traces, time, sample_rate, baseline = 10.0, filter_type = 'butterworth',
  cutoff = 1000, hann_win = 21, butter_order = 4, failure_threshold = 3, latency_window_ms = c(1.0, 15.0),
  rise_limits = c(10, 90), decay_limits = c(90, 10), peak_smooth_win = 5, area_window_ms = c(0, 30), plot_trim_ms = 2
) {

  n_traces <- ncol(traces)

  stim_idx <- which.min(abs(time))
  pre_stim_samples <- as.integer(baseline * sample_rate / 1000)
  baseline_start_idx <- stim_idx - pre_stim_samples

  latency_start_idx <- stim_idx + as.integer(latency_window_ms[1] * sample_rate / 1000)
  latency_end_idx   <- stim_idx + as.integer(latency_window_ms[2] * sample_rate / 1000)

  area_start_idx <- stim_idx + as.integer(area_window_ms[1] * sample_rate / 1000)
  area_end_idx   <- stim_idx + as.integer(area_window_ms[2] * sample_rate / 1000)

  # clamp indices
  baseline_start_idx <- max(1, baseline_start_idx)
  latency_start_idx  <- max(1, latency_start_idx)
  latency_end_idx    <- min(length(time), latency_end_idx)
  area_start_idx     <- max(1, area_start_idx)
  area_end_idx       <- min(length(time), area_end_idx)

  filter_on <- !(is.null(filter_type) || filter_type == 'none')

  # trim filtered response when filter is ON
  trim_n <- if (filter_on) as.integer(plot_trim_ms * sample_rate / 1000) else 0L
  n_samp <- nrow(traces)
  valid <- rep(TRUE, n_samp)
  if (trim_n > 0 && (2 * trim_n) < n_samp) {
    valid[1:trim_n] <- FALSE
    valid[(n_samp - trim_n + 1):n_samp] <- FALSE
  }

  results <- vector('list', n_traces)

  for (i in 1:n_traces) {
    trace <- traces[, i]

    # filter (or not)
    if (filter_type == 'butterworth') {
      filtered <- butter_lowpass_filter(trace, cutoff, sample_rate, butter_order)
    } else if (filter_type == 'hann') {
      filtered <- hann_filter(trace, hann_win)
    } else {
      filtered <- trace
    }

    # baseline indices (time-based window) + optional trim mask
    bidx <- baseline_start_idx:(stim_idx - 1)
    bidx <- bidx[bidx >= 1 & bidx <= n_samp]
    if (filter_on) bidx <- bidx[valid[bidx]]  #  trimmed baseline when filtering

    if (length(bidx) < 2) {
      # avoid crash
      baseline_mean <- 0
      baseline_sd <- NA_real_
    } else {
      # calculate baseline from filtered if filter is ON, else from raw
      if (filter_on) {
        baseline_mean <- mean(filtered[bidx], na.rm = TRUE)
        baseline_sd   <- sd(filtered[bidx], na.rm = TRUE)
      } else {
        baseline_mean <- mean(trace[bidx], na.rm = TRUE)
        baseline_sd   <- sd(trace[bidx], na.rm = TRUE)
      }
    }

    raw_subtracted      <- trace - baseline_mean
    filtered_subtracted <- filtered - baseline_mean

    # peak search always uses filtered_subtracted (even if identical to raw when filter off)
    if (latency_end_idx < latency_start_idx) {
      results[[i]] <- create_psp_result(
        trace_idx = i - 1, baseline_mean = baseline_mean, baseline_sd = baseline_sd,
        is_failure = TRUE, peak_amplitude = 0, peak_time_ms = NA,
        rise_time_ms = NA, decay_time_ms = NA, area_mVms = 0,
        raw_subtracted = raw_subtracted, filtered_subtracted = filtered_subtracted
      )
      next
    }

    search_segment <- filtered_subtracted[latency_start_idx:latency_end_idx]
    if (length(search_segment) == 0) {
      results[[i]] <- create_psp_result(
        trace_idx = i - 1, baseline_mean = baseline_mean, baseline_sd = baseline_sd,
        is_failure = TRUE, peak_amplitude = 0, peak_time_ms = NA,
        rise_time_ms = NA, decay_time_ms = NA, area_mVms = 0,
        raw_subtracted = raw_subtracted, filtered_subtracted = filtered_subtracted
      )
      next
    }

    peak_result <- find_peak_with_smoothing(search_segment, peak_smooth_win)
    peak_idx    <- latency_start_idx + peak_result$peak_idx - 1
    peak_time_ms   <- time[peak_idx]
    peak_amplitude <- filtered_subtracted[peak_idx]

    threshold  <- failure_threshold * baseline_sd
    is_failure <- is.na(threshold) || (peak_amplitude < threshold)

    if (is_failure) {
      max_pos <- max(search_segment, na.rm = TRUE)
      min_neg <- min(search_segment, na.rm = TRUE)
      failure_amplitude <- ifelse(abs(max_pos) >= abs(min_neg), max_pos, min_neg)

      results[[i]] <- create_psp_result(
        trace_idx = i - 1, baseline_mean = baseline_mean, baseline_sd = baseline_sd,
        is_failure = TRUE, peak_amplitude = failure_amplitude, peak_time_ms = peak_time_ms,
        rise_time_ms = NA, decay_time_ms = NA, area_mVms = 0,
        raw_subtracted = raw_subtracted, filtered_subtracted = filtered_subtracted
      )
      next
    }

    rise_result  <- calc_rise_time(filtered_subtracted, peak_idx, stim_idx, rise_limits, sample_rate)
    decay_result <- calc_decay_time(filtered_subtracted, peak_idx, area_end_idx, decay_limits, sample_rate)
    area_mVms    <- calc_area(filtered_subtracted, area_start_idx, area_end_idx, sample_rate)

    results[[i]] <- create_psp_result(
      trace_idx = i - 1, baseline_mean = baseline_mean, baseline_sd = baseline_sd,
      is_failure = FALSE, peak_amplitude = peak_amplitude, peak_time_ms = peak_time_ms,
      rise_time_ms = rise_result$rise_time_ms, decay_time_ms = decay_result$decay_time_ms,
      area_mVms = area_mVms,
      raw_subtracted = raw_subtracted, filtered_subtracted = filtered_subtracted
    )
  }

  results
}

# ==============================================
# summary
# ==============================================

summarise_results <- function(results) {
  successes <- Filter(function(r) !r$is_failure, results)
  failures <- Filter(function(r) r$is_failure, results)
  
  n_total <- length(results)
  n_success <- length(successes)
  n_failure <- length(failures)
  
  summary <- list(
    n_total = n_total,
    n_success = n_success,
    n_failure = n_failure,
    failure_rate = ifelse(n_total > 0, n_failure / n_total, NA),
    success_rate = ifelse(n_total > 0, n_success / n_total, NA)
  )
  
  if (n_success > 0) {
    amps <- sapply(successes, function(r) r$peak_amplitude)
    rises <- sapply(successes, function(r) r$rise_time_ms)
    rises <- rises[!is.na(rises)]
    decays <- sapply(successes, function(r) r$decay_time_ms)
    decays <- decays[!is.na(decays)]
    areas <- sapply(successes, function(r) r$area_mVms)
    
    summary$amplitude_mean <- mean(amps)
    summary$amplitude_sd <- sd(amps)
    summary$rise_time_mean <- ifelse(length(rises) > 0, mean(rises), NA)
    summary$rise_time_sd <- ifelse(length(rises) > 0, sd(rises), NA)
    summary$decay_time_mean <- ifelse(length(decays) > 0, mean(decays), NA)
    summary$decay_time_sd <- ifelse(length(decays) > 0, sd(decays), NA)
    summary$area_mean <- mean(areas)
    summary$area_sd <- sd(areas)
  }
  
  return(summary)
}

print_summary <- function(summary, label = '') {
  cat('\n')
  cat(strrep('=', 60), '\n')
  cat(ifelse(label != '', sprintf('Summary: %s', label), 'Summary'), '\n')
  cat(strrep('=', 60), '\n')
  cat(sprintf('total traces:    %d\n', summary$n_total))
  cat(sprintf('successes:       %d (%.1f%%)\n', summary$n_success, summary$success_rate * 100))
  cat(sprintf('failures:        %d (%.1f%%)\n', summary$n_failure, summary$failure_rate * 100))
  
  if (summary$n_success > 0) {
    cat('\nSuccessful responses:\n')
    cat(sprintf('  amplitude:     %.3f +/- %.3f mV\n', summary$amplitude_mean, summary$amplitude_sd))
    cat(sprintf('  rise time:     %.3f +/- %.3f ms\n', summary$rise_time_mean, summary$rise_time_sd))
    cat(sprintf('  decay time:    %.3f +/- %.3f ms\n', summary$decay_time_mean, summary$decay_time_sd))
    cat(sprintf('  area:          %.3f +/- %.3f mVms\n', summary$area_mean, summary$area_sd))
  }
  cat(strrep('=', 60), '\n')
}

print_individual_results <- function(results, label = '') {
  successes <- Filter(function(r) !r$is_failure, results)
  
  if (label != '') {
    cat(sprintf('\n%s\n', label))
    cat(strrep('-', nchar(label)), '\n')
  }
  
  if (length(successes) == 0) {
    cat('no successful responses.\n')
    return(invisible(NULL))
  }
  
  cat(sprintf('%6s  %10s  %10s  %11s  %13s\n', 'trace', 'amp (mV)', 'rise (ms)', 'decay (ms)', 'area (mVms)'))
  cat(strrep('-', 60), '\n')
  
  for (r in successes) {
    cat(sprintf('%6d  %10.3f  %10.3f  %11.3f  %13.3f\n',
                r$trace_idx, r$peak_amplitude, r$rise_time_ms, r$decay_time_ms, r$area_mVms))
  }
  cat(strrep('-', 60), '\n')
}


# ==============================================
# plot
# ==============================================

plot_stim_traces <- function(
    traces,
    sample_rate,
    baseline = 10.0,
    poststimulation = 100.0,
    extraction_baseline = 20.0,
    plot_individual = TRUE,
    plot_average = FALSE,
    ylim = NULL,
    lwd = 1,
    avg_color = 'slateblue',
    stim_color = 'indianred'
) {
  
  # trim
  samples_to_trim <- as.integer((extraction_baseline - baseline) * sample_rate / 1000)
  if (samples_to_trim > 0) {
    traces <- traces[(samples_to_trim + 1):nrow(traces), , drop = FALSE]
  }
  
  n_traces <- ncol(traces)
  n_samples <- nrow(traces)
  
  dt <- 1000 / sample_rate
  time <- seq(-baseline, by = dt, length.out = n_samples)
  
  xlim <- c(-baseline, poststimulation)
  
  if (is.null(ylim)) {
    ylim <- c(floor(min(traces)), ceiling(max(traces)))
  }
  
  if (plot_individual && plot_average) {
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), tcl = -0.3, mgp = c(2.5, 0.5, 0))
  } else {
    par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), tcl = -0.3, mgp = c(2.5, 0.5, 0))
  }
  
  if (plot_individual) {
    plot(NULL, xlim = xlim, ylim = ylim, xlab = '', ylab = '',
         main = sprintf('individual traces (n=%d)', n_traces),
         bty = 'n', las = 1, xaxs = 'i', yaxs = 'i', axes = FALSE)
    
    # X-axis at bottom, from -baseline to poststimulation
    axis(1, at = pretty(xlim), pos = ylim[1])
    # Y-axis starting at -baseline
    axis(2, at = pretty(ylim), pos = -baseline, las = 1)
    mtext('time (ms)', side = 1, line = 2)
    mtext('mV', side = 2, line = 2.5)
    
    colors <- rainbow(n_traces, alpha = 0.4)
    for (i in 1:n_traces) {
      lines(time, traces[, i], col = colors[i], lwd = lwd)
    }
    abline(v = 0, col = stim_color, lty = 2, lwd = lwd)
  }
  
  if (plot_average) {
    mean_trace <- rowMeans(traces)
    se <- apply(traces, 1, sd) / sqrt(n_traces)
    ci_95 <- qt(0.975, df = n_traces - 1) * se
    
    plot(NULL, xlim = xlim, ylim = ylim, xlab = '', ylab = '',
         main = 'average trace +/- 95% CI',
         bty = 'n', las = 1, xaxs = 'i', yaxs = 'i', axes = FALSE)
    
    axis(1, at = pretty(xlim), pos = ylim[1])
    axis(2, at = pretty(ylim), pos = -baseline, las = 1)
    mtext('time (ms)', side = 1, line = 2)
    mtext('mV', side = 2, line = 2.5)
    
    polygon(c(time, rev(time)),
            c(mean_trace - ci_95, rev(mean_trace + ci_95)),
            col = adjustcolor(avg_color, alpha = 0.25), border = NA)
    
    lines(time, mean_trace, col = avg_color, lwd = 2 * lwd)
    abline(v = 0, col = stim_color, lty = 2, lwd = lwd)
  }
}


# ==============================================
# interactive browser (shiny app)
# ==============================================

browsePSP <- function(results, sample_rate, baseline = 10.0, y_unit = 'mV', latency_window_ms = c(1.0, 10.0), peak_smooth_win = 5,
    rise_limits = c(10, 90), decay_limits = c(90, 10), area_window_ms = c(0, 30), plot_trim_ms = c(2, 2), out_name = 'browse_results_output', 
    lwd = 2, dp = 4) {

  library(shiny)

  n_samples <- length(results[[1]]$filtered_subtracted)
  dt_ms <- 1000 / sample_rate
  time <- seq(from = -baseline, by = dt_ms, length.out = n_samples)
  stim_idx <- which.min(abs(time))
  n_traces <- length(results)

  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
        @media (prefers-color-scheme: dark) {
          body { background-color: #1e1e1e; color: #e0e0e0; }
          .well { background-color: #2d2d2d; border-color: #444; }
          .form-control { background-color: #2d2d2d; color: #c0c0c0 !important; border: 1px solid #555 !important; }
          input[type='number'], input[type='text'] { background-color: #2d2d2d !important; color: #c0c0c0 !important; }
          .selectize-input, .selectize-dropdown { background-color: #ffffff !important; color: #666666 !important; }
          .btn { background-color: #3d3d3d; color: #ffffff; border-color: #555; }
          .btn-primary, .action-button { background-color: #3c8dbc; color: #ffffff; }
          pre, code { background-color: #1a1a1a; color: #f0f0f0; }
        }
      "))
    ),
    titlePanel('PSP Results Browser'),
    sidebarLayout(
      sidebarPanel(
        width = 5,
        tabsetPanel(
          tabPanel('Navigation',
            br(),
            fluidRow(
              column(6, actionButton('prev_btn', '< Prev', class = 'btn-primary', style = 'width: 100%;')),
              column(6, actionButton('next_btn', 'Next >', class = 'btn-primary', style = 'width: 100%;'))
            ),
            hr(),
            actionButton('reclassify_btn', 'Toggle Success/Failure', class = 'btn-default', style = 'width: 100%;'),
            hr(),
            h4('Trace Info'),
            verbatimTextOutput('trace_info'),
            hr(),
            h4('Summary'),
            verbatimTextOutput('status_text')
          ),
          tabPanel('Settings',
            br(),
            numericInput('settings_sample_rate', 'Sample Rate (Hz):', value = sample_rate, min = 1),
            textInput('settings_y_unit', 'Y Unit:', value = y_unit),
            hr(),
            h5('Latency Window (ms)'),
            fluidRow(
              column(6, numericInput('settings_lat_start', 'Start:', value = latency_window_ms[1], step = 0.5)),
              column(6, numericInput('settings_lat_end', 'End:', value = latency_window_ms[2], step = 0.5))
            ),
            numericInput('settings_peak_smooth', 'Peak Smooth Window:', value = peak_smooth_win, min = 1),
            hr(),
            h5('Rise Time Limits (%)'),
            fluidRow(
              column(6, numericInput('settings_rise_lower', 'Lower:', value = rise_limits[1], min = 0, max = 100)),
              column(6, numericInput('settings_rise_upper', 'Upper:', value = rise_limits[2], min = 0, max = 100))
            ),
            h5('Decay Time Limits (%)'),
            fluidRow(
              column(6, numericInput('settings_decay_upper', 'Upper:', value = decay_limits[1], min = 0, max = 100)),
              column(6, numericInput('settings_decay_lower', 'Lower:', value = decay_limits[2], min = 0, max = 100))
            ),
            hr(),
            h5('Area Window (ms)'),
            fluidRow(
              column(6, numericInput('settings_area_start', 'Start:', value = area_window_ms[1], step = 1)),
              column(6, numericInput('settings_area_end', 'End:', value = area_window_ms[2], step = 1))
            ),
            h5('Plot Trim (ms)'),
            fluidRow(
              column(6, numericInput('settings_trim_start', 'Start:', value = plot_trim_ms[1], step = 0.5)),
              column(6, numericInput('settings_trim_end', 'End:', value = plot_trim_ms[2], step = 0.5))
            ),
            hr(),
            numericInput('settings_lwd', 'Line Width:', value = lwd, min = 0.5, max = 5, step = 0.5),
            textInput('settings_out_name', 'Output Name:', value = out_name)
          )
        )
      ),
      mainPanel(
        width = 7,
        plotOutput('trace_plot', height = '500px')
      )
    )
  )

  server <- function(input, output, session) {
    current_idx <- reactiveVal(1)
    results_reactive <- reactiveVal(results)
    local_out_name <- out_name

    plot_params <- reactive({
      t0 <- time[1] + input$settings_trim_start
      t1 <- time[length(time)] - input$settings_trim_end
      pidx <- which(time >= t0 & time <= t1)
      if (length(pidx) == 0) pidx <- seq_along(time)
      all_filt <- do.call(rbind, lapply(results, function(r) r$filtered_subtracted[pidx]))
      all_raw <- do.call(rbind, lapply(results, function(r) r$raw_subtracted[pidx]))
      ymin <- min(c(all_filt, all_raw), na.rm = TRUE)
      ymax <- max(c(all_filt, all_raw), na.rm = TRUE)
      yp <- (ymax - ymin) * 0.1
      if (!is.finite(yp) || yp == 0) yp <- 1e-6
      xl <- c(5 * floor(min(time[pidx]) / 5), 5 * ceiling(max(time[pidx]) / 5))
      list(plot_idx = pidx, xlim = xl, y_min = ymin, y_max = ymax, y_pad = yp)
    })

    analysis_idx <- reactive({
      sr <- input$settings_sample_rate
      list(
        latency_start = stim_idx + as.integer(input$settings_lat_start * sr / 1000),
        latency_end = stim_idx + as.integer(input$settings_lat_end * sr / 1000),
        area_start = stim_idx + as.integer(input$settings_area_start * sr / 1000),
        area_end = stim_idx + as.integer(input$settings_area_end * sr / 1000)
      )
    })

    observeEvent(input$prev_btn, {
      if (current_idx() > 1) current_idx(current_idx() - 1)
    })

    observeEvent(input$next_btn, {
      if (current_idx() < n_traces) current_idx(current_idx() + 1)
    })

    observeEvent(input$reclassify_btn, {
      idx <- current_idx()
      res <- results_reactive()
      r <- res[[idx]]
      aidx <- analysis_idx()
      search_segment <- r$filtered_subtracted[aidx$latency_start:aidx$latency_end]

      if (r$is_failure) {
        peak_result <- find_peak_with_smoothing(search_segment, input$settings_peak_smooth)
        peak_idx <- aidx$latency_start + peak_result$peak_idx - 1
        r$peak_amplitude <- r$filtered_subtracted[peak_idx]
        r$peak_time_ms <- time[peak_idx]
        r$is_failure <- FALSE
        rl <- c(input$settings_rise_lower, input$settings_rise_upper)
        dl <- c(input$settings_decay_upper, input$settings_decay_lower)
        r$rise_time_ms <- calc_rise_time(r$filtered_subtracted, peak_idx, stim_idx, rl, input$settings_sample_rate)$rise_time_ms
        r$decay_time_ms <- calc_decay_time(r$filtered_subtracted, peak_idx, aidx$area_end, dl, input$settings_sample_rate)$decay_time_ms
        r$area_mVms <- sum(r$filtered_subtracted[aidx$area_start:aidx$area_end]) * (1000 / input$settings_sample_rate)
      } else {
        max_pos <- max(search_segment, na.rm = TRUE)
        min_neg <- min(search_segment, na.rm = TRUE)
        r$peak_amplitude <- ifelse(abs(max_pos) >= abs(min_neg), max_pos, min_neg)
        r$is_failure <- TRUE
        r$rise_time_ms <- NA
        r$decay_time_ms <- NA
        r$area_mVms <- 0
      }

      res[[idx]] <- r
      results_reactive(res)
    })

    output$trace_plot <- renderPlot({
      idx <- current_idx()
      r <- results_reactive()[[idx]]
      pp <- plot_params()
      status <- if (r$is_failure) 'FAILURE' else 'SUCCESS'
      current_lwd <- input$settings_lwd

      par(mar = c(5, 6, 4, 2), mgp = c(3.5, 0.7, 0))

      plot(time[pp$plot_idx], r$raw_subtracted[pp$plot_idx],
           type = 'l',
           col = adjustcolor('black', 0.3),
           lwd = current_lwd * 0.6,
           xlim = pp$xlim,
           ylim = c(pp$y_min - pp$y_pad, pp$y_max + pp$y_pad),
           xlab = '', ylab = '',
           bty = 'n',
           xaxs = 'i', yaxs = 'i',
           axes = FALSE,
           main = sprintf('Trace %d / %d  [%s]', idx, n_traces, status),
           cex.main = 1.3)

      lines(time[pp$plot_idx], r$filtered_subtracted[pp$plot_idx],
            col = '#3c8dbc', lwd = current_lwd)

      axis(1, tcl = -0.3, lwd = current_lwd, cex.axis = 1.4)
      axis(2, las = 1, tcl = -0.3, lwd = current_lwd, cex.axis = 1.4)
      mtext('time (ms)', side = 1, line = 2.5, cex = 1.2)
      mtext(input$settings_y_unit, side = 2, line = 4, cex = 1.2)

      abline(h = 0, col = 'grey', lty = 3, lwd = current_lwd * 0.8)
      abline(v = 0, col = '#CD5C5C', lty = 2, lwd = current_lwd * 0.8)

      if (!r$is_failure && !is.na(r$peak_time_ms) &&
          r$peak_time_ms >= pp$xlim[1] && r$peak_time_ms <= pp$xlim[2]) {
        points(r$peak_time_ms, r$peak_amplitude,
               col = '#CD5C5C', pch = 16, cex = 1.5)
      }

      legend('topleft',
             legend = c('raw', 'filtered'),
             col = c(adjustcolor('black', 0.3), '#3c8dbc'),
             lwd = c(current_lwd * 0.6, current_lwd),
             bty = 'n', cex = 1.1)
    })

    output$trace_info <- renderPrint({
      idx <- current_idx()
      r <- results_reactive()[[idx]]
      cat('Trace:', idx, '/', n_traces, '\n')
      cat('Status:', if (r$is_failure) 'Failure' else 'Success', '\n')
      cat('Baseline SD:', round(r$baseline_sd, dp), input$settings_y_unit, '\n')
      cat('Amplitude:', round(r$peak_amplitude, dp), input$settings_y_unit, '\n')
      if (!r$is_failure) {
        if (!is.na(r$rise_time_ms)) cat('Rise time:', round(r$rise_time_ms, dp), 'ms\n')
        if (!is.na(r$decay_time_ms)) cat('Decay time:', round(r$decay_time_ms, dp), 'ms\n')
        cat('Area:', round(r$area_mVms, 2), 'mVms\n')
      }
    })

    output$status_text <- renderPrint({
      res <- results_reactive()
      n_success <- sum(vapply(res, function(r) !r$is_failure, logical(1)))
      n_failure <- n_traces - n_success
      cat('Total:', n_traces, '\n')
      cat('Successes:', n_success, sprintf('(%.1f%%)', 100 * n_success / n_traces), '\n')
      cat('Failures:', n_failure, sprintf('(%.1f%%)', 100 * n_failure / n_traces), '\n')
    })

    session$onSessionEnded(function() {
      assign(local_out_name, isolate(results_reactive()), envir = .GlobalEnv)
    })
  }

  shinyApp(ui, server)
}



# ==============================================
# # MODULAR VERSION OF analysePSP()
# ==============================================

# UI HELPER FUNCTIONS
create_dark_mode_css_psp <- function() {
  tags$head(
    tags$style(HTML("
      @media (prefers-color-scheme: dark) {
        body { background-color: #1e1e1e; color: #e0e0e0; }
        .well { background-color: #2d2d2d; border-color: #444; }
        .form-control { background-color: #2d2d2d; color: #c0c0c0 !important;
                        border: 1px solid #555 !important; }
        input[type='number'], input[type='text'] {
          background-color: #2d2d2d !important; color: #c0c0c0 !important; }
        .selectize-input, .selectize-dropdown {
          background-color: #ffffff !important; color: #666666 !important; }
        .btn { background-color: #3d3d3d; color: #ffffff; border-color: #555; }
        .btn-primary, .action-button { background-color: #3c8dbc; color: #ffffff; }
        pre, code { background-color: #1a1a1a; color: #f0f0f0; }
      }
      .sidebar-well { max-height: 85vh; overflow-y: auto; }
      pre { overflow-x: auto; font-size: 11px; max-height: 50vh; }
    "))
  )
}

create_data_tab_ui_psp <- function(traces, sample_rate, stimulation) {
  tabPanel('Data',
    br(),
    h4('Load Traces'),
    radioButtons('data_source', 'Source:',
      choices = c('From R (passed to function)' = 'r_object',
                  'From file (xlsx/csv)' = 'file'),
      selected = if (!is.null(traces)) 'r_object' else 'file'
    ),
    conditionalPanel(
      condition = "input.data_source == 'file'",
      fileInput('trace_file', 'Choose xlsx or csv:',
                accept = c('.xlsx', '.csv')),
      selectInput('file_sheet', 'Sheet (xlsx):', choices = NULL),
      checkboxInput('file_header', 'Header row', value = TRUE)
    ),
    hr(),
    h4('Trace Selection'),
    uiOutput('trace_range_ui'),
    hr(),
    h4('Sample Rate'),
    numericInput('data_sample_rate', 'Sample rate (Hz):',
      value = if (!is.null(sample_rate)) sample_rate else 20000, min = 1),
    # NOTE: stimulation = ms from trace start
    numericInput('data_stimulation', 'Stimulation time (ms from trace start):',
      value = stimulation, min = 0, step = 1),
    hr(),
    actionButton('load_data_btn', 'Load / Refresh Data',
      class = 'btn-primary', style = 'width: 100%;'),
    br(), br(),
    verbatimTextOutput('data_info_text'),
    br(),
    fileInput("load_rdata_file", "Load saved session (.RData):",
              accept = ".RData")
  )
}

create_analysis_tab_ui_psp <- function() {
  tabPanel('Analysis',
    br(),
    h4('PSPanalysis Settings'),
    # NOTE: default baseline = 20 ms
    numericInput('ana_baseline', 'Baseline (ms):', value = 20.0, min = 0, step = 1),
    # NOTE: filter is optional; default 'none'
    selectInput('ana_filter_type', 'Filter type:',
      choices = c('none', 'butterworth', 'hann'), selected = 'none'),
    conditionalPanel(
      condition = "input.ana_filter_type == 'butterworth'",
      numericInput('ana_cutoff', 'Cutoff (Hz):', value = 1000, min = 100, step = 100),
      numericInput('ana_butter_order', 'Butter order:', value = 4, min = 1, max = 10)
    ),
    conditionalPanel(
      condition = "input.ana_filter_type == 'hann'",
      numericInput('ana_hann_win', 'Hann window:', value = 21, min = 3, step = 2)
    ),
    numericInput('ana_failure_thresh', 'Failure threshold (x SD):',
      value = 5, min = 1, step = 0.5),
    hr(),
    h5('Latency Window (ms)'),
    fluidRow(
      column(6, numericInput('ana_lat_start', 'Start:', value = 1.0, step = 0.5)),
      column(6, numericInput('ana_lat_end', 'End:', value = 10.0, step = 0.5))
    ),
    h5('Rise Limits (%)'),
    fluidRow(
      column(6, numericInput('ana_rise_lo', 'Lower:', value = 10, min = 0, max = 100)),
      column(6, numericInput('ana_rise_hi', 'Upper:', value = 90, min = 0, max = 100))
    ),
    h5('Decay Limits (%)'),
    fluidRow(
      column(6, numericInput('ana_decay_hi', 'Upper:', value = 90, min = 0, max = 100)),
      column(6, numericInput('ana_decay_lo', 'Lower:', value = 10, min = 0, max = 100))
    ),
    numericInput('ana_peak_smooth', 'Peak smooth window:', value = 5, min = 1),
    h5('Area Window (ms)'),
    fluidRow(
      column(6, numericInput('ana_area_start', 'Start:', value = 0, step = 1)),
      column(6, numericInput('ana_area_end', 'End:', value = 30, step = 1))
    ),
    hr(),
    actionButton('run_analysis_btn', 'Run Analysis',
      class = 'btn-primary', style = 'width: 100%;'),
    br(), br(),
    textInput('export_prefix_ana', 'Export file prefix:', value = 'PSP'),
    actionButton('export_analysis_btn', 'Export Analysis (xlsx)',
      icon = icon('download'), class = 'btn-default', style = 'width: 100%;'),
    br(), br(),
    h4('Trace Info'),
    verbatimTextOutput('trace_info_text2'),
    hr(),
    h4('Summary'),
    verbatimTextOutput('browse_summary_text2')
  )
}

create_browse_tab_ui_psp <- function() {
  tabPanel('Browse',
    br(),
    fluidRow(
      column(6, actionButton('browse_prev_btn', '< Prev',
        class = 'btn-primary', style = 'width: 100%;')),
      column(6, actionButton('browse_next_btn', 'Next >',
        class = 'btn-primary', style = 'width: 100%;'))
    ),
    hr(),
    actionButton('reclassify_btn', 'Toggle Success/Failure',
      class = 'btn-default', style = 'width: 100%;'),
    hr(),
    h4('Trace Info'),
    verbatimTextOutput('trace_info_text'),
    hr(),
    h4('Summary'),
    verbatimTextOutput('browse_summary_text')
  )
}

create_fitting_tab_ui_psp <- function() {
  tabPanel('Fitting',
    br(),
    fluidRow(
      column(6, actionButton('fit_prev_btn', '< Prev',
        class = 'btn-primary', style = 'width: 100%;')),
      column(6, actionButton('fit_next_btn', 'Next >',
        class = 'btn-primary', style = 'width: 100%;'))
    ),
    hr(),
    h4('Product1 Fit'),
    numericInput('fit_n_iter', 'Iterations (n):', value = 100,
      min = 1, max = 500, step = 1),
    selectInput('fit_method', 'Method:',
      choices = c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'),
      selected = 'BF.LM'),
    numericInput('fit_limit_ms', 'Fit limit (ms from trace start):',
      value = 40, min = 5, max = 300, step = 1),
    # NOTE: label clarified, default 20 ms
    numericInput('fit_stim_time', 'Stimulation time (ms from trace start):',
      value = 20, min = 0, max = 200, step = 1),
    # NOTE: default baseline 20 ms
    numericInput('fit_baseline', 'Baseline (ms):',
      value = 20, min = 1, max = 100, step = 1),
    selectInput('fit_weight', 'Weights:',
      choices = c('none', '~y_sqrt', '~y'), selected = 'none'),
    hr(),
    checkboxInput('fit_use_filter', 'Filter before fitting', value = FALSE),
    conditionalPanel(
      condition = 'input.fit_use_filter == true',
      numericInput('fit_cutoff_hz', 'Cutoff (Hz):', value = 1000,
        min = 100, max = 10000, step = 100)
    ),
    numericInput('fit_downsample', 'Downsample:', value = 1,
      min = 1, max = 10, step = 1),
    hr(),
    h5('Rise / Decay Interval (fraction)'),
    fluidRow(
      column(6, numericInput('fit_interval_lo', 'Lower:', value = 0.1,
        min = 0.01, max = 0.4, step = 0.01)),
      column(6, numericInput('fit_interval_hi', 'Upper:', value = 0.9,
        min = 0.6, max = 0.99, step = 0.01))
    ),
    hr(),
    actionButton('fit_current_btn', 'Fit Current Trace',
      class = 'btn-primary', style = 'width: 100%;'),
    br(), br(),
    actionButton('fit_all_btn', 'Fit All Successes',
      class = 'btn-default', style = 'width: 100%;'),
    br(), br(),
    textInput('export_prefix_fit', 'Export file prefix:', value = 'PSP'),
    actionButton('export_fits_btn', 'Export Fits (xlsx)',
      icon = icon('download'), class = 'btn-default', style = 'width: 100%;'),
    br(), br(),
    actionButton('export_rdata_btn', 'Export RData',
      icon = icon('download'), class = 'btn-default', style = 'width: 100%;'),
    br(), br(),
    actionButton('export_svg_btn', 'Export SVG Plots',
      icon = icon('download'), class = 'btn-default', style = 'width: 100%;'),
    hr(),
    h4('Fit Results'),
    verbatimTextOutput('fit_results_text'),
    hr(),
    h4('All Fit Summary'),
    verbatimTextOutput('fit_summary_text')
  )
}

create_mle_settings_tab_ui_psp <- function() {
  tabPanel('MLE Settings',
    br(),
    numericInput('mle_iter', 'MLE Iterations:', value = 1000, min = 100, step = 100),
    numericInput('mle_metropolis_scale', 'Metropolis Scale:', value = 1.5,
      min = 0.1, max = 10, step = 0.1),
    numericInput('mle_fit_attempts', 'Fit Attempts:', value = 10, min = 1, max = 500, step = 5),
    checkboxInput('mle_rwm', 'Random Walk Metropolis', value = FALSE),
    selectInput('mle_method', 'MLE Optimisation Method:',
      choices = c('L-BFGS-B', 'Nelder-Mead', 'BFGS', 'CG', 'SANN', 'Brent'),
      selected = 'L-BFGS-B')
  )
}

create_plot_settings_tab_ui_psp <- function() {
  tabPanel('Plot Settings',
    br(),
    numericInput('lwd', 'Line Width:', value = 1.5, min = 0.5, max = 5, step = 0.1),
    numericInput('xbar', 'x-bar Length:', value = 10, min = 1, step = 1),
    numericInput('ybar', 'y-bar Length:', value = 0.5, min = 0.1, step = 0.1),
    textInput('xbar_lab', 'x-axis Units:', 'ms'),
    textInput('ybar_lab', 'y-axis Units:', 'mV'),
    textInput('xlim', 'x limits (e.g., 0,40):', '')
  )
}

create_average_tab_ui_psp <- function() {
  tabPanel('Average',
    br(),
    selectInput('avg_ci_method', 'CI Method:',
      choices = c('Normal' = 'normal', 'Bootstrap' = 'bootstrap'),
      selected = 'normal'),
    numericInput('avg_ci_level', 'CI Level:', value = 0.95,
      min = 0.5, max = 0.99, step = 0.01),
    conditionalPanel(
      condition = "input.avg_ci_method == 'bootstrap'",
      numericInput('avg_nboot', 'Bootstrap samples:', value = 999,
        min = 99, step = 100)
    ),
    numericInput('avg_alpha', 'Shade alpha:', value = 0.3,
      min = 0.05, max = 1, step = 0.05),
    hr(),
    actionButton('compute_avg_btn', 'Compute Average',
      class = 'btn-primary', style = 'width: 100%;'),
    hr(),
    actionButton('export_avg_svg_btn', 'Export Average SVG',
      icon = icon('download'), class = 'btn-default', style = 'width: 100%;'),
    hr(),
    verbatimTextOutput('avg_info_text')
  )
}

create_sidebar_panel_psp <- function(traces, sample_rate, stimulation) {
  sidebarPanel(
    width = 5,
    div(class = 'sidebar-well',
      tabsetPanel(
        id = 'main_tabs',
        create_data_tab_ui_psp(traces, sample_rate, stimulation),
        create_analysis_tab_ui_psp(),
        create_browse_tab_ui_psp(),
        create_fitting_tab_ui_psp(),
        create_mle_settings_tab_ui_psp(),
        create_plot_settings_tab_ui_psp(),
        create_average_tab_ui_psp()
      )
    )
  )
}

create_main_panel_psp <- function() {
  mainPanel(
    width = 7,
    conditionalPanel(
      condition = "input.main_tabs != 'Average'",
      plotOutput('main_plot', height = '700px')
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Average'",
      plotOutput('avg_plot', height = '700px')
    )
  )
}


# SUMMARY HELPER FUNCTION


render_psp_summary <- function(results, dp) {
  n_total <- length(results)
  successes <- Filter(function(r) !r$is_failure, results)
  n_success <- length(successes)
  n_failure <- n_total - n_success
  cat('Summary\n')
  cat(sprintf('Total traces:    %d\n', n_total))
  cat(sprintf('Successes:       %d (%.1f%%)\n', n_success, 100 * n_success / n_total))
  cat(sprintf('Failures:        %d (%.1f%%)\n', n_failure, 100 * n_failure / n_total))
  if (n_success > 0) {
    amps   <- sapply(successes, `[[`, 'peak_amplitude')
    rises  <- na.omit(sapply(successes, `[[`, 'rise_time_ms'))
    decays <- na.omit(sapply(successes, `[[`, 'decay_time_ms'))
    areas  <- sapply(successes, `[[`, 'area_mVms')
    cat('\nAverages:\n')
    cat(sprintf('  Amplitude:     %.3f +/- %.3f mV\n', mean(amps), sd(amps) / sqrt(length(amps))))
    if (length(rises) > 0)
      cat(sprintf('  Rise time:     %.3f +/- %.3f ms\n', mean(rises), sd(rises) / sqrt(length(rises))))
    if (length(decays) > 0)
      cat(sprintf('  Decay time:    %.3f +/- %.3f ms\n', mean(decays), sd(decays) / sqrt(length(decays))))
    cat(sprintf('  Area:          %.3f +/- %.3f mVms\n', mean(areas), sd(areas) / sqrt(length(areas))))
    w <- 60
    cat('\n', strrep('=', w), '\n')
    cat(sprintf('%6s  %10s  %10s  %11s  %13s\n',
      'trace', 'amp (mV)', 'rise (ms)', 'decay (ms)', 'area (mVms)'))
    cat(strrep('-', w), '\n')
    for (r in successes) {
      cat(sprintf('%6d  %10.3f  %10.3f  %11.3f  %13.3f\n',
        r$trace_idx + 1, r$peak_amplitude,
        r$rise_time_ms, r$decay_time_ms, r$area_mVms))
    }
    cat(strrep('-', w), '\n')
  }
}


# FIT HELPER FUNCTIONS

recompute_baseline_response <- function(r, fit_baseline_ms, fit_stim_ms, dt_ms_val) {
  n_pts <- length(r$raw_subtracted)
  if (n_pts < 1 || !is.finite(dt_ms_val) || dt_ms_val <= 0) return(r)

  bl_start <- max(1L, round((fit_stim_ms - fit_baseline_ms) / dt_ms_val) + 1L)
  bl_end <- max(bl_start, round(fit_stim_ms / dt_ms_val))
  bl_end <- min(bl_end, n_pts)

  orig_raw <- r$raw_subtracted + r$baseline_mean
  orig_filt <- r$filtered_subtracted + r$baseline_mean

  new_bl <- mean(orig_raw[bl_start:bl_end], na.rm = TRUE)
  if (!is.finite(new_bl)) return(r)

  r$raw_subtracted <- orig_raw - new_bl
  r$filtered_subtracted <- orig_filt - new_bl
  r$baseline_mean <- new_bl
  r
}

# pass original unsubtracted trace because FITN already handles the baseline
prepare_response_for_fit_psp <- function(r, fit_limit_ms, dt_ms_val) {
  response <- as.numeric(r$raw_subtracted + r$baseline_mean)
  x <- seq(0, (length(response) - 1) * dt_ms_val, by = dt_ms_val)
  response <- response[x < fit_limit_ms]
  return(response)
}

fit_single_trace_psp <- function(r, n_iter, method, weight_method, fit_limit_ms,
    downsample_factor, interval, use_filter, fc,
    stim_time, bl, MLEsettings, MLE.method, dt_ms_val) {
  response <- prepare_response_for_fit_psp(r, fit_limit_ms, dt_ms_val)
  fit_dt <- dt_ms_val
  # guard against NA downsample
  if (!is.na(downsample_factor) && downsample_factor > 1) {
    response <- response[seq(1, length(response), by = downsample_factor)]
    fit_dt <- fit_dt * downsample_factor
  }
  tryCatch({
    nFIT(
      response = response, n = n_iter, N = 1, IEI = 50,
      dt = fit_dt, func = product1N,
      method = method, weight_method = weight_method,
      stimulation_time = stim_time, baseline = bl,
      filter = use_filter, fc = fc, interval = interval,
      MLEsettings = MLEsettings, MLE.method = MLE.method,
      return.output = TRUE, show.output = FALSE, show.plot = FALSE
    )
  }, error = function(e) list(error = conditionMessage(e)))
}

# check fit args for NAs before calling fit
validate_fit_args_psp <- function(args) {
  for (nm in names(args)) {
    val <- args[[nm]]
    if (is.list(val)) {
      if (any(sapply(val, function(x) is.null(x) || any(is.na(x))))) return(FALSE)
    } else {
      if (is.null(val) || any(is.na(val))) return(FALSE)
    }
  }
  return(TRUE)
}


# SAVE HELPER FUNCTIONS


# show file-save modal dialog
save_with_prompt_psp <- function(default_name, save_id) {
  showModal(modalDialog(
    title = 'Save File',
    textInput(paste0(save_id, '_path'), 'Full file path:',
      value = file.path(path.expand('~/Downloads'), default_name)),
    footer = tagList(
      modalButton('Cancel'),
      actionButton(paste0(save_id, '_go'), 'Save', class = 'btn-primary')
    )
  ))
}

# save analysis results to xlsx
do_save_analysis_psp <- function(fpath, res, input, dp) {
  if (is.null(res)) return()
  rows <- lapply(seq_along(res), function(i) {
    r <- res[[i]]
    data.frame(
      Trace = r$trace_idx + 1,
      Status = if (r$is_failure) 'failure' else 'success',
      Amp_mV = round(r$peak_amplitude, dp),
      Rise_ms = if (!r$is_failure) round(r$rise_time_ms, dp) else NA,
      Decay_ms = if (!r$is_failure) round(r$decay_time_ms, dp) else NA,
      Area_mV.ms = if (!r$is_failure) round(r$area_mVms, dp) else NA,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  names(df) <- c('Trace', 'Status', 'Amp (mV)',
    sprintf('Rise %d-%d (ms)', input$ana_rise_lo, input$ana_rise_hi),
    sprintf('Decay %d-%d (ms)', input$ana_decay_hi, input$ana_decay_lo),
    'Area (mVms)')
  n_total <- length(res)
  successes <- Filter(function(r) !r$is_failure, res)
  n_success <- length(successes)
  n_failure <- n_total - n_success
  sum_rows <- list(
    data.frame(Metric = 'Total traces', Value = n_total, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Successes', Value = n_success, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Failures', Value = n_failure, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Success %', Value = round(100 * n_success / n_total, 1), SD = NA, stringsAsFactors = FALSE)
  )
  if (n_success > 0) {
    amps <- sapply(successes, `[[`, 'peak_amplitude')
    rises <- na.omit(sapply(successes, `[[`, 'rise_time_ms'))
    decays <- na.omit(sapply(successes, `[[`, 'decay_time_ms'))
    areas <- sapply(successes, `[[`, 'area_mVms')
    sum_rows <- c(sum_rows, list(
      data.frame(Metric = 'Amplitude (mV)', Value = round(mean(amps), dp), SD = round(sd(amps), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Rise time (ms)', Value = round(mean(rises), dp), SD = round(sd(rises), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Decay time (ms)', Value = round(mean(decays), dp), SD = round(sd(decays), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Area (mVms)', Value = round(mean(areas), dp), SD = round(sd(areas), dp), stringsAsFactors = FALSE)
    ))
  }
  summary_df <- do.call(rbind, sum_rows)
  meta <- data.frame(
    Setting = c('sample_rate', 'stimulation', 'baseline',
      'filter_type', 'cutoff', 'butter_order', 'hann_win',
      'failure_threshold', 'latency_start', 'latency_end',
      'rise_lo', 'rise_hi', 'decay_hi', 'decay_lo',
      'peak_smooth_win', 'area_start', 'area_end'),
    Value = c(input$data_sample_rate, input$data_stimulation, input$ana_baseline,
      input$ana_filter_type,
      if (input$ana_filter_type == 'butterworth') input$ana_cutoff else NA,
      if (input$ana_filter_type == 'butterworth') input$ana_butter_order else NA,
      if (input$ana_filter_type == 'hann') input$ana_hann_win else NA,
      input$ana_failure_thresh,
      input$ana_lat_start, input$ana_lat_end,
      input$ana_rise_lo, input$ana_rise_hi,
      input$ana_decay_hi, input$ana_decay_lo,
      input$ana_peak_smooth, input$ana_area_start, input$ana_area_end),
    stringsAsFactors = FALSE
  )
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 'analysis')
  openxlsx::writeData(wb, 'analysis', df)
  openxlsx::addWorksheet(wb, 'summary')
  openxlsx::writeData(wb, 'summary', summary_df)
  openxlsx::addWorksheet(wb, 'metadata')
  openxlsx::writeData(wb, 'metadata', meta)
  openxlsx::saveWorkbook(wb, fpath, overwrite = TRUE)
  removeModal()
  showNotification(paste('Saved to', fpath), type = 'message')
}

# save fit results to xlsx
do_save_fits_psp <- function(fpath, res, fs, input, dp) {
  if (is.null(res)) { removeModal(); return() }
  rows <- lapply(seq_along(res), function(i) {
    r <- res[[i]]; fo <- fs[[as.character(i)]]
    if (r$is_failure) {
      return(data.frame(Trace = r$trace_idx + 1, Status = 'failure',
        Amp = round(r$peak_amplitude, dp), Rise = NA, Decay = NA,
        Delay = NA, tpeak = NA, HW = NA, Area = NA, stringsAsFactors = FALSE))
    }
    if (!is.null(fo) && is.null(fo$error)) {
      o <- fo$output
      rc <- grep('^r\\d+_\\d+$', names(o), value = TRUE)[1]
      dc <- grep('^d\\d+_\\d+$', names(o), value = TRUE)[1]
      ac <- grep('^area', names(o), value = TRUE)[1]
      return(data.frame(Trace = r$trace_idx + 1, Status = 'success',
        Amp = round(o$A1, dp),
        Rise = if (!is.na(rc)) round(o[[rc]], dp) else NA,
        Decay = if (!is.na(dc)) round(o[[dc]], dp) else NA,
        Delay = round(o$delay, dp), tpeak = round(o$tpeak, dp),
        HW = if (!is.null(o$half_width)) round(o$half_width, dp) else NA,
        Area = if (!is.na(ac)) round(o[[ac]], dp) else NA, stringsAsFactors = FALSE))
    }
    data.frame(Trace = r$trace_idx + 1, Status = 'success',
      Amp = round(r$peak_amplitude, dp), Rise = round(r$rise_time_ms, dp),
      Decay = round(r$decay_time_ms, dp), Delay = NA, tpeak = NA, HW = NA,
      Area = round(r$area_mVms, dp), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  names(df) <- c('Trace', 'Status', 'Amp (mV)',
    sprintf('Rise %d-%d (ms)', round(input$fit_interval_lo * 100), round(input$fit_interval_hi * 100)),
    sprintf('Decay %d-%d (ms)', round(input$fit_interval_hi * 100), round(input$fit_interval_lo * 100)),
    'Delay (ms)', 'tpeak (ms)', 'Half width (ms)', 'Area (mVms)')
  # summary from fitted traces
  n_total <- length(res)
  fitted_rows <- list()
  for (key in names(fs)) {
    fo <- fs[[key]]
    if (!is.null(fo$error) || is.null(fo$output)) next
    o <- fo$output
    rc <- grep('^r\\d+_\\d+$', names(o), value = TRUE)[1]
    dc <- grep('^d\\d+_\\d+$', names(o), value = TRUE)[1]
    ac <- grep('^area', names(o), value = TRUE)[1]
    fitted_rows[[key]] <- list(amp = o$A1,
      rise = if (!is.na(rc)) o[[rc]] else NA,
      decay = if (!is.na(dc)) o[[dc]] else NA,
      area = if (!is.na(ac)) o[[ac]] else NA)
  }
  n_fitted <- length(fitted_rows)
  n_failure <- n_total - n_fitted
  sum_rows <- list(
    data.frame(Metric = 'Total traces', Value = n_total, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Successes', Value = n_fitted, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Failures', Value = n_failure, SD = NA, stringsAsFactors = FALSE),
    data.frame(Metric = 'Success %', Value = round(100 * n_fitted / n_total, 1), SD = NA, stringsAsFactors = FALSE)
  )
  if (length(fitted_rows) > 0) {
    amps <- sapply(fitted_rows, `[[`, 'amp')
    rises <- na.omit(sapply(fitted_rows, `[[`, 'rise'))
    decays <- na.omit(sapply(fitted_rows, `[[`, 'decay'))
    areas <- na.omit(sapply(fitted_rows, `[[`, 'area'))
    sum_rows <- c(sum_rows, list(
      data.frame(Metric = 'Amplitude (mV)', Value = round(mean(amps), dp), SD = round(sd(amps), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Rise time (ms)', Value = round(mean(rises), dp), SD = round(sd(rises), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Decay time (ms)', Value = round(mean(decays), dp), SD = round(sd(decays), dp), stringsAsFactors = FALSE),
      data.frame(Metric = 'Area (mVms)', Value = round(mean(areas), dp), SD = round(sd(areas), dp), stringsAsFactors = FALSE)
    ))
  }
  summary_df <- do.call(rbind, sum_rows)
  meta <- data.frame(
    Setting = c('sample_rate', 'stimulation', 'baseline',
      'filter_type', 'cutoff', 'butter_order', 'hann_win',
      'failure_threshold', 'latency_start', 'latency_end',
      'rise_lo', 'rise_hi', 'decay_hi', 'decay_lo',
      'peak_smooth_win', 'area_start', 'area_end',
      'fit_n_iter', 'fit_method', 'fit_limit_ms',
      'fit_stim_time', 'fit_baseline', 'fit_weight',
      'fit_use_filter', 'fit_cutoff_hz', 'fit_downsample',
      'fit_interval_lo', 'fit_interval_hi',
      'mle_iter', 'mle_metropolis_scale', 'mle_fit_attempts',
      'mle_rwm', 'mle_method'),
    Value = as.character(c(
      input$data_sample_rate, input$data_stimulation, input$ana_baseline,
      input$ana_filter_type,
      if (input$ana_filter_type == 'butterworth') input$ana_cutoff else NA,
      if (input$ana_filter_type == 'butterworth') input$ana_butter_order else NA,
      if (input$ana_filter_type == 'hann') input$ana_hann_win else NA,
      input$ana_failure_thresh,
      input$ana_lat_start, input$ana_lat_end,
      input$ana_rise_lo, input$ana_rise_hi,
      input$ana_decay_hi, input$ana_decay_lo,
      input$ana_peak_smooth, input$ana_area_start, input$ana_area_end,
      input$fit_n_iter, input$fit_method, input$fit_limit_ms,
      input$fit_stim_time, input$fit_baseline, input$fit_weight,
      input$fit_use_filter, input$fit_cutoff_hz, input$fit_downsample,
      input$fit_interval_lo, input$fit_interval_hi,
      input$mle_iter, input$mle_metropolis_scale, input$mle_fit_attempts,
      input$mle_rwm, input$mle_method)),
    stringsAsFactors = FALSE
  )
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 'fits')
  openxlsx::writeData(wb, 'fits', df)
  openxlsx::addWorksheet(wb, 'summary')
  openxlsx::writeData(wb, 'summary', summary_df)
  openxlsx::addWorksheet(wb, 'metadata')
  openxlsx::writeData(wb, 'metadata', meta)
  openxlsx::saveWorkbook(wb, fpath, overwrite = TRUE)
  removeModal()
  showNotification(paste('Saved to', fpath), type = 'message')
}

# save session as RData
do_save_rdata_psp <- function(fpath, analysis_res, fit_store_val, traces_val, sr_val, stim_val, input) {
  psp_data <- list(
    analysis = analysis_res, fits = fit_store_val,
    traces = traces_val, sample_rate = sr_val, stimulation = stim_val,
    plot_settings = list(
      lwd = input$lwd, xbar = input$xbar, ybar = input$ybar,
      xbar_lab = input$xbar_lab, ybar_lab = input$ybar_lab, xlim = input$xlim
    ),
    analysis_settings = list(
      baseline = input$ana_baseline,
      filter_type = input$ana_filter_type,
      cutoff = input$ana_cutoff,
      butter_order = input$ana_butter_order,
      hann_win = input$ana_hann_win,
      failure_thresh = input$ana_failure_thresh,
      lat_start = input$ana_lat_start, lat_end = input$ana_lat_end,
      rise_lo = input$ana_rise_lo, rise_hi = input$ana_rise_hi,
      decay_hi = input$ana_decay_hi, decay_lo = input$ana_decay_lo,
      peak_smooth = input$ana_peak_smooth,
      area_start = input$ana_area_start, area_end = input$ana_area_end
    )
  )
  save(psp_data, file = fpath)
  removeModal()
  showNotification(paste('Saved to', fpath), type = 'message')
}

# save SVG plots to zip
do_save_svg_psp <- function(fpath, res, fs, input, pta_full, stim) {
  if (is.null(res)) { removeModal(); return() }
  lwd <- as.numeric(input$lwd)
  xbar <- as.numeric(input$xbar); ybar <- as.numeric(input$ybar)
  xbar_lab <- input$xbar_lab; ybar_lab <- input$ybar_lab
  xlim_vals <- if (nchar(input$xlim) > 0) {
    v <- as.numeric(unlist(strsplit(input$xlim, ',')))
    if (length(v) == 2 && !any(is.na(v))) v else NULL
  } else NULL
  all_raw <- do.call(c, lapply(res, function(rr) rr$raw_subtracted))
  ymin <- min(all_raw, na.rm = TRUE); ymax <- max(all_raw, na.rm = TRUE)
  ypad <- (ymax - ymin) * 0.1
  if (!is.finite(ypad) || ypad == 0) ypad <- 1e-6
  ylim <- c(ymin - ypad, ymax + ypad)
  if (!is.null(xlim_vals)) { xlim <- xlim_vals
  } else { xlim <- c(5 * floor(min(pta_full) / 5), 5 * ceiling(max(pta_full) / 5)) }
  tmpdir <- tempfile(); dir.create(tmpdir); svg_files <- c()
  for (ii in seq_along(res)) {
    r <- res[[ii]]; fo <- fs[[as.character(ii)]]
    fname <- file.path(tmpdir, sprintf('trace_%d.svg', ii))
    svg(fname, width = 7, height = 5, bg = 'transparent')
    plot(pta_full, r$raw_subtracted, type = 'l', col = 'gray', lwd = lwd,
      xlim = xlim, ylim = ylim, axes = FALSE, xlab = '', ylab = '', bty = 'n')
    if (!is.null(fo) && is.null(fo$error) && !is.null(fo$traces)) {
      fp <- fo$fit_params
      st <- if (!is.null(fp$fit_stim_time)) fp$fit_stim_time else input$fit_stim_time
      bl <- if (!is.null(fp$fit_baseline)) fp$fit_baseline else input$fit_baseline      
      fit_x_orig <- fo$traces$x; fit_dx <- fit_x_orig[2] - fit_x_orig[1]
      ext_start <- xlim[1] - (st - bl); ext_end <- xlim[2] - (st - bl)
      if (ext_end > ext_start && fit_dx > 0) {
        extended_x <- seq(ext_start, ext_end, by = fit_dx)
        fits_adj <- fo$fits; fits_adj[4] <- fits_adj[4] + bl
        extended_y <- product1N(fits_adj, extended_x + fit_dx, N = 1, IEI = 50)
        extended_plot_x <- extended_x + (st - bl)
        in_r <- extended_plot_x >= xlim[1] & extended_plot_x <= xlim[2]
        lines(extended_plot_x[in_r], extended_y[in_r], col = 'indianred', lty = 3, lwd = 2 * lwd)
      }
    }
    usr <- par('usr')
    sb_y <- usr[3] + (usr[4] - usr[3]) / 20
    sb_x <- usr[2] - xbar - (usr[2] - usr[1]) * 0.05
    segments(sb_x, sb_y, sb_x + xbar, sb_y, lwd = lwd)
    segments(sb_x, sb_y, sb_x, sb_y + ybar, lwd = lwd)
    text((2 * sb_x + xbar) / 2, sb_y - (usr[4] - usr[3]) * 0.03,
      paste(xbar, xbar_lab), adj = c(0.5, 1), cex = 0.6)
    text(sb_x - (usr[2] - usr[1]) * 0.02, (2 * sb_y + ybar) / 2,
      paste(ybar, ybar_lab), srt = 90, adj = c(0.5, 0.5), cex = 0.6)
    dev.off()
    svg_files <- c(svg_files, fname)
  }
  if (length(svg_files) > 0) zip(fpath, svg_files, flags = '-j')
  removeModal()
  showNotification(paste('Saved to', fpath), type = 'message')
}

do_save_avg_svg_psp <- function(fpath, ad, pta, res, input, stim) {
  if (is.null(ad) || is.null(pta) || is.null(res)) { removeModal(); return() }

  lwd <- as.numeric(input$lwd)
  xbar <- as.numeric(input$xbar); ybar <- as.numeric(input$ybar)
  xbar_lab <- input$xbar_lab; ybar_lab <- input$ybar_lab
  shade_alpha <- if (!is.null(input$avg_alpha)) input$avg_alpha else 0.3

  xlim_vals <- if (nchar(input$xlim) > 0) {
    v <- as.numeric(unlist(strsplit(input$xlim, ",")))
    if (length(v) == 2 && !any(is.na(v))) v else NULL
  } else NULL
  xlim <- if (!is.null(xlim_vals)) xlim_vals else c(5 * floor(min(pta) / 5), 5 * ceiling(max(pta) / 5))

  all_raw <- do.call(c, lapply(res, function(rr) rr$raw_subtracted))
  ymin <- min(all_raw, na.rm = TRUE); ymax <- max(all_raw, na.rm = TRUE)
  ypad <- (ymax - ymin) * 0.1
  if (!is.finite(ypad) || ypad == 0) ypad <- 1e-6
  shared_ylim <- c(ymin - ypad, ymax + ypad)

  svg(fpath, width = 7, height = 8, bg = 'transparent')
  on.exit(dev.off(), add = TRUE)

  draw_avg_panel <- function(ci_data, title_text) {
    plot(pta, ci_data$avg, type = 'n', xlim = xlim, ylim = shared_ylim,
      xlab = '', ylab = '', bty = 'n', axes = FALSE, main = title_text, cex.main = 1.3)
    polygon(c(pta, rev(pta)), c(ci_data$lo, rev(ci_data$hi)),
      col = adjustcolor('steelblue', alpha.f = shade_alpha), border = NA)
    lines(pta, ci_data$avg, col = 'steelblue', lwd = lwd * 1.5)
    abline(h = 0, col = 'grey', lty = 3, lwd = lwd * 0.8)
    abline(v = stim, col = '#CD5C5C', lty = 2, lwd = lwd * 0.8)

    usr <- par('usr')
    sb_y <- usr[3] + (usr[4] - usr[3]) / 20
    sb_x <- usr[2] - xbar - (usr[2] - usr[1]) * 0.05
    segments(sb_x, sb_y, sb_x + xbar, sb_y, lwd = lwd)
    segments(sb_x, sb_y, sb_x, sb_y + ybar, lwd = lwd)
    text((2 * sb_x + xbar) / 2, sb_y - (usr[4] - usr[3]) * 0.03, paste(xbar, xbar_lab), adj = c(0.5, 1), cex = 0.8)
    text(sb_x - (usr[2] - usr[1]) * 0.02, (2 * sb_y + ybar) / 2, paste(ybar, ybar_lab), srt = 90, adj = c(0.5, 0.5), cex = 0.8)
  }

  par(mfrow = c(2, 1), mar = c(3, 3, 3, 2))
  if (!is.null(ad$success)) {
    draw_avg_panel(ad$success, sprintf('Successes (n=%d) - mean +/- %.0f%% CI', ad$success$n, input$avg_ci_level * 100))
  } else {
    plot.new(); text(0.5, 0.5, 'No successes to average.', cex = 1.2, col = 'grey50')
  }
  draw_avg_panel(ad$all, sprintf('All traces (n=%d) - mean +/- %.0f%% CI', ad$all$n, input$avg_ci_level * 100))

  removeModal()
  showNotification(paste('Saved to', fpath), type = 'message')
}

# AVERAGE HELPER FUNCTION

compute_ci_psp <- function(mat, ci_method, ci_level, nboot) {
  n <- ncol(mat)
  avg <- rowMeans(mat, na.rm = TRUE)
  if (n < 2) return(list(avg = avg, lo = avg, hi = avg, n = n))
  if (ci_method == 'normal') {
    se <- apply(mat, 1, sd, na.rm = TRUE) / sqrt(n)
    z <- qnorm(1 - (1 - ci_level) / 2)
    lo <- avg - z * se
    hi <- avg + z * se
  } else {
    boot_means <- matrix(NA, nrow = nrow(mat), ncol = nboot)
    for (b in seq_len(nboot)) {
      idx <- sample(n, n, replace = TRUE)
      boot_means[, b] <- rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    }
    alpha_q <- (1 - ci_level) / 2
    lo <- apply(boot_means, 1, quantile, probs = alpha_q, na.rm = TRUE)
    hi <- apply(boot_means, 1, quantile, probs = 1 - alpha_q, na.rm = TRUE)
  }
  list(avg = avg, lo = lo, hi = hi, n = n)
}

recompute_baseline_psp <- function(r, fit_baseline_ms, fit_stim_ms, dt_ms_val) {
  n_pts <- length(r$raw_subtracted)
  bl_start <- max(1, round((fit_stim_ms - fit_baseline_ms) / dt_ms_val) + 1)
  bl_end <- max(bl_start, round(fit_stim_ms / dt_ms_val))
  bl_end <- min(bl_end, n_pts)

  # reconstruct original traces, then subtract new baseline from both
  orig_raw <- r$raw_subtracted + r$baseline_mean
  orig_filt <- r$filtered_subtracted + r$baseline_mean
  new_bl <- mean(orig_raw[bl_start:bl_end])

  r$raw_subtracted <- orig_raw - new_bl
  r$filtered_subtracted <- orig_filt - new_bl
  r$baseline_mean <- new_bl
  r
}

# MAIN FUNCTION - analysePSP (Modular Version)

analysePSP <- function(traces = NULL, sample_rate = NULL, stimulation = 20.0, plot_trim_ms = 2, dp = 4,
    color_palette = 'Vik', alpha = 1, colorbar_reverse = TRUE, out_name = 'analysePSP_output') {
  library(shiny)

  # UI
  ui <- fluidPage(
    create_dark_mode_css_psp(),
    titlePanel('PSP Analyser'),
    sidebarLayout(
      create_sidebar_panel_psp(traces, sample_rate, stimulation),
      create_main_panel_psp()
    )
  )

  # SERVER
  server <- function(input, output, session) {
    # reactive state
    loaded_traces    <- reactiveVal(NULL)   # raw trace matrix (samples x traces)
    loaded_sr        <- reactiveVal(NULL)   # sample rate
    loaded_stim      <- reactiveVal(NULL)   # stim ms from extraction (ms from trace start)
    loaded_time      <- reactiveVal(NULL)
    analysis_results <- reactiveVal(NULL)   # list of results from PSPanalysis
    current_idx      <- reactiveVal(1)
    fit_store        <- reactiveVal(list())
    avg_data         <- reactiveVal(NULL)   # average trace CI data

    # sync export prefix between tabs
    observeEvent(input$export_prefix_ana, {
      if (!is.null(input$export_prefix_fit) && input$export_prefix_fit != input$export_prefix_ana)
        updateTextInput(session, 'export_prefix_fit', value = input$export_prefix_ana)
    }, ignoreInit = TRUE)
    observeEvent(input$export_prefix_fit, {
      if (!is.null(input$export_prefix_ana) && input$export_prefix_ana != input$export_prefix_fit)
        updateTextInput(session, 'export_prefix_ana', value = input$export_prefix_fit)
    }, ignoreInit = TRUE)

    # settings
    n_traces <- reactive({
      res <- analysis_results()
      if (is.null(res)) return(0)
      length(res)
    })
    dt_ms <- reactive({
      sr <- loaded_sr()
      if (is.null(sr)) return(0.05)
      1000 / sr
    })
    # determine if analysis filtering is ON
    filter_on <- reactive({
      ft <- input$ana_filter_type
      !is.null(ft) && ft != 'none'
    })
    # internal time axis (0 = stim) and plot axis (0 = trace start)
    int_time <- reactive({
      res <- analysis_results()
      if (is.null(res)) return(NULL)
      n <- length(res[[1]]$filtered_subtracted)
      stim <- loaded_stim()
      seq(from = -stim, by = dt_ms(), length.out = n)
    })
    stim_idx <- reactive({
      ta <- int_time()
      if (is.null(ta)) return(1)
      which.min(abs(ta))
    })
    plot_time <- reactive({
      ta <- int_time()
      if (is.null(ta)) return(NULL)
      ta + loaded_stim()
    })

    # collect fit parameters from inputs
    get_fit_args <- function() {
      list(
        n_iter    = input$fit_n_iter,
        method    = input$fit_method,
        weight_method = input$fit_weight,
        fit_limit_ms  = input$fit_limit_ms,
        downsample_factor = input$fit_downsample,
        interval   = c(input$fit_interval_lo, input$fit_interval_hi),
        use_filter = input$fit_use_filter,
        fc         = if (isTRUE(input$fit_use_filter)) input$fit_cutoff_hz else 1000,
        stim_time  = input$fit_stim_time,
        bl         = input$fit_baseline,
        MLEsettings = list(
          iter             = as.numeric(input$mle_iter),
          metropolis.scale = as.numeric(input$mle_metropolis_scale),
          fit.attempts     = as.numeric(input$mle_fit_attempts),
          RWm              = input$mle_rwm
        ),
        MLE.method = input$mle_method
      )
    }

    # store current fit parameters for later restoration
    build_fit_params <- function(args) {
      list(
        fit_n_iter = args$n_iter, fit_method = args$method,
        fit_limit_ms = args$fit_limit_ms, fit_stim_time = args$stim_time,
        fit_baseline = args$bl, fit_weight = args$weight_method,
        fit_use_filter = args$use_filter, fit_cutoff_hz = args$fc,
        fit_downsample = args$downsample_factor,
        fit_interval_lo = args$interval[1], fit_interval_hi = args$interval[2],
        mle_iter = args$MLEsettings$iter,
        mle_metropolis_scale = args$MLEsettings$metropolis.scale,
        mle_fit_attempts = args$MLEsettings$fit.attempts,
        mle_rwm = args$MLEsettings$RWm, mle_method = args$MLE.method
      )
    }

    # restore fit params into UI inputs
    restore_fit_params <- function(p) {
      updateNumericInput(session, 'fit_n_iter', value = p$fit_n_iter)
      updateSelectInput(session, 'fit_method', selected = p$fit_method)
      updateNumericInput(session, 'fit_limit_ms', value = p$fit_limit_ms)
      updateNumericInput(session, 'fit_stim_time', value = p$fit_stim_time)
      updateNumericInput(session, 'fit_baseline', value = p$fit_baseline)
      updateSelectInput(session, 'fit_weight', selected = p$fit_weight)
      updateCheckboxInput(session, 'fit_use_filter', value = p$fit_use_filter)
      updateNumericInput(session, 'fit_cutoff_hz', value = p$fit_cutoff_hz)
      updateNumericInput(session, 'fit_downsample', value = p$fit_downsample)
      updateNumericInput(session, 'fit_interval_lo', value = p$fit_interval_lo)
      updateNumericInput(session, 'fit_interval_hi', value = p$fit_interval_hi)
      updateNumericInput(session, 'mle_iter', value = p$mle_iter)
      updateNumericInput(session, 'mle_metropolis_scale', value = p$mle_metropolis_scale)
      updateNumericInput(session, 'mle_fit_attempts', value = p$mle_fit_attempts)
      updateCheckboxInput(session, 'mle_rwm', value = p$mle_rwm)
      updateSelectInput(session, 'mle_method', selected = p$mle_method)
    }

    # =========================================================================
    # TAB 1 load data
    # =========================================================================

    # if traces passed directly from R
    observe({
      if (!is.null(traces) && !is.null(sample_rate)) {
        loaded_traces(traces)
        loaded_sr(sample_rate)
        loaded_stim(stimulation)
        n <- nrow(traces)
        dt <- 1000 / sample_rate
        ta <- seq(from = -stimulation, by = dt, length.out = n)
        loaded_time(ta)
      }
    })

    # update sheet selector when xlsx uploaded
    observeEvent(input$trace_file, {
      req(input$trace_file)
      ext <- tools::file_ext(input$trace_file$name)
      if (ext == 'xlsx') {
        sheets <- openxlsx::getSheetNames(input$trace_file$datapath)
        updateSelectInput(session, 'file_sheet', choices = sheets, selected = sheets[1])
      }
    })

    # trace range UI (dynamic based on loaded data)
    output$trace_range_ui <- renderUI({
      tr <- loaded_traces()
      if (is.null(tr)) return(helpText('No traces loaded yet.'))
      nc <- ncol(tr)
      tagList(
        numericInput('trace_col_start', 'First trace column:', value = 1, min = 1, max = nc),
        numericInput('trace_col_end', 'Last trace column:', value = nc, min = 1, max = nc)
      )
    })

    observeEvent(input$load_data_btn, {
      if (input$data_source == 'file') {
        req(input$trace_file)
        ext <- tools::file_ext(input$trace_file$name)
        if (ext == 'csv') {
          tr <- as.matrix(read.csv(input$trace_file$datapath, header = input$file_header))
        } else if (ext == 'xlsx') {
          tr <- as.matrix(openxlsx::read.xlsx(input$trace_file$datapath,
            sheet = input$file_sheet, colNames = input$file_header))
        } else {
          showNotification('Unsupported file type.', type = 'error')
          return()
        }
        loaded_traces(tr)
      }
      loaded_sr(input$data_sample_rate)
      loaded_stim(input$data_stimulation)
      # build time axis for the loaded traces (0 = stim)
      tr <- loaded_traces()
      if (!is.null(tr)) {
        n <- nrow(tr)
        stim <- input$data_stimulation
        sr <- input$data_sample_rate
        dt <- 1000 / sr
        ta <- seq(from = -stim, by = dt, length.out = n)
        loaded_time(ta)
      }
      # auto-populate export prefix
      if (input$data_source == 'file' && !is.null(input$trace_file)) {
        prefix <- tools::file_path_sans_ext(basename(input$trace_file$name))
      } else {
        prefix <- 'PSP'
      }
      updateTextInput(session, 'export_prefix_ana', value = prefix)
      updateTextInput(session, 'export_prefix_fit', value = prefix)
      # reset downstream
      analysis_results(NULL)
      current_idx(1)
      fit_store(list())
      avg_data(NULL)
      showNotification(sprintf('Loaded %d samples x %d traces', nrow(tr), ncol(tr)),
        type = 'message')
    })

    # load saved RData session
    observeEvent(input$load_rdata_file, {
      req(input$load_rdata_file)
      env <- new.env()
      load(input$load_rdata_file$datapath, envir = env)
      prefix <- tools::file_path_sans_ext(basename(input$load_rdata_file$name))
      updateTextInput(session, 'export_prefix_ana', value = prefix)
      updateTextInput(session, 'export_prefix_fit', value = prefix)
      psp_data <- env$psp_data

      loaded_traces(psp_data$traces)
      loaded_sr(psp_data$sample_rate)
      loaded_stim(psp_data$stimulation)

      # rebuild time axis
      if (!is.null(psp_data$traces)) {
        n <- nrow(psp_data$traces)
        dt <- 1000 / psp_data$sample_rate
        ta <- seq(from = -psp_data$stimulation, by = dt, length.out = n)
        loaded_time(ta)
      }

      analysis_results(psp_data$analysis)
      fit_store(psp_data$fits)
      current_idx(1)
      avg_data(NULL)

      # update UI
      updateNumericInput(session, 'data_sample_rate', value = psp_data$sample_rate)
      updateNumericInput(session, 'data_stimulation', value = psp_data$stimulation)

      # restore plot settings
      if (!is.null(psp_data$plot_settings)) {
        ps <- psp_data$plot_settings
        updateNumericInput(session, 'lwd', value = ps$lwd)
        updateNumericInput(session, 'xbar', value = ps$xbar)
        updateNumericInput(session, 'ybar', value = ps$ybar)
        updateTextInput(session, 'xbar_lab', value = ps$xbar_lab)
        updateTextInput(session, 'ybar_lab', value = ps$ybar_lab)
        updateTextInput(session, 'xlim', value = ps$xlim)
      }

      # restore analysis settings
      if (!is.null(psp_data$analysis_settings)) {
        as <- psp_data$analysis_settings
        updateNumericInput(session, 'ana_baseline', value = as$baseline)
        updateSelectInput(session, 'ana_filter_type', selected = as$filter_type)
        updateNumericInput(session, 'ana_cutoff', value = as$cutoff)
        updateNumericInput(session, 'ana_butter_order', value = as$butter_order)
        updateNumericInput(session, 'ana_hann_win', value = as$hann_win)
        updateNumericInput(session, 'ana_failure_thresh', value = as$failure_thresh)
        updateNumericInput(session, 'ana_lat_start', value = as$lat_start)
        updateNumericInput(session, 'ana_lat_end', value = as$lat_end)
        updateNumericInput(session, 'ana_rise_lo', value = as$rise_lo)
        updateNumericInput(session, 'ana_rise_hi', value = as$rise_hi)
        updateNumericInput(session, 'ana_decay_hi', value = as$decay_hi)
        updateNumericInput(session, 'ana_decay_lo', value = as$decay_lo)
        updateNumericInput(session, 'ana_peak_smooth', value = as$peak_smooth)
        updateNumericInput(session, 'ana_area_start', value = as$area_start)
        updateNumericInput(session, 'ana_area_end', value = as$area_end)
      }

      # restore first fitted trace's params if available
      fs <- psp_data$fits
      if (length(fs) > 0) {
        first_key <- names(fs)[1]
        fo <- fs[[first_key]]
        if (!is.null(fo$fit_params)) restore_fit_params(fo$fit_params)
      }

      n_fits <- length(fs)
      n_ana <- if (!is.null(psp_data$analysis)) length(psp_data$analysis) else 0
      showNotification(sprintf('Session restored: %d traces, %d analysed, %d fitted',
        ncol(psp_data$traces), n_ana, n_fits), type = 'message')
    })

    output$data_info_text <- renderPrint({
      tr <- loaded_traces()
      if (is.null(tr)) { cat('No data loaded.\n'); return() }
      cat('Traces:', ncol(tr), '\n')
      cat('Samples:', nrow(tr), '\n')
      cat('Sample rate:', loaded_sr(), 'Hz\n')
      cat('dt:', dt_ms(), 'ms\n')
      cat('Stimulation:', loaded_stim(), 'ms\n')
      dur <- nrow(tr) * dt_ms()
      cat('Duration:', round(dur, 1), 'ms\n')
    })

    # =========================================================================
    # TAB 2 analysis
    # =========================================================================

    observeEvent(input$run_analysis_btn, {
      tr <- loaded_traces()
      ta <- loaded_time()
      sr <- loaded_sr()
      if (is.null(tr) || is.null(ta) || is.null(sr)) {
        showNotification('Load data first.', type = 'warning')
        return()
      }
      c1 <- input$trace_col_start
      c2 <- input$trace_col_end
      if (is.null(c1) || is.null(c2)) { c1 <- 1; c2 <- ncol(tr) }
      c1 <- max(1, c1); c2 <- min(ncol(tr), c2)
      trace_subset <- tr[, c1:c2, drop = FALSE]
      withProgress(message = 'Running PSPanalysis...', value = 0.5, {
        results <- PSPanalysis(
          traces            = trace_subset,
          time              = ta,
          sample_rate       = sr,
          baseline          = input$ana_baseline,
          filter_type       = input$ana_filter_type,
          cutoff            = if (input$ana_filter_type == 'butterworth') input$ana_cutoff else 1000,
          butter_order      = if (input$ana_filter_type == 'butterworth') input$ana_butter_order else 4,
          hann_win          = if (input$ana_filter_type == 'hann') input$ana_hann_win else 21,
          failure_threshold = input$ana_failure_thresh,
          latency_window_ms = c(input$ana_lat_start, input$ana_lat_end),
          rise_limits       = c(input$ana_rise_lo, input$ana_rise_hi),
          decay_limits      = c(input$ana_decay_hi, input$ana_decay_lo),
          peak_smooth_win   = input$ana_peak_smooth,
          area_window_ms    = c(input$ana_area_start, input$ana_area_end),
          plot_trim_ms      = plot_trim_ms
        )
        setProgress(1)
      })
      analysis_results(results)
      current_idx(1)
      fit_store(list())
      avg_data(NULL)
      updateNumericInput(session, 'fit_stim_time', value = loaded_stim())
      updateNumericInput(session, 'fit_baseline', value = input$ana_baseline)
      showNotification(sprintf('Analysis complete: %d traces', length(results)), type = 'message')
    })

    # trace info renderer (shared between Analysis and Browse tabs)
    render_trace_info <- function() {
      res <- analysis_results()
      if (is.null(res)) { cat('Run analysis first.\n'); return() }
      idx <- current_idx()
      r <- res[[idx]]
      cat('Trace:', idx, '/', n_traces(), '\n')
      cat('Status:', if (r$is_failure) 'Failure' else 'Success', '\n')
      cat('Baseline SD:', round(r$baseline_sd, dp), 'mV\n')
      cat('Amplitude:', round(r$peak_amplitude, dp), 'mV\n')
      if (!r$is_failure) {
        if (!is.na(r$rise_time_ms))  cat('Rise time:', round(r$rise_time_ms, dp), 'ms\n')
        if (!is.na(r$decay_time_ms)) cat('Decay time:', round(r$decay_time_ms, dp), 'ms\n')
        cat('Area:', round(r$area_mVms, dp), 'mVms\n')
      }
    }
    output$trace_info_text2 <- renderPrint({ render_trace_info() })
    output$trace_info_text  <- renderPrint({ render_trace_info() })
    output$browse_summary_text2 <- renderPrint({
      res <- analysis_results()
      if (is.null(res)) { cat('Run analysis first.\n'); return() }
      render_psp_summary(res, dp)
    })

    # =========================================================================
    # Navigation (all tabs)
    # =========================================================================

    nav_prev <- function() { if (current_idx() > 1) current_idx(current_idx() - 1) }
    nav_next <- function() {
      nt <- n_traces()
      if (nt > 0 && current_idx() < nt) current_idx(current_idx() + 1)
    }
    observeEvent(input$browse_prev_btn, nav_prev())
    observeEvent(input$browse_next_btn, nav_next())
    observeEvent(input$fit_prev_btn, nav_prev())
    observeEvent(input$fit_next_btn, nav_next())
    observeEvent(current_idx(), {
      idx <- current_idx()
      fs <- fit_store()
      fo <- fs[[as.character(idx)]]
      if (!is.null(fo) && !is.null(fo$fit_params)) restore_fit_params(fo$fit_params)
    })

    # =========================================================================
    # TAB 3 browse
    # =========================================================================

    observeEvent(input$reclassify_btn, {
      res <- analysis_results()
      if (is.null(res)) return()
      idx <- current_idx()
      r <- res[[idx]]
      sr <- loaded_sr()
      ta <- int_time()
      si <- stim_idx()
      lat_start  <- si + as.integer(input$ana_lat_start * sr / 1000)
      lat_end    <- si + as.integer(input$ana_lat_end   * sr / 1000)
      area_start <- si + as.integer(input$ana_area_start * sr / 1000)
      area_end   <- si + as.integer(input$ana_area_end   * sr / 1000)
      lat_start  <- max(1, lat_start)
      lat_end    <- min(length(r$filtered_subtracted), lat_end)
      area_start <- max(1, area_start)
      area_end   <- min(length(r$filtered_subtracted), area_end)
      search_segment <- r$filtered_subtracted[lat_start:lat_end]
      if (r$is_failure) {
        peak_result <- find_peak_with_smoothing(search_segment, input$ana_peak_smooth)
        peak_idx <- lat_start + peak_result$peak_idx - 1
        r$peak_amplitude <- r$filtered_subtracted[peak_idx]
        r$peak_time_ms <- ta[peak_idx]
        r$is_failure <- FALSE
        rl <- c(input$ana_rise_lo, input$ana_rise_hi)
        dl <- c(input$ana_decay_hi, input$ana_decay_lo)
        r$rise_time_ms  <- calc_rise_time(r$filtered_subtracted, peak_idx, si, rl, sr)$rise_time_ms
        r$decay_time_ms <- calc_decay_time(r$filtered_subtracted, peak_idx, area_end, dl, sr)$decay_time_ms
        r$area_mVms     <- sum(r$filtered_subtracted[area_start:area_end]) * (1000 / sr)
      } else {
        max_pos <- max(search_segment, na.rm = TRUE)
        min_neg <- min(search_segment, na.rm = TRUE)
        r$peak_amplitude <- ifelse(abs(max_pos) >= abs(min_neg), max_pos, min_neg)
        r$is_failure <- TRUE
        r$rise_time_ms <- NA; r$decay_time_ms <- NA; r$area_mVms <- 0
      }
      res[[idx]] <- r
      # remove fit if toggled to failure
      if (r$is_failure) {
        fs <- fit_store()
        fs[[as.character(idx)]] <- NULL
        fit_store(fs)
      }
      analysis_results(res)
    })

    output$browse_summary_text <- renderPrint({
      res <- analysis_results()
      if (is.null(res)) return()
      render_psp_summary(res, dp)
    })

    # =========================================================================
    # TAB 4 fitting
    # =========================================================================
    
    update_baseline_after_fit <- function(idx, args) {
      res <- analysis_results()
      if (is.null(res)) return()
      r <- res[[idx]]
      fit_bl <- args$bl
      fit_st <- args$stim_time
      if (is.null(fit_bl) || is.null(fit_st) || is.na(fit_bl) || is.na(fit_st)) return()
      r <- recompute_baseline_psp(r, fit_bl, fit_st, dt_ms())
      res[[idx]] <- r
      analysis_results(res)
    }

    observeEvent(input$fit_current_btn, {
      res <- analysis_results()
      if (is.null(res)) { showNotification('Run analysis first.', type = 'warning'); return() }
      idx <- current_idx()
      r <- res[[idx]]
      if (r$is_failure) { showNotification('Cannot fit a failure.', type = 'warning'); return() }
      args <- get_fit_args()
      # guard against NA fit parameters (user mid-edit)
      if (!validate_fit_args_psp(args)) {
        showNotification('Fit parameters incomplete - check inputs.', type = 'warning')
        return()
      }
      withProgress(message = sprintf('Fitting trace %d ...', idx), value = 0.5, {
        fit_out <- fit_single_trace_psp(r,
          args$n_iter, args$method, args$weight_method, args$fit_limit_ms,
          args$downsample_factor, args$interval, args$use_filter, args$fc,
          args$stim_time, args$bl, args$MLEsettings, args$MLE.method, dt_ms())
        setProgress(1)
      })
      if (!is.null(fit_out$error)) {
        fs <- fit_store()
        fs[[as.character(idx)]] <- NULL
        fit_store(fs)
        showNotification(paste('Fit failed:', fit_out$error), type = 'error'); return()
      }
      fs <- fit_store()
      fs[[as.character(idx)]] <- fit_out
      fs[[as.character(idx)]]$fit_params <- build_fit_params(args)
      fit_store(fs)
      update_baseline_after_fit(idx, args)
      showNotification(sprintf('Trace %d fitted (GoF=%.4f)', idx, fit_out$gof), type = 'message')
    })

    observeEvent(input$fit_all_btn, {
      res <- analysis_results()
      if (is.null(res)) { showNotification('Run analysis first.', type = 'warning'); return() }
      success_idx <- which(sapply(res, function(r) !r$is_failure))
      if (length(success_idx) == 0) {
        showNotification('No successes to fit.', type = 'warning'); return()
      }
      args <- get_fit_args()
      if (!validate_fit_args_psp(args)) {
        showNotification('Fit parameters incomplete - check inputs.', type = 'warning')
        return()
      }
      fs <- fit_store()
      withProgress(message = 'Fitting all successes...', value = 0, {
        for (i in seq_along(success_idx)) {
          ii <- success_idx[i]
          setProgress(i / length(success_idx),
            detail = sprintf('Trace %d (%d/%d)', ii, i, length(success_idx)))
          fit_out <- fit_single_trace_psp(res[[ii]],
            args$n_iter, args$method, args$weight_method, args$fit_limit_ms,
            args$downsample_factor, args$interval, args$use_filter, args$fc,
            args$stim_time, args$bl, args$MLEsettings, args$MLE.method, dt_ms())
          if (is.null(fit_out$error)) {
            fs[[as.character(ii)]] <- fit_out
            fs[[as.character(ii)]]$fit_params <- build_fit_params(args)
          }
        }
      })
      fit_store(fs)
      # update baselines for all fitted traces
      res <- analysis_results()
      for (i in seq_along(success_idx)) {
        ii <- success_idx[i]
        if (!is.null(fs[[as.character(ii)]])) {
          fit_bl <- args$bl
          fit_st <- args$stim_time
          if (!is.null(fit_bl) && !is.null(fit_st) && !is.na(fit_bl) && !is.na(fit_st)) {
            r <- recompute_baseline_psp(res[[ii]], fit_bl, fit_st, dt_ms())
            res[[ii]] <- r
          }
        }
      }

      analysis_results(res)
      showNotification(sprintf('Fitted %d / %d', length(fs), length(success_idx)), type = 'message')
    })

    output$fit_results_text <- renderPrint({
      res <- analysis_results(); if (is.null(res)) { cat('Run analysis first.\n'); return() }
      idx <- current_idx()
      fs <- fit_store()
      fit_out <- fs[[as.character(idx)]]
      if (is.null(fit_out)) { cat('No fit for trace', idx, '\n'); return() }
      if (!is.null(fit_out$error)) { cat('Fit error:', fit_out$error, '\n'); return() }
      fo <- fit_out$output
      cat(sprintf('Trace %d - Product1 Fit\n', idx))
      cat(strrep('-', 40), '\n')
      cat('A1:        ', fo$A1, 'mV\n')
      cat('tau_rise:  ', fo[['τrise']], 'ms\n')
      cat('tau_decay: ', fo[['τdecay']], 'ms\n')
      cat('tpeak:     ', fo$tpeak, 'ms\n')
      rc <- grep('^r\\d+_\\d+$', names(fo), value = TRUE)
      if (length(rc) > 0) cat(sprintf('%-11s %s ms\n', paste0(rc[1], ':'), fo[[rc[1]]]))
      dc <- grep('^d\\d+_\\d+$', names(fo), value = TRUE)
      if (length(dc) > 0) cat(sprintf('%-11s %s ms\n', paste0(dc[1], ':'), fo[[dc[1]]]))
      cat('delay:     ', fo$delay, 'ms\n')
      ac <- grep('^area', names(fo), value = TRUE)
      if (length(ac) > 0) cat('area:      ', fo[[ac[1]]], 'mVms\n')
      if (!is.null(fo$half_width)) cat('half_width:', fo$half_width, 'ms\n')
      cat(strrep('-', 40), '\n')
      cat('gof:', round(fit_out$gof, 6), '\n')
    })

    output$fit_summary_text <- renderPrint({
      res <- analysis_results()
      fs <- fit_store()
      if (is.null(res)) { cat('Run analysis first.\n'); return() }
      if (length(fs) == 0) { cat('No fits yet.\n'); return() }
      n_total <- length(res)
      rows <- list()
      for (key in names(fs)) {
        fo <- fs[[key]]
        if (!is.null(fo$error)) next
        o  <- fo$output
        rc <- grep('^r\\d+_\\d+$', names(o), value = TRUE)[1]
        dc <- grep('^d\\d+_\\d+$', names(o), value = TRUE)[1]
        ac <- grep('^area', names(o), value = TRUE)[1]
        rows[[key]] <- list(
          trace = as.integer(key), amp = o$A1,
          rise  = if (!is.na(rc)) o[[rc]] else NA,
          decay = if (!is.na(dc)) o[[dc]] else NA,
          delay = o$delay, tpeak = o$tpeak,
          hw    = if ('half_width' %in% names(o)) o$half_width else NA,
          area  = if (!is.na(ac)) o[[ac]] else NA
        )
      }
      if (length(rows) == 0) { cat('No successful fits.\n'); return() }
      n_fitted <- length(rows)
      n_failure <- n_total - n_fitted
      cat('Summary\n')
      cat(sprintf('Total traces:    %d\n', n_total))
      cat(sprintf('Successes:       %d (%.1f%%)\n', n_fitted, 100 * n_fitted / n_total))
      cat(sprintf('Failures:        %d (%.1f%%)\n', n_failure, 100 * n_failure / n_total))
      amps   <- sapply(rows, `[[`, 'amp')
      rises  <- na.omit(sapply(rows, `[[`, 'rise'))
      decays <- na.omit(sapply(rows, `[[`, 'decay'))
      areas  <- na.omit(sapply(rows, `[[`, 'area'))
      cat('\nAverages:\n')
      cat(sprintf('  Amplitude:     %.3f +/- %.3f mV\n', mean(amps), sd(amps) / sqrt(length(amps))))
      if (length(rises) > 0)
        cat(sprintf('  Rise time:     %.3f +/- %.3f ms\n', mean(rises), sd(rises) / sqrt(length(rises))))
      if (length(decays) > 0)
        cat(sprintf('  Decay time:    %.3f +/- %.3f ms\n', mean(decays), sd(decays) / sqrt(length(decays))))
      if (length(areas) > 0)
        cat(sprintf('  Area:          %.3f +/- %.3f mVms\n', mean(areas), sd(areas) / sqrt(length(areas))))
      w <- 98
      cat('\n', strrep('=', w), '\n')
      cat(sprintf('%6s  %9s  %9s  %10s  %10s  %10s  %11s  %12s\n',
        'trace', 'amp (mV)', 'rise (ms)', 'decay (ms)',
        'delay (ms)', 'tpeak (ms)', 'half width', 'area (mVms)'))
      cat(strrep('-', w), '\n')
      ordered <- rows[order(sapply(rows, `[[`, 'trace'))]
      for (r in ordered) {
        cat(sprintf('%6d  %9.3f  %9.3f  %10.3f  %10.3f  %10.3f  %11.3f  %12.3f\n',
          r$trace, r$amp, r$rise, r$decay, r$delay, r$tpeak, r$hw, r$area))
      }
      cat(strrep('-', w), '\n')
    })
    # =========================================================================
    # Exports
    # =========================================================================

    # export Analysis
    observeEvent(input$export_analysis_btn, {
      save_with_prompt_psp(paste0(input$export_prefix_ana, '_analysis.xlsx'), 'save_ana')
    })
    observeEvent(input$save_ana_go, {
      fpath <- input$save_ana_path
      if (file.exists(fpath)) {
        showModal(modalDialog(title = 'File exists', paste0('Overwrite ', basename(fpath), '?'),
          footer = tagList(modalButton('Cancel'),
            actionButton('save_ana_overwrite', 'Overwrite', class = 'btn-primary'))))
      } else { do_save_analysis_psp(fpath, analysis_results(), input, dp) }
    })
    observeEvent(input$save_ana_overwrite, {
      do_save_analysis_psp(input$save_ana_path, analysis_results(), input, dp)
    })

    # export Fits
    observeEvent(input$export_fits_btn, {
      save_with_prompt_psp(paste0(input$export_prefix_fit, '_fits.xlsx'), 'save_fits')
    })
    observeEvent(input$save_fits_go, {
      fpath <- input$save_fits_path
      if (file.exists(fpath)) {
        showModal(modalDialog(title = 'File exists', paste0('Overwrite ', basename(fpath), '?'),
          footer = tagList(modalButton('Cancel'),
            actionButton('save_fits_overwrite', 'Overwrite', class = 'btn-primary'))))
      } else { do_save_fits_psp(fpath, analysis_results(), fit_store(), input, dp) }
    })
    observeEvent(input$save_fits_overwrite, {
      do_save_fits_psp(input$save_fits_path, analysis_results(), fit_store(), input, dp)
    })

    # export RData
    observeEvent(input$export_rdata_btn, {
      save_with_prompt_psp(paste0(input$export_prefix_fit, '.RData'), 'save_rdata')
    })
    observeEvent(input$save_rdata_go, {
      fpath <- input$save_rdata_path
      if (file.exists(fpath)) {
        showModal(modalDialog(title = 'File exists', paste0('Overwrite ', basename(fpath), '?'),
          footer = tagList(modalButton('Cancel'),
            actionButton('save_rdata_overwrite', 'Overwrite', class = 'btn-primary'))))
      } else {
        do_save_rdata_psp(fpath, analysis_results(), fit_store(),
          loaded_traces(), loaded_sr(), loaded_stim(), input)
      }
    })
    observeEvent(input$save_rdata_overwrite, {
      do_save_rdata_psp(input$save_rdata_path, analysis_results(), fit_store(),
        loaded_traces(), loaded_sr(), loaded_stim(), input)
    })

    # export SVG
    observeEvent(input$export_svg_btn, {
      save_with_prompt_psp(paste0(input$export_prefix_fit, '_fit_plots.zip'), 'save_svg')
    })
    observeEvent(input$save_svg_go, {
      fpath <- input$save_svg_path
      if (file.exists(fpath)) {
        showModal(modalDialog(title = 'File exists', paste0('Overwrite ', basename(fpath), '?'),
          footer = tagList(modalButton('Cancel'),
            actionButton('save_svg_overwrite', 'Overwrite', class = 'btn-primary'))))
      } else {
        do_save_svg_psp(fpath, analysis_results(), fit_store(), input, plot_time(), loaded_stim())
      }
    })
    observeEvent(input$save_svg_overwrite, {
      do_save_svg_psp(input$save_svg_path, analysis_results(), fit_store(), input, plot_time(), loaded_stim())
    })

    # export Average SVG
    observeEvent(input$export_avg_svg_btn, {
      save_with_prompt_psp(paste0(input$export_prefix_fit, '_average.svg'), 'save_avg_svg')
    })
    observeEvent(input$save_avg_svg_go, {
      fpath <- input$save_avg_svg_path
      if (file.exists(fpath)) {
        showModal(modalDialog(title = 'File exists', paste0('Overwrite ', basename(fpath), '?'),
          footer = tagList(modalButton('Cancel'),
            actionButton('save_avg_svg_overwrite', 'Overwrite', class = 'btn-primary'))))
      } else {
        do_save_avg_svg_psp(fpath, avg_data(), plot_time(), analysis_results(), input, loaded_stim())
      }
    })
    observeEvent(input$save_avg_svg_overwrite, {
      do_save_avg_svg_psp(input$save_avg_svg_path, avg_data(), plot_time(), analysis_results(), input, loaded_stim())
    })

    # =========================================================================
    # main plot
    # =========================================================================

    output$main_plot <- renderPlot({
      res <- analysis_results()
      lwd <- as.numeric(input$lwd)
      xbar <- as.numeric(input$xbar)
      ybar <- as.numeric(input$ybar)
      xbar_lab <- input$xbar_lab
      ybar_lab <- input$ybar_lab
      xlim_vals <- if (nchar(input$xlim) > 0) {
        v <- as.numeric(unlist(strsplit(input$xlim, ',')))
        if (length(v) == 2 && !any(is.na(v))) v else NULL
      } else NULL

      # no analysis so raw traces with AXES
      if (is.null(res)) {
        tr <- loaded_traces()
        if (is.null(tr)) {
          plot.new()
          text(0.5, 0.5, 'Load data in the Data tab, then run Analysis.',
            cex = 1.5, col = 'grey50')
          return()
        }
        n <- nrow(tr); sr <- loaded_sr()
        if (is.null(sr)) sr <- 20000
        dt <- 1000 / sr
        x <- seq(0, by = dt, length.out = n)
        current_lwd <- 2
        xlim <- c(20 * floor(min(x) / 20), 20 * ceiling(max(x) / 20))
        yr <- range(tr, na.rm = TRUE)
        ylim <- c(floor(yr[1]), ceiling(yr[2]))
        par(mar = c(5, 6, 4, 2), mgp = c(3.5, 0.7, 0))
        plot(NULL, xlim = xlim, ylim = ylim,
          xlab = '', ylab = '',
          main = sprintf('raw loaded traces (n=%d)', ncol(tr)),
          bty = 'n', xaxs = 'i', yaxs = 'i', axes = FALSE, cex.main = 1.3)
        axis(1, tcl = -0.3, lwd = current_lwd, cex.axis = 1.4)
        axis(2, las = 1, tcl = -0.3, lwd = current_lwd, cex.axis = 1.4)
        mtext('time (ms)', side = 1, line = 2.5, cex = 1.2)
        mtext('mV', side = 2, line = 4, cex = 1.2)
        base_cols <- .palette_cols(color_palette, ncol(tr), alpha, colorbar_reverse)
        cols <- adjustcolor(base_cols, alpha.f = 0.6)
        for (i in 1:ncol(tr)) lines(x, tr[, i], col = cols[i], lwd = current_lwd * 0.6)
        abline(h = 0, col = 'grey', lty = 3, lwd = current_lwd * 0.8)
        return()
      }

      # individual trace view (browse + fit) with scale bars
      idx <- current_idx()
      r <- res[[idx]]
      pta_full <- plot_time()
      stim <- loaded_stim()

      # set y range
      all_raw <- do.call(c, lapply(res, function(rr) rr$raw_subtracted))
      if (isTRUE(filter_on())) {
        t0 <- pta_full[1] + plot_trim_ms
        t1 <- pta_full[length(pta_full)] - plot_trim_ms
        fidx <- which(pta_full >= t0 & pta_full <= t1)
        if (length(fidx) == 0) fidx <- seq_along(pta_full)
        all_filt <- do.call(c, lapply(res, function(rr) rr$filtered_subtracted[fidx]))
        ymin <- min(c(all_raw, all_filt), na.rm = TRUE)
        ymax <- max(c(all_raw, all_filt), na.rm = TRUE)
      } else {
        fidx <- seq_along(pta_full)
        ymin <- min(all_raw, na.rm = TRUE)
        ymax <- max(all_raw, na.rm = TRUE)
      }
      ypad <- (ymax - ymin) * 0.1
      if (!is.finite(ypad) || ypad == 0) ypad <- 1e-6

      # FIXED x range
      if (!is.null(xlim_vals)) {
        xlim <- xlim_vals
      } else {
        xlim <- c(5 * floor(min(pta_full) / 5), 5 * ceiling(max(pta_full) / 5))
      }
      status <- if (r$is_failure) 'failure' else 'success'

      display_trace <- r$raw_subtracted

      # draw full raw trace
      plot(pta_full, display_trace,
        type = 'l', col = adjustcolor('black', 0.3), lwd = lwd * 0.6,
        xlim = xlim, ylim = c(ymin - ypad, ymax + ypad),
        xlab = '', ylab = '', bty = 'n', axes = FALSE,
        main = sprintf('trace %d / %d  [%s]', idx, n_traces(), status), cex.main = 1.3)

      if (isTRUE(filter_on())) {
        lines(pta_full[fidx], r$filtered_subtracted[fidx], col = '#3c8dbc', lwd = lwd)
      }

      abline(h = 0, col = 'grey', lty = 3, lwd = lwd * 0.8)
      abline(v = stim, col = '#CD5C5C', lty = 2, lwd = lwd * 0.8)

      # peak marker
      pk_x <- NULL; pk_y <- NULL
      if (!r$is_failure) {
        pk_x <- stim + r$peak_time_ms
        pk_y <- r$peak_amplitude
      }

      # get fit data BEFORE referencing it
      fs <- fit_store()
      fit_out <- fs[[as.character(idx)]]

      # legend setup
      leg_lab <- 'raw'; leg_col <- adjustcolor('black', 0.3)
      leg_lwd <- lwd; leg_lty <- 1
      if (isTRUE(filter_on())) {
        leg_lab <- c(leg_lab, 'filtered'); leg_col <- c(leg_col, '#3c8dbc')
        leg_lwd <- c(leg_lwd, lwd); leg_lty <- c(leg_lty, 1)
      }

      # fit overlay + peak marker (single block)
      if (!is.null(fit_out) && !is.null(fit_out$traces)) {
        fp <- fit_out$fit_params
        st <- if (!is.null(fp$fit_stim_time)) fp$fit_stim_time else input$fit_stim_time
        bl <- if (!is.null(fp$fit_baseline)) fp$fit_baseline else input$fit_baseline          
        fit_x_orig <- fit_out$traces$x
        fit_dx <- fit_x_orig[2] - fit_x_orig[1]
        ext_start <- xlim[1] - (st - bl)
        ext_end <- xlim[2] - (st - bl)
        if (ext_end > ext_start && fit_dx > 0) {
          extended_x <- seq(ext_start, ext_end, by = fit_dx)
          fits_adj <- fit_out$fits
          fits_adj[4] <- fits_adj[4] + bl
          extended_y <- product1N(fits_adj, extended_x + fit_dx, N = 1, IEI = 50)
          extended_plot_x <- extended_x + (st - bl)
          in_r <- extended_plot_x >= xlim[1] & extended_plot_x <= xlim[2]
          lines(extended_plot_x[in_r], extended_y[in_r], col = 'indianred', lty = 3, lwd = lwd * 2)
          pk_idx <- which.max(extended_y[in_r])
          pk_x <- extended_plot_x[in_r][pk_idx]
          pk_y <- extended_y[in_r][pk_idx]
          leg_lab <- c(leg_lab, 'fit'); leg_col <- c(leg_col, 'indianred')
          leg_lwd <- c(leg_lwd, lwd); leg_lty <- c(leg_lty, 3)
        }
      }
      if (!is.null(pk_x)) points(pk_x, pk_y, pch = 16, col = 'indianred', cex = 1.5)
      legend('topleft', legend = leg_lab, col = leg_col,
        lwd = as.numeric(leg_lwd), lty = as.numeric(leg_lty), bty = 'n', cex = 1.1)

      # area window horizontal line
      usr <- par('usr')
      text(stim + (usr[2] - usr[1]) * 0.01, usr[4] - (usr[4] - usr[3]) * 0.02,
        'stim', col = 'indianred', cex = 1, adj = c(0, 1))
      text(input$fit_limit_ms + (usr[2] - usr[1]) * 0.01, usr[4] - (usr[4] - usr[3]) * 0.02,
        'fit limit', col = 'indianred', cex = 1, adj = c(0, 1))
      aw_y <- usr[3] + (usr[4] - usr[3]) * 0.08
      aw_x0 <- stim + input$ana_area_start
      aw_x1 <- stim + input$ana_area_end
      segments(aw_x0, aw_y, aw_x1, aw_y, col = 'steelblue', lwd = lwd * 1.5)
      text((aw_x0 + aw_x1) / 2, aw_y - (usr[4] - usr[3]) * 0.03,
        'area window', col = 'steelblue', cex = 1, adj = c(0.5, 1))
      text((aw_x0 + aw_x1) / 2, aw_y - (usr[4] - usr[3]) * 0.06,
        '(main analysis)', col = 'steelblue', cex = 1, adj = c(0.5, 1))

      # fit limit vertical line
      abline(v = input$fit_limit_ms, col = 'indianred', lty = 2, lwd = lwd * 0.8)

      # scale bars
      sb_y <- usr[3] + (usr[4] - usr[3]) / 20
      sb_x <- usr[2] - xbar - (usr[2] - usr[1]) * 0.05
      segments(sb_x, sb_y, sb_x + xbar, sb_y, lwd = lwd)
      segments(sb_x, sb_y, sb_x, sb_y + ybar, lwd = lwd)
      text((2 * sb_x + xbar) / 2, sb_y - (usr[4] - usr[3]) * 0.03,
        paste(xbar, xbar_lab), adj = c(0.5, 1), cex = 0.8)
      text(sb_x - (usr[2] - usr[1]) * 0.02, (2 * sb_y + ybar) / 2,
        paste(ybar, ybar_lab), srt = 90, adj = c(0.5, 0.5), cex = 0.8)
    })

    # =========================================================================
    # TAB 7 Average
    # =========================================================================

    observeEvent(input$compute_avg_btn, {
      res <- analysis_results()
      if (is.null(res)) { showNotification('Run analysis first.', type = 'warning'); return() }
      fs <- fit_store()
      n_pts <- length(res[[1]]$raw_subtracted)

      # build matrix of baseline-subtracted traces
      # each trace recomputed with its own fit baseline if available
      all_mat <- matrix(NA, nrow = n_pts, ncol = length(res))
      for (i in seq_along(res)) {
        r <- res[[i]]
        fo <- fs[[as.character(i)]]
        if (!is.null(fo) && !is.null(fo$fit_params)) {
          bl_ms <- fo$fit_params$fit_baseline
          st_ms <- fo$fit_params$fit_stim_time
          bl_start <- max(1, round((st_ms - bl_ms) / dt_ms()) + 1)
          bl_end <- max(bl_start, round(st_ms / dt_ms()))
          bl_end <- min(bl_end, n_pts)
          orig <- r$raw_subtracted + r$baseline_mean
          bl_mean <- mean(orig[bl_start:bl_end])
          all_mat[, i] <- orig - bl_mean
        } else {
          all_mat[, i] <- r$raw_subtracted
        }
      }

      success_idx <- which(sapply(res, function(r) !r$is_failure))
      success_mat <- if (length(success_idx) > 0) all_mat[, success_idx, drop = FALSE] else NULL

      ci_level <- input$avg_ci_level
      ci_method <- input$avg_ci_method
      nboot <- if (!is.null(input$avg_nboot)) input$avg_nboot else 999

      withProgress(message = 'Computing average...', value = 0.5, {
        all_ci <- compute_ci_psp(all_mat, ci_method, ci_level, nboot)
        success_ci <- if (!is.null(success_mat) && ncol(success_mat) > 1) {
          compute_ci_psp(success_mat, ci_method, ci_level, nboot)
        } else NULL
        setProgress(1)
      })
      avg_data(list(all = all_ci, success = success_ci))
      showNotification('Average computed.', type = 'message')
    })

    output$avg_info_text <- renderPrint({
      ad <- avg_data()
      if (is.null(ad)) { cat('Click Compute Average.\n'); return() }
      cat(sprintf('All traces: n = %d\n', ad$all$n))
      if (!is.null(ad$success)) cat(sprintf('Successes:  n = %d\n', ad$success$n))
      cat(sprintf('CI method: %s (%.0f%%)\n', input$avg_ci_method, input$avg_ci_level * 100))
    })

    output$avg_plot <- renderPlot({
      ad <- avg_data()
      if (is.null(ad)) {
        plot.new()
        text(0.5, 0.5, 'Click Compute Average.', cex = 1.5, col = 'grey50')
        return()
      }
      pta <- plot_time()
      lwd <- as.numeric(input$lwd)
      xbar <- as.numeric(input$xbar)
      ybar <- as.numeric(input$ybar)
      xbar_lab <- input$xbar_lab
      ybar_lab <- input$ybar_lab
      shade_alpha <- if (!is.null(input$avg_alpha)) input$avg_alpha else 0.3
      xlim_vals <- if (nchar(input$xlim) > 0) {
        v <- as.numeric(unlist(strsplit(input$xlim, ',')))
        if (length(v) == 2 && !any(is.na(v))) v else NULL
      } else NULL
      if (!is.null(xlim_vals)) { xlim <- xlim_vals
      } else { xlim <- c(5 * floor(min(pta) / 5), 5 * ceiling(max(pta) / 5)) }

      # use same y scale as individual trace plots
      res <- analysis_results()
      all_raw <- do.call(c, lapply(res, function(rr) rr$raw_subtracted))
      ymin <- min(all_raw, na.rm = TRUE)
      ymax <- max(all_raw, na.rm = TRUE)
      ypad <- (ymax - ymin) * 0.1
      if (!is.finite(ypad) || ypad == 0) ypad <- 1e-6
      shared_ylim <- c(ymin - ypad, ymax + ypad)

      # helper to draw one average panel
      draw_avg_panel <- function(ci_data, title_text) {
        avg <- ci_data$avg
        lo <- ci_data$lo
        hi <- ci_data$hi
        
        ylim <- shared_ylim
        
        plot(pta, avg, type = 'n', xlim = xlim, ylim = ylim,
          xlab = '', ylab = '', bty = 'n', axes = FALSE,
          main = title_text, cex.main = 1.3)
        polygon(c(pta, rev(pta)), c(lo, rev(hi)),
          col = adjustcolor('steelblue', alpha.f = shade_alpha), border = NA)
        lines(pta, avg, col = 'steelblue', lwd = lwd * 1.5)
        abline(h = 0, col = 'grey', lty = 3, lwd = lwd * 0.8)
        abline(v = loaded_stim(), col = '#CD5C5C', lty = 2, lwd = lwd * 0.8)
        # scale bars
        usr <- par('usr')
        sb_y <- usr[3] + (usr[4] - usr[3]) / 20
        sb_x <- usr[2] - xbar - (usr[2] - usr[1]) * 0.05
        segments(sb_x, sb_y, sb_x + xbar, sb_y, lwd = lwd)
        segments(sb_x, sb_y, sb_x, sb_y + ybar, lwd = lwd)
        text((2 * sb_x + xbar) / 2, sb_y - (usr[4] - usr[3]) * 0.03,
          paste(xbar, xbar_lab), adj = c(0.5, 1), cex = 0.8)
        text(sb_x - (usr[2] - usr[1]) * 0.02, (2 * sb_y + ybar) / 2,
          paste(ybar, ybar_lab), srt = 90, adj = c(0.5, 0.5), cex = 0.8)
      }

      par(mfrow = c(2, 1), mar = c(3, 3, 3, 2))

      # top panel: successes only
      if (!is.null(ad$success)) {
        draw_avg_panel(ad$success,
          sprintf('Successes (n=%d) - mean +/- %.0f%% CI',
            ad$success$n, input$avg_ci_level * 100))
      } else {
        plot.new()
        text(0.5, 0.5, 'No successes to average.', cex = 1.2, col = 'grey50')
      }

      # bottom panel: all traces
      draw_avg_panel(ad$all,
        sprintf('All traces (n=%d) - mean +/- %.0f%% CI',
          ad$all$n, input$avg_ci_level * 100))
    })

    # =========================================================================
    # save on session end
    # =========================================================================

    session$onSessionEnded(function() {
      final <- list(
        results     = isolate(analysis_results()),
        fit_store   = isolate(fit_store()),
        traces      = isolate(loaded_traces()),
        sample_rate = isolate(loaded_sr()),
        stimulation = isolate(loaded_stim())
      )
      assign(out_name, final, envir = .GlobalEnv)
    })
  }

  shinyApp(ui, server)
}


# ==============================================
# EXPORT TO EXCEL
# ==============================================

export_results_to_excel <- function(results, time, output_path, filename) {
  
  traces_df <- data.frame(
    lapply(seq_along(results), function(i) {
      results[[i]]$raw_subtracted
    })
  )
  names(traces_df) <- paste0('trace', seq_along(results))
  
  success_indices <- which(sapply(results, function(r) !r$is_failure))
  success_df <- data.frame(success_trace_idx = success_indices)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 'raw_traces')
  openxlsx::addWorksheet(wb, 'success_map')
  
  openxlsx::writeData(wb, 'raw_traces', traces_df)
  openxlsx::writeData(wb, 'success_map', success_df)
  
  full_path <- file.path(output_path, paste0(filename, '.xlsx'))
  openxlsx::saveWorkbook(wb, full_path, overwrite = TRUE)
  
  cat(sprintf('Saved to %s\n', full_path))
  cat(sprintf('  - raw_traces: %d traces\n', length(results)))
  cat(sprintf('  - success_map: %d successful indices\n', length(success_indices)))
}

.alphaize <- function(cols, alpha=1) {
  m <- grDevices::col2rgb(cols)
  grDevices::rgb(m[1,], m[2,], m[3,], alpha*255, maxColorValue=255)
}

.palette_cols <- function(name, n, alpha=1, reverse=FALSE) {
  valid <- c('Vik','jet','viridis','cividis','PuOr','BrBG','Roma','Batlow','Berlin')
  idx <- pmatch(tolower(name[1]), tolower(valid), nomatch = NA)
  if (is.na(idx)) stop('color_palette must be one of: ', paste(valid, collapse=', '))
  key_lc <- tolower(valid[idx])
  funs <- list(
    jet   = jet.colors,    viridis = viridis.colors, cividis = cividis.colors,
    puor  = puor.colors,   brbg    = brbg.colors,    roma    = roma.colors,
    vik   = vik.colors,    batlow  = batlow.colors,  berlin  = berlin.colors
  )
  fun <- funs[[key_lc]]
  if (!is.function(fun)) stop('palette function not found for: ', valid[idx])
  cols <- fun(n, alpha = alpha)
  if (reverse) rev(cols) else cols
}

jet.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow','#FF7F00','red','#7F0000')
  )(n)
  .alphaize(cols, alpha)
}

viridis.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#440154','#3b528b','#21908d','#5ec962','#fde725')
  )(n)
  .alphaize(cols, alpha)
}

cividis.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#00204c','#414487','#7e03a8','#a5f86b','#ffffc0')
  )(n)
  .alphaize(cols, alpha)
}

puor.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#7f3b08','#fdb863','#f7f7f7','#b2abd2','#5e3c99')
  )(n)
  .alphaize(cols, alpha)
}

brbg.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#543005','#bf812d','#f6e8c3','#c7eae5','#35978f','#003c30')
  )(n)
  .alphaize(cols, alpha)
}

roma.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#7E1700','#8C360A','#984F13','#A3651E','#AD7B27','#B99235','#C4AB4A','#CEC56C','#D1DC94','#C9E8B5',
      '#B4E9CC','#95E0D6','#73CED5','#54B8D0','#3DA2C9','#2F8CBF','#2677B7','#1F61AD','#164AA3','#023198')
  )(n)
  .alphaize(cols, alpha)
}

vik.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#001260','#01276E','#023B7B','#055189','#166898','#3983AB','#629FBD','#8CB9CF','#B7D3E1','#DFE5E9',
      '#EDDCD2','#E5C2AE','#DAA78A','#CF8F69','#C6774A','#BA5F2C','#A5400F','#882506','#6F1107','#590007')
  )(n)
  .alphaize(cols, alpha)
}

batlow.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#001959','#0B2D5D','#103D5E','#134B61','#195661','#25625F','#366A58','#49714E','#5F7842','#757E37',
      '#8E852D','#A88B2B','#C29037','#D9954A','#ED9A63','#F9A382','#FDADA1','#FDB7BE','#FCC1DB','#F9CCF9')
  )(n)
  .alphaize(cols, alpha)
}

berlin.colors <- function(n, alpha=1) {
  cols <- grDevices::colorRampPalette(
    c('#9EB0FF','#7FABF0','#5DA5DD','#4093C0','#307A9E','#25617E','#1C4960','#153342','#101F27','#121112',
      '#200A03','#2F0E00','#421300','#571B06','#742C16','#904430','#AB5D4E','#C6776C','#E2928C','#FFACAC')
  )(n)
  .alphaize(cols, alpha)
}

plot_abf_channels <- function(
  abf, sweep = 1, channels = NULL, channel_colors = NULL,
  height_ratios = NULL, ylims = NULL, xlim = c(0, 50),
  font_family = 'Helvetica Neue', font_size = 10, line_width = 1,
  figsize = c(10, 9), trim_axes = TRUE,
  save = FALSE, filename = 'abf_channels', save_path = NULL
) {

  # resolve channels (1-based)
  n_channels <- length(abf$channelNames)
  if (is.null(channels)) channels <- seq_len(n_channels)
  n <- length(channels)

  # pick sweep matrix (samples × channels)
  d <- abf$data[[sweep]]

  # time axis
  n_samples <- nrow(d)
  x <- seq(0, by = abf$samplingIntervalInSec, length.out = n_samples)

  # defaults
  default_cols <- c('#1a1a1a', '#4682B4', '#6A5ACD', '#E8524A')
  if (is.null(channel_colors)) {
    channel_colors <- default_cols[((seq_len(n) - 1) %% length(default_cols)) + 1]
  }
  if (is.null(ylims)) ylims <- vector('list', n)
  if (is.null(height_ratios)) height_ratios <- rep(1, n)

  hr <- height_ratios / sum(height_ratios)

  gap <- 0.03
  margin_bottom <- 0.08
  margin_top <- 0.02
  margin_left <- 0.10
  margin_right <- 0.03
  total_plot <- 1 - margin_bottom - margin_top - gap * (n - 1)
  panel_heights <- hr * total_plot

  if (save) {
    full_path <- if (!is.null(save_path)) {
      file.path(save_path, paste0(filename, '.svg'))
    } else {
      paste0(filename, '.svg')
    }
    svg(full_path, width = figsize[1], height = figsize[2],
        bg = 'transparent', family = font_family)
    on.exit(dev.off(), add = TRUE)
  }

  par(family = font_family, cex = font_size / 10, lwd = line_width)

  for (i in seq_len(n)) {
    ch <- channels[i]
    y <- d[, ch] 

    y_top <- 1 - margin_top - sum(panel_heights[seq_len(i - 1)]) - gap * (i - 1)
    y_bot <- y_top - panel_heights[i]

    par(fig = c(margin_left, 1 - margin_right, y_bot, y_top),
        new = if (i > 1) TRUE else FALSE,
        mar = c(0, 4, 0, 0),
        mgp = c(2.5, 0.5, 0),
        tcl = -0.15)

    if (!is.null(ylims[[i]])) {
      yl <- ylims[[i]]
    } else {
      yr <- range(y, na.rm = TRUE)
      ypad <- diff(yr) * 0.05
      yl <- yr + c(-ypad, ypad)
    }

    xl <- if (!is.null(xlim)) xlim else range(x)
    in_range <- x >= xl[1] & x <= xl[2]

    plot(x[in_range], y[in_range],
         type = 'l', col = channel_colors[i], lwd = line_width,
         xlim = xl, ylim = yl,
         xlab = '', ylab = '',
         axes = FALSE, xaxs = 'i', yaxs = 'i')

    mtext(paste0(abf$channelNames[ch], '\n(', abf$channelUnits[ch], ')'),
          side = 2, line = 2.5, cex = font_size / 12)

    y_ticks <- pretty(yl)
    y_ticks <- y_ticks[y_ticks >= yl[1] & y_ticks <= yl[2]]
    axis(2, at = y_ticks, las = 1, lwd = 0, lwd.ticks = line_width,
         cex.axis = font_size / 12)
    if (trim_axes && length(y_ticks) > 0) {
      segments(par('usr')[1], min(y_ticks),
               par('usr')[1], max(y_ticks),
               lwd = line_width, xpd = TRUE)
    } else {
      axis(2, at = yl, labels = FALSE, lwd = line_width, lwd.ticks = 0)
    }

    if (i == n) {
      x_ticks <- pretty(xl)
      x_ticks <- x_ticks[x_ticks >= xl[1] & x_ticks <= xl[2]]
      axis(1, at = x_ticks, lwd = 0, lwd.ticks = line_width,
           cex.axis = font_size / 12)
      mtext('Time (s)', side = 1, line = 2, cex = font_size / 12)
      if (trim_axes && length(x_ticks) > 0) {
        segments(min(x_ticks), par('usr')[3],
                 max(x_ticks), par('usr')[3],
                 lwd = line_width, xpd = TRUE)
      }
    }
  }

  if (save) message('saved to ', full_path)
  invisible(NULL)
}


# Multivariate outlier detection
# Rousseeuw, P. J., & van Zomeren, B. C. (1990). Unmasking multivariate outliers and leverage points. Journal of the American Statistical Association, 85(411), 633–639
# Rousseeuw, P. J., & Van Driessen, K. (1999). A fast algorithm for the minimum covariance determinant estimator. Technometrics, 41(3), 212–223
# 'MD'  - classical Mahalanobis distance (masked by outliers)
# 'MVE' - Robust distance from Minimum Volume Ellipsoid
# 'MCD' - Robust distance from Minimum Covariance Determinant
# Under multivariate normality, MD^2 ~ chi^2_p, giving the cutoff sqrt(qchisq(quant, df=p))

mv_outliers <- function(mat, method=c('MCD', 'MVE', 'MD'), alpha=0.5, quant=0.975, plot=FALSE, 
                        type=c('xy','dd','both'), palette=c('Roma','viridis','jet','cividis','PuOr','BrBG','Vik','Batlow','Berlin'),
                        xlab=NULL, ylab=NULL, lwd=0.8, width=5, height=5, bg='transparent', filename='outlier_plot.svg', save=FALSE) {

  method  <- match.arg(method)
  type    <- match.arg(type)
  palette <- match.arg(palette)
  if (method %in% c('MCD', 'MVE') &&
      !requireNamespace('robustbase', quietly=TRUE))
    stop("Package 'robustbase' is required for MCD / MVE methods")

  is_num  <- vapply(mat, is.numeric, logical(1))
  xm      <- as.matrix(mat[, is_num, drop=FALSE])
  grp_col <- mat[, !is_num, drop=FALSE]
  group   <- if (ncol(grp_col) > 0) grp_col[, 1] else NULL

  if (anyNA(xm)) stop("x contains NA values - remove or impute first")
  n  <- nrow(xm); p <- ncol(xm)
  if (n <= 2 * p) warning("n <= 2*p: robust estimates may be unstable")

  nms    <- if (!is.null(rownames(xm))) rownames(xm) else as.character(seq_len(n))
  cutoff <- sqrt(qchisq(quant, df=p))

  MD <- sqrt(mahalanobis(xm, center=colMeans(xm), cov=cov(xm)))

  RD <- NULL
  if (method == 'MCD') {
    fit <- tryCatch(robustbase::covMcd(xm, alpha=alpha),
                    error=function(e) NULL)
    if (is.null(fit)) stop("covMcd() failed - check sample size and collinearity")
    RD <- sqrt(mahalanobis(xm, fit$center, fit$cov))
  } else if (method == 'MVE') {
    fit <- tryCatch(robustbase::covMve(xm, alpha=alpha),
                    error=function(e) NULL)
    if (is.null(fit)) stop("covMve() failed - check sample size and collinearity")
    RD <- sqrt(mahalanobis(xm, fit$center, fit$cov))
  }

  dist_used <- if (method == 'MD') MD else RD
  out <- dist_used > cutoff
  names(out) <- nms

  attr(out, 'method')  <- method
  attr(out, 'cutoff')  <- cutoff
  attr(out, 'quant')   <- quant
  attr(out, 'MD')      <- setNames(MD, nms)
  attr(out, 'RD')      <- if (!is.null(RD)) setNames(RD, nms)
  attr(out, 'x')       <- xm
  attr(out, 'group')   <- if (!is.null(group)) setNames(as.character(group), nms)
  attr(out, 'palette') <- palette

  if (plot) {
    if (type == 'xy') {
      xy_plot(out, xlab=xlab, ylab=ylab, lwd=lwd,
              width=width, height=height, bg=bg, filename=filename, save=save)
    } else if (type == 'dd') {
      dd_plot(out, lwd=lwd,
              width=width, height=height, bg=bg, filename=filename, save=save)
    } else {
      dev.new(width=width*2, height=height, noRStudioGD=TRUE)
      par(bg=NA, mfrow=c(1,2))
      .draw_xy(out, xlab=xlab, ylab=ylab, lwd=lwd)
      .draw_dd(out, lwd=lwd)
      par(mfrow=c(1,1))
      if (save) dev.copy(svg, file=filename, width=width*2, height=height)
      if (save) dev.off()
    }
  }
  
  out
}

.group_cols <- function(outlier) {
  group   <- attr(outlier, 'group')
  palette <- attr(outlier, 'palette')
  if (is.null(group)) return(list(pt_col=rep('black', length(outlier)), lvls=NULL, cols=NULL))
  lvls   <- unique(group)
  cols   <- .palette_cols(palette, length(lvls))
  pt_col <- cols[match(group, lvls)]
  list(pt_col=pt_col, lvls=lvls, cols=cols)
}

.draw_xy <- function(outlier, xlab=NULL, ylab=NULL, lwd=0.8) {
  xm     <- attr(outlier, 'x')
  method <- attr(outlier, 'method')
  group  <- attr(outlier, 'group')
  xv     <- xm[, 1]
  yv     <- xm[, 2]
  xlab   <- if (is.null(xlab)) colnames(xm)[1] else xlab
  ylab   <- if (is.null(ylab)) colnames(xm)[2] else ylab

  ticks_x <- pretty(range(xv), n=5)
  ticks_y <- pretty(range(yv), n=5)
  lim_x   <- range(ticks_x)
  lim_y   <- range(ticks_y)

  gc     <- .group_cols(outlier)
  pt_col <- gc$pt_col

  plot(xv, yv, type='n', xlim=lim_x, ylim=lim_y, xlab=xlab, ylab=ylab, main=method,
       bty='n', axes=FALSE, las=1)
  abline(0, 1, col='darkgray', lty=3, lwd=lwd)

  points(xv[!outlier], yv[!outlier], pch=19, cex=0.5,
         col=adjustcolor(pt_col[!outlier], alpha.f=0.6))
  points(xv[ outlier], yv[ outlier], pch=19, cex=0.75,
         col=adjustcolor(pt_col[outlier], alpha.f=0.8))
  points(xv[ outlier], yv[ outlier], pch=1,  cex=1.1,
         col='indianred', lwd=1.2)
  if (any(outlier))
    text(xv[outlier], yv[outlier], labels=names(outlier)[outlier],
         pos=4, col='indianred', cex=0.75)

  if (!is.null(gc$lvls))
    text(rep(lim_x[1], length(gc$lvls)),
         lim_y[2] - (seq_along(gc$lvls) - 1) * diff(lim_y) * 0.06,
         labels=gc$lvls, col=gc$cols, adj=0, cex=0.75)

  axis(1, at=ticks_x, tcl=-0.2, las=1, lwd=lwd)
  axis(2, at=ticks_y, tcl=-0.2, las=1, lwd=lwd)
}

.draw_dd <- function(outlier, lwd=0.8) {
  MD     <- attr(outlier, 'MD')
  RD     <- attr(outlier, 'RD')
  cut    <- attr(outlier, 'cutoff')

  ticks  <- pretty(c(MD, RD, cut), n=5)
  lim    <- range(ticks)

  gc     <- .group_cols(outlier)
  pt_col <- gc$pt_col

  plot(MD, RD, type='n', xlim=lim, ylim=lim,
       xlab='Mahalanobis distance', ylab='Robust distance', main='Distance-Distance',
       bty='n', axes=FALSE, las=1)
  abline(0, 1,  col='darkgray',  lty=3, lwd=lwd)
  abline(h=cut, col='indianred', lty=3, lwd=lwd)
  abline(v=cut, col='indianred', lty=3, lwd=lwd)

  points(MD[!outlier], RD[!outlier], pch=19, cex=0.5,
         col=adjustcolor(pt_col[!outlier], alpha.f=0.6))
  points(MD[ outlier], RD[ outlier], pch=19, cex=0.75,
         col=adjustcolor(pt_col[outlier], alpha.f=0.8))
  points(MD[ outlier], RD[ outlier], pch=1,  cex=1.1,
         col='indianred', lwd=1.2)
  if (any(outlier))
    text(MD[outlier], RD[outlier], labels=names(outlier)[outlier],
         pos=4, col='indianred', cex=0.75)

  if (!is.null(gc$lvls))
    text(rep(lim[1], length(gc$lvls)),
         lim[2] - (seq_along(gc$lvls) - 1) * diff(lim) * 0.06,
         labels=gc$lvls, col=gc$cols, adj=0, cex=0.75)

  axis(1, at=ticks, tcl=-0.2, las=1, lwd=lwd)
  axis(2, at=ticks, tcl=-0.2, las=1, lwd=lwd)
}


# Scatter plot of data coloured by outlier status (Rousseeuw & van Zomeren 1990)
xy_plot <- function(outlier, xlab=NULL, ylab=NULL, lwd=0.8,
                    width=5, height=5, bg='transparent', filename='outlier_plot.svg', save=FALSE) {
  dev.new(width=width, height=height, noRStudioGD=TRUE)
  par(bg=NA)
  .draw_xy(outlier, xlab=xlab, ylab=ylab, lwd=lwd)
  if (save) dev.copy(svg, file=filename, width=width, height=height)
  if (save) dev.off()
  invisible(NULL)
}



# Distance-distance plot: robust vs classical Mahalanobis distance
dd_plot <- function(outlier, lwd=0.8,
                    width=5, height=5, bg='transparent', filename='dd_plot.svg', save=FALSE) {
  if (is.null(attr(outlier, 'RD')))
    stop("No robust distances; rerun with method = 'MCD' or 'MVE'")
  dev.new(width=width, height=height, noRStudioGD=TRUE)
  par(bg=NA)
  .draw_dd(outlier, lwd=lwd)
  if (save) dev.copy(svg, file=filename, width=width, height=height)
  if (save) dev.off()
  invisible(NULL)
}

# gets animal id from file name
animal_id_fun <- function(data_list) {
  animal_offset <- 0
  animal_map <- lapply(names(data_list), function(data_name) {
    expt.id <- colnames(data_list[[data_name]])
    # remove the final three-digit experiment number
    animal_date <- substr(expt.id, 1, nchar(expt.id) - 3)
    animal.id <- match(animal_date, unique(animal_date)) + animal_offset
    animal_offset <<- max(animal.id)
    data.frame(
      data=data_name,
      expt.id=expt.id,
      animal.id=animal.id
    )
  })
  animal_map <- do.call(rbind, animal_map)
  rownames(animal_map) <- NULL
  animal_map
}

# simple output for CR2 stats test
CR2output <- function(test_result, formula, df, cluster='animal') {

  var.names <- all.vars(formula)

  parameter <- var.names[1]
  group <- var.names[2]

  group.levels <- levels(droplevels(factor(df[[group]])))

  if (length(group.levels)!=2) {
    stop('CR2output requires exactly two groups')
  }

  cr2.row <- test_result[test_result$Coef!='(Intercept)',]

  if (nrow(cr2.row)!=1) {
    stop('CR2output requires exactly one group coefficient')
  }

  group.n <- table(factor(df[[group]], levels=group.levels))

  output <- data.frame(
    parameter=parameter,
    comparison=paste0('within ', group, ' (clustered by ', cluster, ')'),
    contrast=paste(group.levels, collapse=' vs '),
    n=paste(as.integer(group.n), collapse=' vs '),
    test='CR2 cluster-robust test with Satterthwaite correction',
    alternative='two.sided',
    `test stat`='t',
    stat=cr2.row$tstat,
    `p value`=cr2.row$p_Satt,
    `p adjusted`=cr2.row$p_Satt,
    check.names=FALSE
  )

  rownames(output) <- NULL

  output
}


# rlmer.bootstrap() fits a robust linear mixed-effects model and performs a stratified cluster bootstrap. 
# Complete clusters are resampled with replacement within their original two-level group, preserving the observed 
# group sizes and all observations belonging to each cluster.
# H0=TRUE: removes the fitted group effect, generates a null bootstrap distribution of the studentized coefficient, 
# and returns a two-sided bootstrap p-value.
# H0=FALSE: resamples the observed data without removing the group effect and returns the group difference, 
# percentile confidence interval and sign-tail probability.
# The function automatically interprets the factor contrasts, so the reported parameter is the second factor level 
# minus the first. Failed model fits are excluded and reported as Nfailed; successful bootstrap fits are returned in bootstrap.values.

rlmer.cluster.bootstrap <- function(formula, longdata, cluster, group, H0=TRUE,B=1e4,
              n_cores=1, seed=42, CI=0.95, data.name=''){

  group.contrasts <- contrasts(longdata[[group]])
  longdata[[group]] <- droplevels(longdata[[group]])
  longdata[[cluster]] <- droplevels(longdata[[cluster]])
  group.levels <- levels(longdata[[group]])

  group.contrasts <- group.contrasts[group.levels,,drop=FALSE]
  contrasts(longdata[[group]]) <- group.contrasts

  if (length(group.levels) != 2) {
    stop('group must contain exactly two levels')
  }

  group.contrasts <- contrasts(longdata[[group]])

  if (is.null(group.contrasts) || ncol(group.contrasts) != 1) {
    stop('group must have one two-level contrast')
  }

  rfit <- rlmer(formula, data=longdata)
  coef_name <- setdiff(names(fixef(rfit)), '(Intercept)')

  if (length(coef_name) != 1) {
    stop('the model must contain one non-intercept fixed-effect coefficient')
  }

  contrast.multiplier <- unname(group.contrasts[2,1]-group.contrasts[1,1])
  contrast.name <- paste(group.levels[2], '-', group.levels[1])

  if (H0) {
    X <- model.matrix(formula(rfit, fixed.only=TRUE), data=longdata)

    longdata.bootstrap <- longdata
    longdata.bootstrap[[all.vars(formula)[1]]] <- longdata[[all.vars(formula)[1]]]-
      drop(X[,coef_name,drop=FALSE]%*%fixef(rfit)[coef_name])
  } else {
    longdata.bootstrap <- longdata
  }

  cluster_list <- split(longdata.bootstrap, longdata.bootstrap[[cluster]])

  cluster_group <- vapply(cluster_list, function(x){
    value <- unique(as.character(x[[group]]))

    if (length(value) != 1) {
      stop('each cluster must belong to one group')
    }

    value
  }, character(1))

  cluster_list <- split(cluster_list, cluster_group)
  n_clusters <- lengths(cluster_list)

  if (any(n_clusters < 2)) {
    stop('each group must contain at least two clusters')
  }

  if (any(n_clusters < 10)) {
    warning(
      paste0('bootstrap calibration may be unreliable with fewer than 10 clusters in a group: ', paste(names(n_clusters), n_clusters, sep='=', collapse=', ')),
      call.=FALSE
    )
  }

  cluster.summary <- paste(names(n_clusters), n_clusters, sep='=', collapse=', ')

  bootstrap.fun <- function(ii){
    samp <- unlist(lapply(seq_along(cluster_list), function(kk){
      sample_idx <- sample(seq_along(cluster_list[[kk]]), size=n_clusters[kk], replace=TRUE)
      lapply(seq_along(sample_idx), function(j){
        tmp <- cluster_list[[kk]][[sample_idx[j]]]
        tmp[[cluster]] <- paste0('boot_', ii, '*', kk, '*', j)
        tmp
      })
    }), recursive=FALSE)

    dboot <- do.call(rbind, samp)
    dboot[[group]] <- factor(dboot[[group]], levels=group.levels)
    contrasts(dboot[[group]]) <- group.contrasts

    fit <- try(rlmer(formula, data=dboot), silent=TRUE)

    if (inherits(fit, 'try-error')) {
      return(c(estimate=NA_real_, statistic=NA_real_))
    }

    estimate <- contrast.multiplier*unname(fixef(fit)[coef_name])

    if (H0) {
      se <- sqrt(diag(vcov(fit)))[coef_name]

      if (!is.finite(se) || se <= 0) {
        return(c(estimate=NA_real_, statistic=NA_real_))
      }

      statistic <- sign(contrast.multiplier)*unname(fixef(fit)[coef_name]/se)
    } else {
      statistic <- NA_real_
    }

    c(estimate=estimate, statistic=statistic)
  }

  old.kind <- RNGkind()
  had.seed <- exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)
  old.seed <- if (had.seed) get('.Random.seed', envir=.GlobalEnv, inherits=FALSE) else NULL

  on.exit({
    do.call(RNGkind, as.list(old.kind))

    if (had.seed) {
      assign('.Random.seed', old.seed, envir=.GlobalEnv)
    } else if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)) {
      rm('.Random.seed', envir=.GlobalEnv)
    }
  }, add=TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  if (n_cores == 1) {
    bootstrap.results <- lapply(seq_len(B), bootstrap.fun)
  } else {
    bootstrap.results <- parallel::mclapply(seq_len(B), bootstrap.fun, mc.cores=n_cores)
  }

  bootstrap.results <- do.call(rbind, bootstrap.results)

  if (H0) {
    successful <- is.finite(bootstrap.results[, 'estimate'])&
      is.finite(bootstrap.results[, 'statistic'])
  } else {
    successful <- is.finite(bootstrap.results[, 'estimate'])
  }

  bootstrap.results <- bootstrap.results[successful,,drop=FALSE]
  bootstrap.estimates <- bootstrap.results[, 'estimate']

  if (H0) {
    bootstrap.values <- bootstrap.results[, 'statistic']
  } else {
    bootstrap.values <- bootstrap.estimates
  }

  observed.estimate <- contrast.multiplier*unname(fixef(rfit)[coef_name])
  target <- if (H0) 0 else observed.estimate
  alpha <- 1-CI
  ci <- quantile(bootstrap.estimates, probs=c(alpha/2, 1-alpha/2))

  summary.boot <- data.frame(
    parameter=contrast.name,
    target=target,
    mean=mean(bootstrap.estimates),
    bias=mean(bootstrap.estimates)-target,
    se=sd(bootstrap.estimates),
    ci.lower=unname(ci[1]),
    ci.upper=unname(ci[2]),
    Nboot=length(bootstrap.estimates),
    Nfailed=B-length(bootstrap.estimates),
    check.names=FALSE
  )

  rownames(summary.boot) <- NULL
  model.name <- paste(deparse(formula), collapse=' ')

  if (H0) {
    se.observed <- sqrt(diag(vcov(rfit)))[coef_name]
    statistic.observed <- sign(contrast.multiplier)*unname(fixef(rfit)[coef_name]/se.observed)

    pval <- (sum(abs(bootstrap.values) >= abs(statistic.observed))+1)/(length(bootstrap.values)+1)

    summary <- data.frame(
      model=paste0('rlmer null-centred stratified cluster bootstrap: ', model.name),
      data=data.name,
      parameter=contrast.name,
      statistic=statistic.observed,
      'p value'=pval,
      clusters=cluster.summary,
      Nboot=length(bootstrap.values),
      Nfailed=B-length(bootstrap.values),
      n_cores=n_cores,
      check.names=FALSE
    )
  } else {
    p_negative <- mean(bootstrap.estimates < 0)
    p_positive <- mean(bootstrap.estimates > 0)
    tail.probability <- 2*min(p_negative, p_positive)

    summary <- data.frame(
      model=paste0('rlmer stratified cluster bootstrap: ', model.name),
      data=data.name,
      parameter=contrast.name,
      estimate=observed.estimate,
      'CI (2.5)'=unname(ci[1]),
      'CI (97.5)'=unname(ci[2]),
      'sign-tail probability'=tail.probability,
      clusters=cluster.summary,
      Nboot=length(bootstrap.estimates),
      Nfailed=B-length(bootstrap.estimates),
      n_cores=n_cores,
      check.names=FALSE
    )
  }

  rownames(summary) <- NULL

  list(
    fit=rfit,
    summary=summary,
    summary.boot=summary.boot,
    bootstrap.values=bootstrap.values,
    bootstrap.estimates=bootstrap.estimates,
    n.clusters=n_clusters,
    H0=H0
  )
}

# Exact stratified cluster permutation test.
# rlmer.cluster.permutation() fits a robust linear mixed-effects model and performs a cluster-level permutation test.
# Complete clusters are reassigned between the two group levels while preserving the number of clusters in each group
# and, when supplied, preserving the group allocation within each permutation stratum.
# The function refits the model for every permitted allocation and compares the absolute observed studentized coefficient
# with the corresponding permutation distribution to calculate a two-sided permutation p-value.
# Exact enumeration is used when the number of allocations is manageable; otherwise, group allocations are sampled.
# The function automatically interprets the factor contrasts, so the reported parameter is the second factor level
# minus the first. The number of permutations, minimum attainable p-value and null distribution are returned.

rlmer.cluster.permutation <- function(formula, longdata, cluster, group, strata=NULL,
              n_cores=1, CI=0.95, max.permutations=1e4, data.name=''){

  group.contrasts <- contrasts(longdata[[group]])
  longdata[[group]] <- droplevels(longdata[[group]])
  longdata[[cluster]] <- droplevels(longdata[[cluster]])
  group.levels <- levels(longdata[[group]])

  group.contrasts <- group.contrasts[group.levels,,drop=FALSE]
  contrasts(longdata[[group]]) <- group.contrasts

  if (length(group.levels) != 2) {
    stop('group must contain exactly two levels')
  }

  group.contrasts <- contrasts(longdata[[group]])

  if (is.null(group.contrasts) || ncol(group.contrasts) != 1) {
    stop('group must have one two-level contrast')
  }

  rfit <- rlmer(formula, data=longdata)
  coef_name <- setdiff(names(fixef(rfit)), '(Intercept)')

  if (length(coef_name) != 1) {
    stop('the model must contain one non-intercept fixed-effect coefficient')
  }

  contrast.multiplier <- unname(group.contrasts[2,1]-group.contrasts[1,1])
  contrast.name <- paste(group.levels[2], '-', group.levels[1])

  cluster_list <- split(longdata, longdata[[cluster]])
  cluster.ids <- names(cluster_list)

  cluster_group <- vapply(cluster_list, function(x){
    value <- unique(as.character(x[[group]]))

    if (length(value) != 1) {
      stop('each cluster must belong to one group')
    }

    value
  }, character(1))

  n_clusters <- table(factor(cluster_group, levels=group.levels))

  if (any(n_clusters < 1)) {
    stop('each group must contain at least one cluster')
  }

  cluster.summary <- paste(names(n_clusters), as.numeric(n_clusters), sep='=', collapse=', ')

  if (is.null(strata)) {
    cluster_stratum <- rep('all', length(cluster.ids))
    names(cluster_stratum) <- cluster.ids
    strata.summary <- 'none'
  } else {
    cluster_stratum <- vapply(cluster_list, function(x){
      value <- unique(as.character(x[[strata]]))

      if (length(value) != 1) {
        stop('each cluster must belong to one permutation stratum')
      }

      value
    }, character(1))

    strata.summary <- strata
  }

  stratum.indices <- split(seq_along(cluster.ids), cluster_stratum)

  permutation.counts <- vapply(stratum.indices, function(index){
    n.second <- sum(cluster_group[index] == group.levels[2])
    choose(length(index), n.second)
  }, numeric(1))

  Npermutations <- prod(permutation.counts)

  if (Npermutations < 2) {
    stop('the specified strata permit no exchange of group labels')
  }

  if (Npermutations > max.permutations) {
    stop(
      paste0(
        'the exact permutation design contains ',
        Npermutations,
        ' allocations and exceeds max.permutations=',
        max.permutations
      )
    )
  }

  permutation.parts <- lapply(stratum.indices, function(index){
    n.second <- sum(cluster_group[index] == group.levels[2])

    if (n.second == 0) {
      return(list(integer(0)))
    }

    if (n.second == length(index)) {
      return(list(index))
    }

    combn(index, n.second, simplify=FALSE)
  })

  permutation.grid <- do.call(
    expand.grid,
    lapply(permutation.parts, seq_along)
  )

  permutations <- lapply(seq_len(nrow(permutation.grid)), function(ii){
    unlist(lapply(seq_along(permutation.parts), function(jj){
      permutation.parts[[jj]][[permutation.grid[ii,jj]]]
    }), use.names=FALSE)
  })

  permutation.fun <- function(second.group){
    permuted.group <- rep(group.levels[1], length(cluster.ids))
    permuted.group[second.group] <- group.levels[2]
    names(permuted.group) <- cluster.ids

    dperm <- longdata
    dperm[[group]] <- factor(
      permuted.group[as.character(dperm[[cluster]])],
      levels=group.levels
    )
    contrasts(dperm[[group]]) <- group.contrasts

    fit <- try(rlmer(formula, data=dperm), silent=TRUE)

    if (inherits(fit, 'try-error')) {
      return(c(estimate=NA_real_, statistic=NA_real_))
    }

    se <- sqrt(diag(vcov(fit)))[coef_name]

    if (!is.finite(se) || se <= 0) {
      return(c(estimate=NA_real_, statistic=NA_real_))
    }

    estimate <- contrast.multiplier*unname(fixef(fit)[coef_name])
    statistic <- sign(contrast.multiplier)*unname(fixef(fit)[coef_name]/se)

    c(estimate=estimate, statistic=statistic)
  }

  if (n_cores == 1) {
    permutation.results <- lapply(permutations, permutation.fun)
  } else {
    permutation.results <- parallel::mclapply(permutations, permutation.fun, mc.cores=n_cores)
  }

  permutation.results <- do.call(rbind, permutation.results)

  successful <- is.finite(permutation.results[, 'estimate'])&
    is.finite(permutation.results[, 'statistic'])

  if (!all(successful)) {
    stop(
      paste0(
        sum(!successful),
        ' permutation models failed; the exact permutation p value is unavailable'
      )
    )
  }

  permutation.estimates <- permutation.results[, 'estimate']
  permutation.statistics <- permutation.results[, 'statistic']

  se.observed <- sqrt(diag(vcov(rfit)))[coef_name]
  estimate.observed <- contrast.multiplier*unname(fixef(rfit)[coef_name])
  statistic.observed <- sign(contrast.multiplier)*unname(fixef(rfit)[coef_name]/se.observed)

  tolerance <- sqrt(.Machine$double.eps)*max(1, abs(statistic.observed))

  if (!any(abs(permutation.statistics-statistic.observed) <= tolerance)) {
    stop('the observed allocation is missing from the permutation distribution')
  }

  pval <- mean(
    abs(permutation.statistics) >= abs(statistic.observed)-tolerance
  )

  max.statistic <- max(abs(permutation.statistics))
  maximum.tolerance <- sqrt(.Machine$double.eps)*max(1, max.statistic)

  minimum.p <- mean(
    abs(permutation.statistics) >= max.statistic-maximum.tolerance
  )

  alpha <- 1-CI
  null.interval <- quantile(
    permutation.estimates,
    probs=c(alpha/2, 1-alpha/2)
  )

  summary.permutation <- data.frame(
    parameter=contrast.name,
    target=0,
    null.mean=mean(permutation.estimates),
    se=sd(permutation.estimates),
    null.lower=unname(null.interval[1]),
    null.upper=unname(null.interval[2]),
    minimum.p=minimum.p,
    Npermutations=length(permutation.estimates),
    check.names=FALSE
  )

  rownames(summary.permutation) <- NULL
  model.name <- paste(deparse(formula), collapse=' ')

  summary <- data.frame(
    model=paste0('rlmer exact stratified cluster permutation: ', model.name),
    data=data.name,
    parameter=contrast.name,
    estimate=estimate.observed,
    statistic=statistic.observed,
    'p value'=pval,
    minimum.p=minimum.p,
    clusters=cluster.summary,
    strata=strata.summary,
    method='exact',
    Npermutations=length(permutation.statistics),
    n_cores=n_cores,
    check.names=FALSE
  )

  rownames(summary) <- NULL

  list(
    fit=rfit,
    summary=summary,
    summary.permutation=summary.permutation,
    permutation.estimates=permutation.estimates,
    permutation.statistics=permutation.statistics,
    n.clusters=n_clusters,
    strata=strata,
    method='exact'
  )
}

# bayesian.cluster.analysis() fits a Bayesian hierarchical model and performs prior and posterior predictive analyses.
# The model includes a fixed two-level group effect and a random intercept for the clustering variable, preserving the
# dependence among observations from the same cluster.
# A prior-only model is used to assess whether the specified priors generate plausible outcomes before the observed
# outcome values are used to estimate the posterior distribution.
# The population-level group difference is calculated on the original outcome scale from posterior expected values,
# with the reported parameter defined as the second factor level minus the first.
# The function returns the posterior estimate, credible interval, directional posterior probabilities, latent-scale ICC
# and the probability that one new observation from the second group exceeds one new observation from the first group.
# It also returns prior and posterior group-difference distributions, posterior predictions for new clusters and,
# when requested, a graph comparing the prior and posterior contrast and the two posterior predictive distributions.

bayesian.cluster.graphs <- function(result, pred_xlim=c(0, 1000), plotsave=FALSE,
              svg_path=NULL, filename='bayesian_summary.svg',
              width=9, height=4.5){

  group.levels <- result$group.levels
  prior_diff <- result$prior.difference
  post_diff <- result$posterior.difference
  preds <- result$future.predictions

  p_positive <- mean(post_diff > 0)
  p_negative <- mean(post_diff < 0)
  p_two_sided <- min(1, 2*min(p_positive, p_negative))
  p_future <- mean(preds[,2] > preds[,1])

  diff_xlim <- range(
    unname(quantile(prior_diff, c(0.001, 0.999))),
    unname(quantile(post_diff, c(0.001, 0.999)))
    )

  dev.new(width=width, height=height, noRStudioGD=TRUE)

  par(las=1, mfrow=c(1, 2), mar=c(5.5, 4, 3, 1),
    cex=1, lwd=1, xaxs='i', yaxs='i', tcl=-0.2)

  # A. Prior and posterior group difference

  prior_dens <- density(prior_diff, adjust=1.25, n=2048,
    from=diff_xlim[1], to=diff_xlim[2])

  dens <- density(post_diff, adjust=1.25, n=2048,
    from=diff_xlim[1], to=diff_xlim[2])

  diff_ylim <- c(0, max(c(prior_dens$y, dens$y))*1.2)

  plot(dens, main='Group difference',
    xlab=paste0('Difference (', group.levels[2], ' - ', group.levels[1], ')'),
    xlim=diff_xlim, ylim=diff_ylim,
    lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
    bty='n', axes=FALSE)

  axis(1, lwd=1, cex.axis=0.85)
  axis(2, lwd=1, cex.axis=0.85)
  box(bty='n')

  polygon(c(dens$x[1], dens$x, dens$x[length(dens$x)]),
    c(0, dens$y, 0),
    col=rgb(106/255, 90/255, 205/255, 0.6),
    border='slateblue', lwd=1)

  lines(prior_dens$x, prior_dens$y, col='black', lty=3, lwd=1.5)

  abline(v=0, lty=2, lwd=1)

  mtext('A', side=3, line=1.5, adj=0, font=2, cex=1.5)

  legend('topright', legend=c('Prior', 'Posterior'), col=c('black', 'slateblue'),
    lty=c(3, 1), lwd=c(1.5, 1), bty='n', box.lwd=0, cex=0.85)

  mtext(paste0('Pr(', group.levels[2], '>', group.levels[1], ')=',
    round(p_positive, 3), '; Pr(2-tail)=', round(p_two_sided, 3)),
    side=1, line=4.25, cex=0.75)

  # B. Posterior predictive

  if (any(preds <= 0)) {
    stop('posterior predictions must be positive for the log-density graph')
  }

  log_dens1 <- density(log(preds[,1]), adjust=1.25, n=2048)
  log_dens2 <- density(log(preds[,2]), adjust=1.25, n=2048)

  dens1.x <- exp(log_dens1$x)
  dens2.x <- exp(log_dens2$x)

  dens1.y <- log_dens1$y/dens1.x
  dens2.y <- log_dens2$y/dens2.x

  pred_ylim <- c(0, max(c(dens1.y, dens2.y))*1.2)

  plot(dens1.x, dens1.y, type='l',
    xlab='Predicted outcome for one new observation',
    main='Posterior predictive',
    ylim=pred_ylim, xlim=pred_xlim,
    lwd=1, cex.axis=0.85, cex.lab=0.85, cex.main=0.95,
    bty='n', axes=FALSE)

  axis(1, lwd=1, cex.axis=0.85)
  axis(2, lwd=1, cex.axis=0.85)
  box(bty='n')

  polygon(c(pred_xlim[1], dens1.x, pred_xlim[2]),
    c(0, dens1.y, 0),
    col=rgb(106/255, 90/255, 205/255, 0.6),
    border='slateblue', lwd=1)

  polygon(c(pred_xlim[1], dens2.x, pred_xlim[2]),
    c(0, dens2.y, 0),
    col=rgb(205/255, 92/255, 92/255, 0.6),
    border='indianred', lwd=1)

  mtext('B', side=3, line=1.5, adj=0, font=2, cex=1.5)

  legend('topright',
    fill=c(rgb(106/255, 90/255, 205/255, 0.6),
      rgb(205/255, 92/255, 92/255, 0.6)),
    legend=group.levels, bty='n', box.lwd=0, cex=0.85)

  mtext(paste0('Pr(new ', group.levels[2], ' observation>new ',
    group.levels[1], ' observation)=', round(p_future, 3)),
    side=1, line=4.25, cex=0.75)

  if (plotsave) {
    if (is.null(svg_path)) {
      stop('svg_path must be supplied when plotsave=TRUE')
    }

    save_graph(svg_path=svg_path, filename=filename,
      width=width, height=height, bg='transparent')
  }

  invisible(result)
}


bayesian.cluster.analysis <- function(formula, longdata, cluster, group, prior_spec,
              family=lognormal(), prior.iter=4000, iter=8000,
              chains=4, seed=42, CI=0.95, control=list(adapt_delta=0.99),
              pred_xlim=c(0, 1000), plot.graph=TRUE, plotsave=FALSE,
              svg_path=NULL, filename='bayesian_summary.svg',
              width=9, height=4.5, data.name=''){

  if (!is.factor(longdata[[group]])) {
    stop('group must be a factor')
  }

  longdata[[group]] <- droplevels(longdata[[group]])
  longdata[[cluster]] <- droplevels(factor(longdata[[cluster]]))

  group.levels <- levels(longdata[[group]])

  if (length(group.levels) != 2) {
    stop('group must contain exactly two levels')
  }

  group.contrasts <- contrasts(longdata[[group]])

  if (is.null(group.contrasts) || ncol(group.contrasts) != 1) {
    stop('group must have one two-level contrast')
  }

  if (!nzchar(data.name)) {
    data.name <- deparse(substitute(longdata))
  }

  alpha <- 1-CI

  # Prior predictive model

  prior_model <- brm(bf(formula), data=longdata, family=family,
    prior=prior_spec, sample_prior='only',
    iter=prior.iter, chains=chains, seed=seed)

  prior_predictions <- posterior_predict(prior_model)

  prior.predictive.summary <- data.frame(
    minimum=apply(prior_predictions, 1, min),
    median=apply(prior_predictions, 1, median),
    maximum=apply(prior_predictions, 1, max)
    )

  prior.predictive.quantiles <- rbind(
    minimum=quantile(prior.predictive.summary$minimum, c(alpha/2, 0.5, 1-alpha/2)),
    median=quantile(prior.predictive.summary$median, c(alpha/2, 0.5, 1-alpha/2)),
    maximum=quantile(prior.predictive.summary$maximum, c(alpha/2, 0.5, 1-alpha/2))
    )

  prior.proportion.above.limit <- mean(prior.predictive.summary$maximum > pred_xlim[2])
  prior.upper <- unname(quantile(prior.predictive.summary$maximum, 1-alpha/2))

  # Fit the Bayesian model

  fitted_bayesian_model <- brm(bf(formula), data=longdata, family=family,
    prior=prior_spec, iter=iter, chains=chains, seed=seed,
    control=control)

  fixed.effects <- fixef(fitted_bayesian_model)
  coef_name <- setdiff(rownames(fixed.effects), 'Intercept')

  if (length(coef_name) != 1) {
    stop('the model must contain one non-intercept fixed-effect coefficient')
  }

  coef_name_b <- paste0('b_', coef_name)

  coefficient.summary <- posterior_summary(fitted_bayesian_model, variable=coef_name_b)
  draws <- as_draws_df(fitted_bayesian_model)

  # Population-level group difference

  newdata.mean <- longdata[rep(1, 2),,drop=FALSE]

  newdata.mean[[group]] <- factor(group.levels, levels=group.levels)
  newdata.mean[[cluster]] <- factor(rep(levels(longdata[[cluster]])[1], 2),
    levels=levels(longdata[[cluster]]))

  contrasts(newdata.mean[[group]]) <- group.contrasts

  epred <- posterior_epred(fitted_bayesian_model, newdata=newdata.mean, re_formula=NA)

  # Difference is the second factor level minus the first factor level.

  post_diff <- epred[,2]-epred[,1]

  post_ci <- quantile(post_diff, c(alpha/2, 1-alpha/2))

  p_positive <- mean(post_diff > 0)
  p_negative <- mean(post_diff < 0)
  p_two_sided <- min(1, 2*min(p_positive, p_negative))

  # Prior distribution of the same group difference

  prior_epred <- posterior_epred(prior_model, newdata=newdata.mean, re_formula=NA)

  prior_diff <- prior_epred[,2]-prior_epred[,1]

  prior_ci <- quantile(prior_diff, c(alpha/2, 1-alpha/2))

  # Latent-scale ICC

  cluster.sd.name <- grep('^sd_.*__Intercept$', names(draws), value=TRUE)

  if (length(cluster.sd.name) != 1 || !'sigma' %in% names(draws)) {
    stop('the model must contain one random-intercept standard deviation and sigma')
  }

  icc <- draws[[cluster.sd.name]]^2/(draws[[cluster.sd.name]]^2+draws$sigma^2)

  icc_ci <- quantile(icc, c(alpha/2, 1-alpha/2))

  # Predictions for new observations from new clusters

  new.cluster.levels <- paste0('new_', seq_along(group.levels))

  newdata.future <- longdata[rep(1, 2),,drop=FALSE]

  newdata.future[[group]] <- factor(group.levels, levels=group.levels)
  newdata.future[[cluster]] <- factor(new.cluster.levels,
    levels=c(levels(longdata[[cluster]]), new.cluster.levels))

  contrasts(newdata.future[[group]]) <- group.contrasts

  set.seed(seed)

  preds <- posterior_predict(fitted_bayesian_model, newdata=newdata.future,
    re_formula=NULL, allow_new_levels=TRUE, sample_new_levels='gaussian')

  p_future <- mean(preds[,2] > preds[,1])

  # Summary

  model.name <- paste0(
    'brm: ',
    paste(deparse(formula), collapse=' '),
    ', family = ',
    stats::family(fitted_bayesian_model)$family,
    '()'
    )

  parameter.name <- paste(group.levels[2], '-', group.levels[1])

  bayes_summary <- data.frame(
    model=model.name,
    data=data.name,
    parameter=parameter.name,
    model_parameter=coef_name_b,
    estimate=mean(post_diff),
    'CI (2.5)'=unname(post_ci[1]),
    'CI (97.5)'=unname(post_ci[2]),
    p_negative,
    p_positive,
    p_two_sided,
    p_future,
    check.names=FALSE
    )

  names(bayes_summary)[8:11] <- c(
    paste0('Pr(', group.levels[1], '>', group.levels[2], ')'),
    paste0('Pr(', group.levels[2], '>', group.levels[1], ')'),
    'Pr(2-tail)',
    paste0('Pr(new ', group.levels[2], ' obs>new ',group.levels[1], ' obs)')
    )

  rownames(bayes_summary) <- NULL

  icc.summary <- data.frame(
    median=median(icc),
    'CI (2.5%)'=unname(icc_ci[1]),
    'CI (97.5%)'=unname(icc_ci[2]),
    check.names=FALSE
    )

  prior.difference.summary <- data.frame(
    mean=mean(prior_diff),
    median=median(prior_diff),
    'CI (2.5%)'=unname(prior_ci[1]),
    'CI (97.5%)'=unname(prior_ci[2]),
    check.names=FALSE
    )

  result <- list(
    prior.model=prior_model,
    fitted.model=fitted_bayesian_model,
    summary=bayes_summary,
    coefficient.summary=coefficient.summary,
    prior.predictive.summary=prior.predictive.summary,
    prior.predictive.quantiles=prior.predictive.quantiles,
    prior.proportion.above.limit=prior.proportion.above.limit,
    prior.upper=prior.upper,
    prior.difference.summary=prior.difference.summary,
    ICC=icc.summary,
    group.levels=group.levels,
    prior.difference=prior_diff,
    posterior.difference=post_diff,
    future.predictions=preds
    )

  if (plot.graph) {
    bayesian.cluster.graphs(result=result, pred_xlim=pred_xlim,
      plotsave=plotsave, svg_path=svg_path, filename=filename,
      width=width, height=height)
  }

  result
}

#' p.adjust.ff
p.adjust.ff <- function (p, method = c("Sidak-Holm", "Bonf", "Sidak-Bonf", "Holm-Bonf", "Hochberg", "Hommel", "BH", "BY", "BH.a", "BH.ts"), alpha = 0.05) {
    if (is.null(p)) p <- numeric(0L)
    if (!is.numeric(p) || anyNA(p) || any(!is.finite(p)) ||
            any(p < 0 | p > 1))
        stop('p must contain only finite numbers between zero and one', call.=FALSE)
    if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
            alpha <= 0 || alpha >= 1)
        stop('alpha must be one number between zero and one', call.=FALSE)
    # Only use match.arg if method is not NULL
    if (!is.null(method)) {
        method <- match.arg(method)
        if (!length(p)) return(data.frame(
            p=numeric(0L), sig.level=numeric(0L), p.adj=numeric(0L),
            sig=logical(0L)))
        if (length(p)==1L) return(data.frame(
            p=p, sig.level=alpha, p.adj=p, sig=p<=alpha))
    }    
    # Proceed only if method is not NULL
    if (!is.null(method)) {
        if (method == "BH") {
            return(bh.ff(p=p, alpha=alpha))
        } else if (method == "BY") {
            return(by.ff(p=p, alpha=alpha))
        } else if (method == "BH.a") {
            return(bh.a.ff(p=p, alpha=alpha))
        } else if (method == "BH.ts") {
            return(bh.ts.ff(p=p, alpha=alpha))
        } else if (method == "Bonf") {
            return(bonferroni.ff(p=p, alpha=alpha))
        } else if (method == "Sidak-Bonf") {
            return(sidak.bonferroni.ff(p=p, alpha=alpha))
        } else if (method == "Holm-Bonf") {
            return(holm.bonferroni.ff(p=p, alpha=alpha))
        } else if (method == "Sidak-Holm") {
            return(sidak.holm.ff(p=p, alpha=alpha))
        } else if (method == "Hochberg") {
            return(hochberg.ff(p=p, alpha=alpha))
        } else if (method == "Hommel") {
            return(hommel.ff(p=p, alpha=alpha))
        }
    } else {
        # handling if method is NULL
        alpha.levels <- rep(alpha, length(p))
        sig <- p <= alpha.levels
        out <- as.data.frame(cbind(p, alpha = alpha.levels))
        out$sig <- sig
        return(out)
    }
}

#' bh.ff
bh.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    o <- order(p)
    ro <- order(o)
    p.ord <-  p[o]
    alpha.levels <- i/n * alpha 
    sig <- p.ord < alpha.levels
    if ( length(which(sig)) != 0 ){
        ind <- max( which(sig) )
        sig <- i <= ind
    }
    p.adj <- pmin( 1, rev( cummin(rev(alpha/alpha.levels * p.ord ) ) ) )
    p.adj <- p.adj[ro]  
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]
    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj)) 
    out$sig=sig
    out
}

#' bh.a.ff
bh.a.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    o <- order(p)
    ro <- order(o)
    p.ord <- p[o]

    p.adj.initial <- pmin(1, rev(cummin(rev(n/i*p.ord))))
    n.rejected <- sum(p.adj.initial <= alpha)

    if (n.rejected==0 || n.rejected==n) {
        alpha.levels <- i/n*alpha
        p.adj <- p.adj.initial
    } else {
        h0.sequence <- (n+1-i)/(1-p.ord)
        increase <- which(diff(h0.sequence) > 0)

        if (length(increase)==0) {
            h0.a <- n
        } else {
            stop.index <- increase[1]+1
            h0.a <- min(ceiling(h0.sequence[stop.index]), n)
        }

        alpha.levels <- pmin(1, i/h0.a*alpha)
        p.adj <- pmin(1, rev(cummin(rev(h0.a/i*p.ord))))
    }

    p.adj <- p.adj[ro]
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]

    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj))
    out$sig <- sig
    out
}

#' bh.ts.ff
bh.ts.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    o <- order(p)
    ro <- order(o)
    p.ord <- p[o]

    alpha.first <- alpha/(1+alpha)
    p.adj.first <- pmin(1, rev(cummin(rev(n/i*p.ord))))
    n.rejected <- sum(p.adj.first <= alpha.first)

    if (n.rejected==0 || n.rejected==n) {
        alpha.levels <- i/n*alpha.first
        p.adj <- pmin(1, alpha/alpha.first*p.adj.first)
    } else {
        h0.TSBH <- n-n.rejected
        alpha.levels <- pmin(1, i/h0.TSBH*alpha.first)
        p.adj <- pmin(
            1,
            rev(cummin(rev(alpha/alpha.first*h0.TSBH/i*p.ord)))
        )
    }

    p.adj <- p.adj[ro]
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]

    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj))
    out$sig <- sig
    out
}

#' by.ff
by.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    q <- sum(1/i)
    o <- order(p)
    ro <- order(o)
    p.ord <-  p[o]
    alpha.levels <- i/n * alpha/q   #' Theorem 1.3
    sig <- p.ord < alpha.levels 
    if ( length(which(sig)) != 0 ){
        ind <- max( which(sig) )
        sig <- i <= ind
    }
    p.adj <- pmin( 1, rev( cummin(rev(alpha/alpha.levels * p.ord ) ) ) )
    p.adj <- p.adj[ro]  
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]
    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj)) 
    out$sig=sig
    out
}

#' bonferroni.ff
bonferroni.ff <- function(p, alpha){
    n <- length(p)
    alpha.levels <- rep(alpha/n, n)
    sig <- p <= alpha.levels
    p.adj <- pmin(1, alpha/alpha.levels * p)
    out <- as.data.frame( cbind( p, sig.level=alpha.levels,p.adj=p.adj) )
    out$sig=sig
    out
}

#' sidak.bonferroni.ff
sidak.bonferroni.ff <- function(p, alpha){
    n <- length(p)
    alpha.levels <- 1 - (1-alpha)^(1/n)
    alpha.levels <- rep(alpha.levels, n)
    sig <- p <= alpha.levels
    p.adj <- pmin(1, -expm1(n * log1p(-p)))
    out <- as.data.frame( cbind( p, sig.level=alpha.levels,p.adj=p.adj) )
    out$sig=sig
    out
}


#' holm.bonferroni.ff
holm.bonferroni.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    alpha.levels <- alpha/(n-i+1)
    o <- order(p)
    ro <- order(o)
    p.ord <-  p[o]
    p.adj <- pmin(1, cummax(alpha/alpha.levels * p.ord))
    sig <- p.adj <= alpha
    p.adj <- p.adj[ro]  
    sig <- sig[ro]
    alpha.levels <- alpha.levels[ro]
    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj)) 
    out$sig=sig
    out
}

#' sidak.holm.ff
sidak.holm.ff <- function(p, alpha){
    n <- length(p)
    remaining <- n:1
    alpha.levels <- -expm1(log1p(-alpha)/remaining)
    o <- order(p)
    ro <- order(o)
    p.ord <-  p[o]
    p.adj <- pmin(1, cummax(-expm1(remaining * log1p(-p.ord))))
    sig <- p.adj <= alpha
    p.adj <- p.adj[ro]
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]
    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj)) 
    out$sig=sig
    out
}

#' hochberg.ff
hochberg.ff <- function(p, alpha){
    n <- length(p)
    i <- n:1
    alpha.levels <- alpha/(n-i+1)
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    p.ord <-  p[o]
    sig <- p.ord < alpha.levels
    if ( length(which(sig)) != 0 ){
        ind <- min( which(sig) )
        sig <- ind <= 1:n
    }
    p.adj <- rev( pmin( 1, rev( cummin(alpha/alpha.levels * p.ord ) ) ) )
    p.adj <- p.adj[ro]  
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]
    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj)) 
    out$sig=sig
    out
}

#' hommel.ff
hommel.ff <- function(p, alpha){
    n <- length(p)
    i <- seq_len(n)
    o <- order(p)
    ro <- order(o)
    p.ord <- p[o]

    j.finder <- sapply(seq_len(n), function(j){
        k <- seq_len(j)
        sum(p.ord[n-j+k] > k*alpha/j) == j
    })

    if (sum(j.finder)==0) {
        alpha.levels <- rep(alpha, n)
    } else {
        j <- max(which(j.finder))
        alpha.levels <- rep(alpha/j, n)
    }

    p.adj <- stats::p.adjust(p, method="hommel")
    sig <- p.adj <= alpha
    alpha.levels <- alpha.levels[ro]

    out <- as.data.frame(cbind(p, sig.level=alpha.levels, p.adj=p.adj))
    out$sig <- sig
    out
}

# summary for nonparametric Brunner-Munzel nonparametric Behrens-Fisher test
nparcomp.summary <- function(object, paired.method=c('BM', 'PERM', 'both')){

  model.data <- model.frame(formula=object$input$formula, data=object$input$data, na.action=na.omit)

  response.name <- names(model.data)[1]
  group.name <- names(model.data)[2]
  group.levels <- as.character(object$Info$Sample)
  group.n <- object$Info$Size

  CI.probabilities <- switch(object$input$alternative[1],
    two.sided=c((1-object$input$conf.level)/2, 1-(1-object$input$conf.level)/2),
    less=c(0, object$input$conf.level),
    greater=c(1-object$input$conf.level, 1)
    )

  lower.name <- paste0('CI (', format(CI.probabilities[1], trim=TRUE, scientific=FALSE), ')')
  upper.name <- paste0('CI (', format(CI.probabilities[2], trim=TRUE, scientific=FALSE), ')')

  if (inherits(object, 'nparttest')){

    result.df <- NA_real_

    if (!is.null(object$AsyMethod) && grepl(' with .* DF$', object$AsyMethod)){
      result.df <- as.numeric(sub(' DF$', '', sub('.* with ', '', object$AsyMethod)))
    }

    stats_summary <- data.frame(
      parameter=response.name,
      comparison=paste('within', group.name, '(unpaired)'),
      contrast=paste(group.levels[2], 'vs', group.levels[1]),
      n=paste(group.n[2], 'vs', group.n[1]),
      test='Brunner-Munzel nonparametric Behrens-Fisher test',
      alternative=object$input$alternative[1],
      'test stat'='T',
      df=result.df,
      stat=object$Analysis$T,
      estimand=object$Analysis$Effect,
      'relative effect'=object$Analysis$Estimator,
      lower=object$Analysis$Lower,
      upper=object$Analysis$Upper,
      'p value'=object$Analysis$p.Value,
      check.names=FALSE
      )

  } else if (inherits(object, 'nparttestpaired')){

    paired.method <- match.arg(paired.method)

    if (paired.method=='both'){
      result.rows <- rownames(object$Analysis)
    } else {
      result.rows <- paired.method
    }

    result.analysis <- object$Analysis[result.rows,,drop=FALSE]

    test.names <- c(
      BM='Paired Brunner-Munzel test',
      PERM='Paired studentized permutation test'
      )

    stats_summary <- data.frame(
      parameter=rep(response.name, nrow(result.analysis)),
      comparison=rep(paste('within', group.name, '(paired)'), nrow(result.analysis)),
      contrast=rep(paste(group.levels[2], 'vs', group.levels[1]), nrow(result.analysis)),
      n=rep(group.n[1], nrow(result.analysis)),
      test=unname(test.names[rownames(result.analysis)]),
      alternative=rep(object$input$alternative[1], nrow(result.analysis)),
      'test stat'=rep('T', nrow(result.analysis)),
      df=ifelse(rownames(result.analysis)=='BM', group.n[1]-1, NA_real_),
      stat=result.analysis[, 'T'],
      estimand=rep(paste0('p(', group.levels[1], ',', group.levels[2], ')'), nrow(result.analysis)),
      'relative effect'=result.analysis[, 'p.hat'],
      lower=result.analysis[, 'Lower'],
      upper=result.analysis[, 'Upper'],
      'p value'=result.analysis[, 'p.value'],
      check.names=FALSE
      )

  } else {
    stop('object must be returned by npar.t.test() or npar.t.test.paired()')
  }

  names(stats_summary)[names(stats_summary)=='lower'] <- lower.name
  names(stats_summary)[names(stats_summary)=='upper'] <- upper.name

  rownames(stats_summary) <- NULL

  stats_summary
}

# add adjusted p values to a combined nparcomp summary
nparcomp.adjust <- function(x){

  if (!is.data.frame(x)){
    stop('x must be a data frame')
  }

  if (!'p value' %in% names(x)){
    stop("x must contain a 'p value' column")
  }

  npar.adjusted <- p.adjust.ff(x[['p value']], method='Sidak-Holm', alpha=0.05)
  x[['p adjusted']] <- npar.adjusted[,'p.adj']

  x
}

nparcomp2w <- function(formula, longdata, subject, alternative='two.sided', CI=0.95, alpha=0.05){

  formula.variables <- all.vars(formula)
  response.name <- formula.variables[1]
  predictor.names <- formula.variables[-1]

  if (length(predictor.names)!=2){
    stop('formula must contain exactly two predictors')
  }

  required.names <- c(response.name, predictor.names, subject)

  if (!all(required.names %in% names(longdata))){
    stop('formula or subject variables are missing from longdata')
  }

  analysis.data <- longdata[,required.names,drop=FALSE]
  analysis.data <- analysis.data[complete.cases(analysis.data),,drop=FALSE]

  if (!nrow(analysis.data)){
    stop('no complete observations are available')
  }

  analysis.data[[subject]] <- droplevels(factor(analysis.data[[subject]]))

  for (predictor in predictor.names){
    analysis.data[[predictor]] <- droplevels(factor(analysis.data[[predictor]]))

    if (nlevels(analysis.data[[predictor]])!=2){
      stop('each predictor must contain exactly two levels')
    }
  }

  subject.rows <- split(seq_len(nrow(analysis.data)), analysis.data[[subject]])

  varies.within.subject <- vapply(predictor.names, function(predictor){
    any(vapply(subject.rows, function(i){
      length(unique(analysis.data[[predictor]][i]))>1
    }, logical(1)))
  }, logical(1))

  if (sum(varies.within.subject)!=1){
    stop('one predictor must vary within subject and one must vary between subjects')
  }

  within.name <- predictor.names[varies.within.subject]
  between.name <- predictor.names[!varies.within.subject]

  within.levels <- levels(analysis.data[[within.name]])
  between.levels <- levels(analysis.data[[between.name]])

  between.per.subject <- vapply(subject.rows, function(i){
    length(unique(analysis.data[[between.name]][i]))
  }, integer(1))

  if (any(between.per.subject!=1)){
    stop('the between-subject factor must be constant within subject')
  }

  result.summary <- vector('list', length(within.levels)+length(between.levels))
  result.index <- 1

  for (within.level in within.levels){

    test.data <- analysis.data[analysis.data[[within.name]]==within.level,,drop=FALSE]
    test.data <- droplevels(test.data)

    if (anyDuplicated(test.data[[subject]])){
      stop('independent comparisons require one observation per subject and factor level')
    }

    result <- nparcomp::npar.t.test(formula=reformulate(between.name, response=response.name), data=test.data,
      method='t.app', alternative=alternative, conf.level=CI, info=FALSE, rounds=Inf)

    result.summary[[result.index]] <- transform(nparcomp.summary(result),
      parameter=paste(within.level, response.name))

    result.index <- result.index+1
  }

  for (between.level in between.levels){

    test.data <- analysis.data[analysis.data[[between.name]]==between.level,,drop=FALSE]
    test.data <- droplevels(test.data)
    test.data <- test.data[order(test.data[[subject]], test.data[[within.name]]),,drop=FALSE]

    pair.table <- table(test.data[[subject]], test.data[[within.name]])

    if (any(pair.table!=1)){
      stop('paired comparisons require one complete observation per subject and factor level')
    }

    result <- nparcomp::npar.t.test.paired(formula=reformulate(within.name, response=response.name), data=test.data,
      alternative=alternative, conf.level=CI, info=FALSE, plot.simci=FALSE, rounds=Inf)

    result.summary[[result.index]] <- transform(nparcomp.summary(result, paired.method='BM'),
      parameter=paste(between.level, response.name))

    result.index <- result.index+1
  }

  stats_summary <- do.call(rbind, result.summary)
  rownames(stats_summary) <- NULL

  npar.adjusted <- p.adjust.ff(stats_summary$'p value', method='Sidak-Holm', alpha=alpha)
  stats_summary$'p adjusted' <- npar.adjusted[,'p.adj']

  stats_summary
}

MCnpartest <- function(formula, longdata, subject=NULL, contr=NULL, conf.level=0.95,
  alternative=c('two.sided', 'less', 'greater'),
  independent.method=c('t.app', 'logit', 'probit', 'normal', 'permu'),
  paired.method=c('BM', 'PERM', 'both'), rounds=Inf, plot.simci=FALSE,
  info=FALSE, nperm=10000, p.adjust.method='Sidak-Holm', alpha=0.05){

  alternative <- match.arg(alternative)
  independent.method <- match.arg(independent.method)
  paired.method <- match.arg(paired.method)

  formula.variables <- all.vars(formula)
  response.name <- formula.variables[1]
  predictor.names <- formula.variables[-1]

  if (!(length(predictor.names) %in% c(1,2))){
    stop('formula must contain one or two predictors')
  }

  if (length(predictor.names)==2 && is.null(subject)){
    stop('subject must be specified for a two-factor design')
  }

  required.names <- c(response.name, predictor.names, subject)

  if (!all(required.names %in% names(longdata))){
    stop('formula or subject variables are missing from longdata')
  }

  analysis.data <- longdata[,required.names,drop=FALSE]
  analysis.data <- analysis.data[complete.cases(analysis.data),,drop=FALSE]

  if (!nrow(analysis.data)){
    stop('no complete observations are available')
  }

  for (predictor in predictor.names){
    analysis.data[[predictor]] <- droplevels(factor(analysis.data[[predictor]]))

    if (nlevels(analysis.data[[predictor]])<2){
      stop('each predictor must contain at least two levels')
    }
  }

  if (!is.null(subject)){
    analysis.data[[subject]] <- droplevels(factor(analysis.data[[subject]]))
  }

  if (length(predictor.names)==1){
    cell <- analysis.data[[predictor.names]]
  } else {
    cell <- do.call(interaction, c(analysis.data[predictor.names],
      list(sep=':', drop=TRUE, lex.order=FALSE)))
  }

  cell <- droplevels(factor(cell))
  cell.levels <- levels(cell)
  cell.grid <- analysis.data[match(cell.levels, as.character(cell)),predictor.names,drop=FALSE]
  rownames(cell.grid) <- cell.levels

  if (is.null(contr)){

    result.contrasts <- list()

    if (length(predictor.names)==1){
      predictor.levels <- levels(analysis.data[[predictor.names]])

      for (i in seq_len(length(predictor.levels)-1)){
        contrast.row <- rep(0, length(cell.levels))
        contrast.row[match(predictor.levels[i], cell.levels)] <- -1
        contrast.row[match(predictor.levels[i+1], cell.levels)] <- 1
        result.contrasts[[length(result.contrasts)+1]] <- contrast.row
      }

    } else {

      subject.rows <- split(seq_len(nrow(analysis.data)), analysis.data[[subject]])

      varies.within.subject <- vapply(predictor.names, function(predictor){
        any(vapply(subject.rows, function(i){
          length(unique(analysis.data[[predictor]][i]))>1
        }, logical(1)))
      }, logical(1))

      if (sum(varies.within.subject)!=1){
        stop('one predictor must vary within subject and one must vary between subjects')
      }

      predictor.order <- c(predictor.names[!varies.within.subject],
        predictor.names[varies.within.subject])

      for (target.name in predictor.order){

        conditioning.name <- setdiff(predictor.names, target.name)
        conditioning.levels <- levels(analysis.data[[conditioning.name]])
        target.levels <- levels(analysis.data[[target.name]])

        for (conditioning.level in conditioning.levels){
          for (i in seq_len(length(target.levels)-1)){

            negative.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i]
              ]

            positive.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i+1]
              ]

            if (length(negative.cell)==1 && length(positive.cell)==1){
              contrast.row <- rep(0, length(cell.levels))
              contrast.row[match(negative.cell, cell.levels)] <- -1
              contrast.row[match(positive.cell, cell.levels)] <- 1
              result.contrasts[[length(result.contrasts)+1]] <- contrast.row
            }
          }
        }
      }
    }

    contr <- do.call(rbind, result.contrasts)
    colnames(contr) <- cell.levels

  } else {

    if (is.numeric(contr) && is.null(dim(contr))){
      contr <- matrix(contr, nrow=1)
    }

    if (!is.matrix(contr) || !is.numeric(contr) || !nrow(contr) ||
      anyNA(contr) || any(!is.finite(contr)) || ncol(contr)!=length(cell.levels)){
      stop('contr must be a finite numeric matrix with one column per factor-level combination')
    }

    if (!is.null(colnames(contr))){
      if (!setequal(colnames(contr), cell.levels)){
        stop(paste('contr column names must match:', paste(cell.levels, collapse=', ')))
      }

      contr <- contr[,cell.levels,drop=FALSE]
    } else {
      colnames(contr) <- cell.levels
    }
  }

  pairwise.contrast <- apply(contr, 1, function(x){
    nonzero <- sort(as.numeric(x[x!=0]))
    length(nonzero)==2 && all(nonzero==c(-1,1))
  })

  if (any(!pairwise.contrast)){
    stop('each row of contr must contain one -1, one 1 and zeros elsewhere')
  }

  CI.probabilities <- switch(alternative,
    two.sided=c((1-conf.level)/2, 1-(1-conf.level)/2),
    less=c(0, conf.level),
    greater=c(1-conf.level, 1)
    )

  lower.name <- paste0('CI (', format(CI.probabilities[1], trim=TRUE, scientific=FALSE), ')')
  upper.name <- paste0('CI (', format(CI.probabilities[2], trim=TRUE, scientific=FALSE), ')')
  results <- list()

  for (i in seq_len(nrow(contr))){

    negative.cell <- cell.levels[contr[i,]==-1]
    positive.cell <- cell.levels[contr[i,]==1]
    keep <- cell %in% c(negative.cell, positive.cell)

    test.data <- analysis.data[keep,,drop=FALSE]

    if (!is.null(subject)){
      test.data[[subject]] <- droplevels(test.data[[subject]])
    }

    test.data$.response <- test.data[[response.name]]
    test.data$.cell <- factor(as.character(cell[keep]), levels=c(negative.cell, positive.cell))

    negative.data <- test.data[test.data$.cell==negative.cell,,drop=FALSE]
    positive.data <- test.data[test.data$.cell==positive.cell,,drop=FALSE]

    paired <- FALSE

    if (!is.null(subject)){
      negative.subjects <- unique(as.character(negative.data[[subject]]))
      positive.subjects <- unique(as.character(positive.data[[subject]]))
      common.subjects <- intersect(negative.subjects, positive.subjects)

      if (setequal(negative.subjects, positive.subjects)){
        paired <- TRUE
      } else if (length(common.subjects)){
        stop('a contrast cannot contain a mixture of paired and unpaired subjects')
      }
    }

    if (paired){
      pair.table <- table(test.data[[subject]], test.data$.cell)

      if (any(pair.table!=1)){
        stop('paired comparisons require one complete observation per subject and contrast level')
      }

      test.data <- test.data[order(test.data[[subject]], test.data$.cell),,drop=FALSE]

      n.pairs <- nrow(positive.data)

      if (paired.method!='BM' && n.pairs>13){
        if (!is.numeric(nperm) || length(nperm)!=1 || is.na(nperm) ||
          !is.finite(nperm) || nperm!=10000){
          stop("paired.method='PERM' or 'both' requires nperm=10000 when there are more than 13 pairs")
        }

        paired.nperm <- nperm
      } else {
        paired.nperm <- 10000
      }

      result <- nparcomp::npar.t.test.paired(.response ~ .cell, data=test.data,
        conf.level=conf.level, alternative=alternative, nperm=paired.nperm, rounds=rounds,
        info=info, plot.simci=plot.simci)

      result.rows <- if (paired.method=='both') rownames(result$Analysis) else paired.method
      result.analysis <- result$Analysis[result.rows,,drop=FALSE]

      test.names <- c(
        BM='Paired Brunner-Munzel test',
        PERM='Paired studentized permutation test'
        )

      result.summary <- data.frame(
        test=unname(test.names[rownames(result.analysis)]),
        'test stat'=rep('T', nrow(result.analysis)),
        df=ifelse(rownames(result.analysis)=='BM', nrow(positive.data)-1, NA_real_),
        stat=result.analysis[, 'T'],
        estimand=rep(paste0('p(', negative.cell, ',', positive.cell, ')'), nrow(result.analysis)),
        'relative effect'=result.analysis[, 'p.hat'],
        lower=result.analysis[, 'Lower'],
        upper=result.analysis[, 'Upper'],
        'p value'=result.analysis[, 'p.value'],
        adjustment.method=rownames(result.analysis),
        check.names=FALSE
        )

      n.output <- as.character(nrow(positive.data))

    } else {

      if (!is.null(subject) &&
        (anyDuplicated(negative.data[[subject]]) || anyDuplicated(positive.data[[subject]]))){
        stop('unpaired comparisons require one observation per subject and contrast level')
      }

      result <- nparcomp::npar.t.test(.response ~ .cell, data=test.data,
        conf.level=conf.level, alternative=alternative, rounds=rounds,
        method=independent.method, plot.simci=plot.simci, info=info, nperm=nperm)

      if (independent.method=='permu'){

        method.names <- c(
          id='Studentized permutation test',
          logit='Studentized permutation test (logit)',
          probit='Studentized permutation test (probit)'
          )

        result.summary <- data.frame(
          test=unname(method.names[rownames(result$Analysis)]),
          'test stat'=rep('T', nrow(result$Analysis)),
          df=rep(NA_real_, nrow(result$Analysis)),
          stat=result$Analysis[, 'Statistic'],
          estimand=rep(paste0('p(', negative.cell, ',', positive.cell, ')'), nrow(result$Analysis)),
          'relative effect'=result$Analysis[, 'Estimator'],
          lower=result$Analysis[, 'Lower'],
          upper=result$Analysis[, 'Upper'],
          'p value'=result$Analysis[, 'p.value'],
          adjustment.method=rownames(result$Analysis),
          check.names=FALSE
          )

      } else {

        method.names <- c(
          't.app'='Brunner-Munzel nonparametric Behrens-Fisher test',
          logit='Nonparametric Behrens-Fisher logit approximation',
          probit='Nonparametric Behrens-Fisher probit approximation',
          normal='Nonparametric Behrens-Fisher normal approximation'
          )

        result.df <- NA_real_

        if (independent.method=='t.app' && grepl(' with .* DF$', result$AsyMethod)){
          result.df <- as.numeric(sub(' DF$', '', sub('.* with ', '', result$AsyMethod)))
        }

        result.summary <- data.frame(
          test=unname(method.names[independent.method]),
          'test stat'='T',
          df=result.df,
          stat=result$Analysis[, 'T'],
          estimand=result$Analysis[, 'Effect'],
          'relative effect'=result$Analysis[, 'Estimator'],
          lower=result$Analysis[, 'Lower'],
          upper=result$Analysis[, 'Upper'],
          'p value'=result$Analysis[, 'p.Value'],
          adjustment.method=independent.method,
          check.names=FALSE
          )
      }

      n.output <- paste(nrow(positive.data), 'vs', nrow(negative.data))
    }

    differing <- predictor.names[
      vapply(predictor.names, function(x){
        as.character(cell.grid[negative.cell,x])!=as.character(cell.grid[positive.cell,x])
      }, logical(1))
      ]

    same <- setdiff(predictor.names, differing)

    if (length(predictor.names)==1){
      parameter <- response.name
      comparison <- paste('within', predictor.names, if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(as.character(cell.grid[positive.cell,differing]), 'vs',
        as.character(cell.grid[negative.cell,differing]))
    } else if (length(differing)==1){
      parameter <- response.name
      comparison <- paste('within', same, as.character(cell.grid[positive.cell,same]),
        if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(as.character(cell.grid[positive.cell,differing]), 'vs',
        as.character(cell.grid[negative.cell,differing]))
    } else {
      parameter <- response.name
      comparison <- paste('between cells', if (paired) '(paired)' else '(unpaired)')
      contrast.name <- paste(positive.cell, 'vs', negative.cell)
    }

    result.summary <- data.frame(
      parameter=rep(parameter, nrow(result.summary)),
      comparison=rep(comparison, nrow(result.summary)),
      contrast=rep(contrast.name, nrow(result.summary)),
      n=rep(n.output, nrow(result.summary)),
      test=result.summary$test,
      alternative=rep(alternative, nrow(result.summary)),
      'test stat'=result.summary$'test stat',
      df=result.summary$df,
      stat=result.summary$stat,
      estimand=result.summary$estimand,
      'relative effect'=result.summary$'relative effect',
      lower=result.summary$lower,
      upper=result.summary$upper,
      'p value'=result.summary$'p value',
      family=rep(if (paired) 'paired' else 'unpaired', nrow(result.summary)),
      adjustment.method=result.summary$adjustment.method,
      check.names=FALSE
      )

    names(result.summary)[names(result.summary)=='lower'] <- lower.name
    names(result.summary)[names(result.summary)=='upper'] <- upper.name

    results[[length(results)+1]] <- result.summary
  }

  stats_summary <- do.call(rbind, results)
  rownames(stats_summary) <- NULL

  stats_summary$'p adjusted' <- NA_real_

  adjustment.family <- interaction(stats_summary$family, stats_summary$adjustment.method,
    drop=TRUE)

  for (family in levels(adjustment.family)){
    index <- which(adjustment.family==family)
    adjusted <- p.adjust.ff(stats_summary$'p value'[index], method=p.adjust.method, alpha=alpha)
    stats_summary$'p adjusted'[index] <- adjusted[,'p.adj']
  }

  stats_summary$family <- NULL
  stats_summary$adjustment.method <- NULL

  stats_summary
}

LMsummary <- function(fit, between, data.name, p.adjust.method='Sidak-Holm',
  alpha=0.05, CI=95){

  CI.level <- CI/100
  analysis.data <- model.frame(fit)
  response <- all.vars(formula(fit))[1]
  results <- NULL
  comparison <- NULL
  n <- NULL

  for (i in seq_along(between)){

    other <- setdiff(between, between[i])
    specs <- if (length(other)==0){
      as.formula(paste('~', between[i]))
    } else {
      as.formula(paste('~', between[i], '|', paste(other, collapse='*')))
    }

    emm <- emmeans::emmeans(fit, specs=specs)
    comparisons <- emmeans::contrast(emm, method='revpairwise', adjust='none')
    result <- as.data.frame(summary(comparisons, infer=c(TRUE,TRUE),
      by=NULL, adjust='none', level=CI.level))

    contrast.levels <- strsplit(as.character(result$contrast), ' - ', fixed=TRUE)

    result.n <- vapply(seq_len(nrow(result)), function(j){

      positive <- as.character(analysis.data[[between[i]]])==contrast.levels[[j]][1]
      negative <- as.character(analysis.data[[between[i]]])==contrast.levels[[j]][2]

      if (length(other)>0){
        for (k in other){
          positive <- positive & as.character(analysis.data[[k]])==as.character(result[[k]][j])
          negative <- negative & as.character(analysis.data[[k]])==as.character(result[[k]][j])
        }
      }

      paste(sum(positive), 'vs', sum(negative))
    }, character(1))

    result.comparison <- if (length(other)==0){
      rep(paste('within', between[i], '(unpaired)'), nrow(result))
    } else {
      apply(result[, other, drop=FALSE], 1, function(x){
        paste('within', paste(other, x, collapse=' '), '(unpaired)')
      })
    }

    results <- rbind(results, result[, c('contrast', 'estimate', 'SE', 'df',
      't.ratio', 'p.value', 'lower.CL', 'upper.CL')])
    comparison <- c(comparison, result.comparison)
    n <- c(n, result.n)
  }

  adjusted <- p.adjust.ff(results$p.value, method=p.adjust.method, alpha=alpha)

  model.name <- paste0(
    'lm: ',
    paste(deparse(formula(fit)), collapse=' '),
    '; residual degrees-of-freedom inference'
    )

  CI.probability <- c((1-CI.level)/2, 1-(1-CI.level)/2)
  CI.names <- paste0('CI (', formatC(CI.probability, format='f', digits=3), ')')

  lm_summary <- data.frame(
    model=model.name,
    data=data.name,
    parameter=rep(response, nrow(results)),
    comparison=comparison,
    contrast=gsub(' - ', ' vs ', results$contrast, fixed=TRUE),
    n=n,
    test='Linear model with residual degrees-of-freedom inference',
    alternative='two.sided',
    'test stat'='t',
    df=results$df,
    stat=results$t.ratio,
    estimate=results$estimate,
    se=results$SE,
    lower=results$lower.CL,
    upper=results$upper.CL,
    'p value'=results$p.value,
    'p adjusted'=adjusted[,'p.adj'],
    check.names=FALSE
    )

  names(lm_summary)[names(lm_summary)=='lower'] <- CI.names[1]
  names(lm_summary)[names(lm_summary)=='upper'] <- CI.names[2]

  rownames(lm_summary) <- NULL

  lm_summary
}


LMEsummary <- function(fit, between, within, subject, data.name,
  df=c('kenward-roger', 'satterthwaite'), p.adjust.method='Sidak-Holm',
  alpha=0.05, CI=95){

  df <- match.arg(df)
  CI.level <- CI/100

  df.label <- if (df=='kenward-roger'){
    'Kenward-Roger'
  } else {
    'Satterthwaite'
  }

  analysis.data <- model.frame(fit)
  response <- all.vars(formula(fit))[1]

  if (is.null(between)){

    emm.paired <- emmeans::emmeans(fit,
      specs=as.formula(paste('~', within)), lmer.df=df)

    comparisons.paired <- emmeans::contrast(emm.paired,
      method='revpairwise', adjust='none')

    results.paired <- as.data.frame(summary(comparisons.paired,
      infer=c(TRUE,TRUE), by=NULL, adjust='none', level=CI.level))

    adjusted.paired <- p.adjust.ff(results.paired$p.value,
      method=p.adjust.method, alpha=alpha)

    subject.values <- as.character(analysis.data[[subject]])
    within.values <- as.character(analysis.data[[within]])
    paired.levels <- strsplit(as.character(results.paired$contrast),
      ' - ', fixed=TRUE)

    n.paired <- vapply(seq_len(nrow(results.paired)), function(i){

      positive.subjects <- unique(subject.values[
        within.values==paired.levels[[i]][1]
        ])

      negative.subjects <- unique(subject.values[
        within.values==paired.levels[[i]][2]
        ])

      as.character(length(intersect(positive.subjects, negative.subjects)))
    }, character(1))

    model.name <- paste0(
      'lmer (REML): ',
      paste(deparse(formula(fit)), collapse=' '),
      '; ',
      df.label,
      ' inference'
      )

    CI.probability <- c((1-CI.level)/2, 1-(1-CI.level)/2)
    CI.names <- paste0('CI (',
      formatC(CI.probability, format='f', digits=3), ')')

    lme_summary <- data.frame(
      model=model.name,
      data=data.name,
      parameter=rep(response, nrow(results.paired)),
      comparison=rep(paste('within', within, '(paired)'), nrow(results.paired)),
      contrast=gsub(' - ', ' vs ', results.paired$contrast, fixed=TRUE),
      n=n.paired,
      test=paste('Linear mixed-effects model with', df.label, 'inference'),
      alternative='two.sided',
      'test stat'='t',
      df=results.paired$df,
      stat=results.paired$t.ratio,
      estimate=results.paired$estimate,
      se=results.paired$SE,
      lower=results.paired$lower.CL,
      upper=results.paired$upper.CL,
      'p value'=results.paired$p.value,
      'p adjusted'=adjusted.paired[,'p.adj'],
      check.names=FALSE
      )

    names(lme_summary)[names(lme_summary)=='lower'] <- CI.names[1]
    names(lme_summary)[names(lme_summary)=='upper'] <- CI.names[2]

    rownames(lme_summary) <- NULL

    return(lme_summary)
  }

  emm.unpaired <- emmeans::emmeans(fit,
    specs=as.formula(paste('~', between, '|', within)), lmer.df=df)

  comparisons.unpaired <- emmeans::contrast(emm.unpaired,
    method='revpairwise', adjust='none')

  results.unpaired <- as.data.frame(summary(comparisons.unpaired,
    infer=c(TRUE,TRUE), by=NULL, adjust='none', level=CI.level))

  adjusted.unpaired <- p.adjust.ff(results.unpaired$p.value,
    method=p.adjust.method, alpha=alpha)

  results.unpaired$'p adjusted' <- adjusted.unpaired[,'p.adj']

  emm.paired <- emmeans::emmeans(fit,
    specs=as.formula(paste('~', within, '|', between)), lmer.df=df)

  comparisons.paired <- emmeans::contrast(emm.paired,
    method='revpairwise', adjust='none')

  results.paired <- as.data.frame(summary(comparisons.paired,
    infer=c(TRUE,TRUE), by=NULL, adjust='none', level=CI.level))

  adjusted.paired <- p.adjust.ff(results.paired$p.value,
    method=p.adjust.method, alpha=alpha)

  results.paired$'p adjusted' <- adjusted.paired[,'p.adj']

  subject.values <- as.character(analysis.data[[subject]])
  between.values <- as.character(analysis.data[[between]])
  within.values <- as.character(analysis.data[[within]])

  unpaired.levels <- strsplit(as.character(results.unpaired$contrast),
    ' - ', fixed=TRUE)

  n.unpaired <- vapply(seq_len(nrow(results.unpaired)), function(i){

    within.level <- as.character(results.unpaired[[within]][i])

    positive.subjects <- unique(subject.values[
      within.values==within.level &
        between.values==unpaired.levels[[i]][1]
      ])

    negative.subjects <- unique(subject.values[
      within.values==within.level &
        between.values==unpaired.levels[[i]][2]
      ])

    paste(length(positive.subjects), 'vs', length(negative.subjects))
  }, character(1))

  paired.levels <- strsplit(as.character(results.paired$contrast),
    ' - ', fixed=TRUE)

  n.paired <- vapply(seq_len(nrow(results.paired)), function(i){

    between.level <- as.character(results.paired[[between]][i])

    positive.subjects <- unique(subject.values[
      between.values==between.level &
        within.values==paired.levels[[i]][1]
      ])

    negative.subjects <- unique(subject.values[
      between.values==between.level &
        within.values==paired.levels[[i]][2]
      ])

    as.character(length(intersect(positive.subjects, negative.subjects)))
  }, character(1))

  model.name <- paste0(
    'lmer (REML): ',
    paste(deparse(formula(fit)), collapse=' '),
    '; ',
    df.label,
    ' inference'
    )

  CI.probability <- c((1-CI.level)/2, 1-(1-CI.level)/2)
  CI.names <- paste0('CI (',
    formatC(CI.probability, format='f', digits=3), ')')

  lme_summary <- data.frame(
    model=model.name,
    data=data.name,
    parameter=rep(response, nrow(results.unpaired)+nrow(results.paired)),
    comparison=c(
      paste('within', within, results.unpaired[[within]], '(unpaired)'),
      paste('within', between, results.paired[[between]], '(paired)')
      ),
    contrast=gsub(' - ', ' vs ',
      c(results.unpaired$contrast, results.paired$contrast), fixed=TRUE),
    n=c(n.unpaired, n.paired),
    test=paste('Linear mixed-effects model with', df.label, 'inference'),
    alternative='two.sided',
    'test stat'='t',
    df=c(results.unpaired$df, results.paired$df),
    stat=c(results.unpaired$t.ratio, results.paired$t.ratio),
    estimate=c(results.unpaired$estimate, results.paired$estimate),
    se=c(results.unpaired$SE, results.paired$SE),
    lower=c(results.unpaired$lower.CL, results.paired$lower.CL),
    upper=c(results.unpaired$upper.CL, results.paired$upper.CL),
    'p value'=c(results.unpaired$p.value, results.paired$p.value),
    'p adjusted'=c(
      results.unpaired$'p adjusted',
      results.paired$'p adjusted'
      ),
    check.names=FALSE
    )

  names(lme_summary)[names(lme_summary)=='lower'] <- CI.names[1]
  names(lme_summary)[names(lme_summary)=='upper'] <- CI.names[2]

  rownames(lme_summary) <- NULL

  lme_summary
}

BayesMCtest <- function(formula, df, na_rm_subjects=TRUE, contr=NULL, mu=0,
  rscale='medium'){

  formula.string <- paste(deparse(formula), collapse='')
  has.error <- grepl('Error', formula.string)

  if (has.error){
    error.part <- sub('.*Error\\((.*)\\).*', '\\1', formula.string)
    subject.name <- strsplit(error.part, '/')[[1]][1]
    subject.name <- gsub('[[:space:]]', '', subject.name)
    main.formula.string <- sub('\\+\\s*Error\\(.*\\)', '', formula.string)
    main.formula <- as.formula(main.formula.string)
  } else {
    subject.name <- NULL
    main.formula <- formula
  }

  formula.variables <- all.vars(main.formula)
  response.name <- formula.variables[1]
  predictor.names <- formula.variables[-1]

  if (!(length(predictor.names) %in% c(1,2))){
    stop('formula must contain one or two predictors')
  }

  required.names <- c(response.name, predictor.names, subject.name)

  if (!all(required.names %in% names(df))){
    stop('formula or subject variables are missing from df')
  }

  if (na_rm_subjects && !is.null(subject.name)){
    remove.subject <- ave(is.na(df[[response.name]]), df[[subject.name]], FUN=any)
    df <- df[!remove.subject,,drop=FALSE]
  }

  analysis.data <- df[,required.names,drop=FALSE]
  analysis.data <- analysis.data[complete.cases(analysis.data),,drop=FALSE]

  if (!nrow(analysis.data)){
    stop('no complete observations are available')
  }

  for (predictor in predictor.names){
    analysis.data[[predictor]] <- droplevels(factor(analysis.data[[predictor]]))

    if (nlevels(analysis.data[[predictor]])<2){
      stop('each predictor must contain at least two levels')
    }
  }

  if (!is.null(subject.name)){
    analysis.data[[subject.name]] <- droplevels(factor(analysis.data[[subject.name]]))
  }

  if (length(predictor.names)==1){
    cell <- analysis.data[[predictor.names]]
  } else {
    cell <- do.call(interaction, c(analysis.data[predictor.names],
      list(sep=':', drop=TRUE, lex.order=FALSE)))
  }

  cell <- droplevels(factor(cell))
  cell.levels <- levels(cell)

  cell.grid <- analysis.data[
    match(cell.levels, as.character(cell)),
    predictor.names,
    drop=FALSE
    ]

  rownames(cell.grid) <- cell.levels


  # Create default adjacent-level contrasts

  if (is.null(contr)){

    result.contrasts <- list()

    if (length(predictor.names)==1){

      predictor.levels <- levels(analysis.data[[predictor.names]])

      for (i in seq_len(length(predictor.levels)-1)){
        contrast.row <- rep(0, length(cell.levels))
        contrast.row[match(predictor.levels[i], cell.levels)] <- -1
        contrast.row[match(predictor.levels[i+1], cell.levels)] <- 1

        result.contrasts[[length(result.contrasts)+1]] <- contrast.row
      }

    } else {

      predictor.order <- predictor.names

      if (!is.null(subject.name)){
        subject.rows <- split(seq_len(nrow(analysis.data)),
          analysis.data[[subject.name]])

        varies.within.subject <- vapply(predictor.names, function(predictor){
          any(vapply(subject.rows, function(i){
            length(unique(analysis.data[[predictor]][i]))>1
            }, logical(1)))
          }, logical(1))

        if (sum(varies.within.subject)==1){
          predictor.order <- c(
            predictor.names[!varies.within.subject],
            predictor.names[varies.within.subject]
            )
        }
      }

      for (target.name in predictor.order){

        conditioning.name <- setdiff(predictor.names, target.name)
        conditioning.levels <- levels(analysis.data[[conditioning.name]])
        target.levels <- levels(analysis.data[[target.name]])

        for (conditioning.level in conditioning.levels){

          for (i in seq_len(length(target.levels)-1)){

            negative.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i]
              ]

            positive.cell <- rownames(cell.grid)[
              cell.grid[[conditioning.name]]==conditioning.level &
                cell.grid[[target.name]]==target.levels[i+1]
              ]

            if (length(negative.cell)==1 && length(positive.cell)==1){
              contrast.row <- rep(0, length(cell.levels))
              contrast.row[match(negative.cell, cell.levels)] <- -1
              contrast.row[match(positive.cell, cell.levels)] <- 1

              result.contrasts[[length(result.contrasts)+1]] <- contrast.row
            }
          }
        }
      }
    }

    contr <- do.call(rbind, result.contrasts)
    colnames(contr) <- cell.levels

  } else {

    if (is.numeric(contr) && is.null(dim(contr))){
      contr <- matrix(contr, nrow=1)
    }

    if (!is.matrix(contr) || !is.numeric(contr) || !nrow(contr) ||
      anyNA(contr) || any(!is.finite(contr)) ||
      ncol(contr)!=length(cell.levels)){
      stop('contr must be a finite numeric matrix with one column per factor-level combination')
    }

    if (!is.null(colnames(contr))){

      if (!setequal(colnames(contr), cell.levels)){
        stop(paste('contr column names must match:',
          paste(cell.levels, collapse=', ')))
      }

      contr <- contr[,cell.levels,drop=FALSE]

    } else {
      colnames(contr) <- cell.levels
    }
  }

  pairwise.contrast <- apply(contr, 1, function(x){
    nonzero <- sort(as.numeric(x[x!=0]))
    length(nonzero)==2 && all(nonzero==c(-1,1))
    })

  if (any(!pairwise.contrast)){
    stop('each row of contr must contain one -1, one 1 and zeros elsewhere')
  }


  # Perform Bayesian paired or independent-samples t-tests

  results <- vector('list', nrow(contr))
  BayesFactor.results <- vector('list', nrow(contr))

  for (i in seq_len(nrow(contr))){

    negative.cell <- cell.levels[contr[i,]==-1]
    positive.cell <- cell.levels[contr[i,]==1]
    keep <- cell %in% c(negative.cell, positive.cell)

    test.data <- analysis.data[keep,,drop=FALSE]

    if (!is.null(subject.name)){
      test.data[[subject.name]] <- droplevels(test.data[[subject.name]])
    }

    test.data$.cell <- factor(as.character(cell[keep]),
      levels=c(negative.cell, positive.cell))

    negative.data <- test.data[
      test.data$.cell==negative.cell,
      ,
      drop=FALSE
      ]

    positive.data <- test.data[
      test.data$.cell==positive.cell,
      ,
      drop=FALSE
      ]

    paired <- FALSE

    if (!is.null(subject.name)){
      negative.subjects <- unique(as.character(
        negative.data[[subject.name]]
        ))

      positive.subjects <- unique(as.character(
        positive.data[[subject.name]]
        ))

      common.subjects <- intersect(negative.subjects, positive.subjects)

      if (setequal(negative.subjects, positive.subjects)){
        paired <- TRUE
      } else if (length(common.subjects)){
        stop('a contrast cannot contain a mixture of paired and unpaired subjects')
      }
    }

    if (paired){

      pair.table <- table(test.data[[subject.name]], test.data$.cell)

      if (any(pair.table!=1)){
        stop('paired comparisons require one complete observation per subject and contrast level')
      }

      negative.data <- negative.data[
        order(negative.data[[subject.name]]),
        ,
        drop=FALSE
        ]

      positive.data <- positive.data[
        order(positive.data[[subject.name]]),
        ,
        drop=FALSE
        ]

      n.output <- as.character(nrow(positive.data))

    } else {

      if (!is.null(subject.name) &&
        (anyDuplicated(negative.data[[subject.name]]) ||
        anyDuplicated(positive.data[[subject.name]]))){
        stop('unpaired comparisons require one observation per subject and contrast level')
      }

      n.output <- paste(nrow(positive.data), 'vs', nrow(negative.data))
    }

    positive.values <- positive.data[[response.name]]
    negative.values <- negative.data[[response.name]]

    test <- BayesFactor::ttestBF(
      x=positive.values,
      y=negative.values,
      paired=paired,
      mu=mu,
      rscale=rscale
      )

    BF10 <- unname(BayesFactor::extractBF(test, onlybf=TRUE)[1])

    differing <- predictor.names[
      vapply(predictor.names, function(x){
        as.character(cell.grid[negative.cell,x])!=
          as.character(cell.grid[positive.cell,x])
        }, logical(1))
      ]

    same <- setdiff(predictor.names, differing)

    if (length(predictor.names)==1){

      parameter <- response.name

      comparison <- paste(
        'within',
        predictor.names,
        if (paired) '(paired)' else '(unpaired)'
        )

      contrast.name <- paste(
        as.character(cell.grid[positive.cell,differing]),
        'vs',
        as.character(cell.grid[negative.cell,differing])
        )

    } else if (length(differing)==1){

      parameter <- response.name

      comparison <- paste(
        'within',
        same,
        as.character(cell.grid[positive.cell,same]),
        if (paired) '(paired)' else '(unpaired)'
        )

      contrast.name <- paste(
        as.character(cell.grid[positive.cell,differing]),
        'vs',
        as.character(cell.grid[negative.cell,differing])
        )

    } else {

      parameter <- response.name

      comparison <- paste(
        'between cells',
        if (paired) '(paired)' else '(unpaired)'
        )

      contrast.name <- paste(positive.cell, 'vs', negative.cell)
    }

    prior.name <- if (is.character(rscale)){
      paste('Cauchy,', rscale, 'scale')
    } else {
      paste('Cauchy, rscale =', rscale)
    }

    results[[i]] <- data.frame(
      parameter=parameter,
      comparison=comparison,
      contrast=contrast.name,
      n=as.character(n.output),
      test=if (paired){
        'Bayesian paired t-test'
      } else {
        'Bayesian independent-samples t-test'
      },
      alternative='two.sided',
      prior=prior.name,
      BF10=BF10,
      family=if (paired) 'paired' else 'unpaired',
      check.names=FALSE
      )

    BayesFactor.results[[i]] <- test
  }

  bayes_summary <- do.call(rbind, results)
  rownames(bayes_summary) <- NULL

  family.order <- match(bayes_summary$family, c('unpaired','paired'))
  result.order <- order(family.order)

  bayes_summary <- bayes_summary[result.order,,drop=FALSE]
  BayesFactor.results <- BayesFactor.results[result.order]

  bayes_summary$family <- NULL
  rownames(bayes_summary) <- NULL

  attr(bayes_summary, 'BayesFactor') <- BayesFactor.results

  bayes_summary
}
