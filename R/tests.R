# Martingale Tests and Monte Carlo convergence --------------------------------------------------

#' Convergence of Monte Carlo prices
#' 
#' @description
#' Analyzes the convergence of Monte Carlo prices by computing average prices and confidence intervals
#' as the number of simulations increases.
#' 
#' @param r Interest rate time series or constant
#' @param X Price time series
#' @param maturity Maturity time point
#' @param plot Logical indicating whether to plot the convergence (default: TRUE)
#' @param ... Additional arguments passed to plotting functions
#' 
#' @details
#' The function computes the average price and confidence intervals for increasing numbers of simulations
#' to analyze the convergence of Monte Carlo estimates. It discounts the price series using the provided
#' interest rate and evaluates at the specified maturity.
#' 
#' If plot=TRUE, it creates a plot showing:
#' \itemize{
#'   \item Confidence intervals as shaded area
#'   \item Average price as a blue line
#'   \item X-axis showing number of simulations
#'   \item Y-axis showing Monte Carlo estimated prices
#' }
#' 
#' @return A list containing:
#' \itemize{
#'   \item avg.price: Vector of average prices for different numbers of simulations
#'   \item conf.int: Matrix of confidence intervals (lower and upper bounds)
#' }
#' 
#' @examples
#' # Generate sample data
#' r <- 0.05
#' X <- ts(matrix(rnorm(1000), 100, 10), start = 0, deltat = 0.1)
#' 
#' # Analyze convergence
#' convergence <- esgmccv(r, X, maturity = 1)
#' 
#' @export

esgmccv <- function(r, X, maturity, plot = TRUE, ...)
{
  if(missing(r) || missing(X) || missing(maturity))
    stop("'r', 'X', and 'maturity' must be provided")
  
  bool.X <- (is.ts(X) && tsp(X)[1] > 0)
  
  if(bool.X)
    maturity <- maturity - deltat(X)
  
  Y <- esgdiscountfactor(r, X)
  Z <- window(Y, start = maturity, end = maturity)
  N <- length(Z)  
  
  avg.price <- sapply(2:N, function(x) mean(Z[1:x]))
  conf.int <- t(sapply(2:N, function(x) t.test(Z[1:x])$conf.int[1:2]))
  colnames(conf.int) <- c("lower bound", "upper bound")
  
  if (plot == TRUE)
  {
    x <- 2:N
    matplot(x, conf.int, type = 'l', xlab = "number of simulations", 
            ylab = "monte carlo estim. price", lty = c(1, 1), lwd = 2, ...)
    polygon(c(x, rev(x)), 
            c(as.vector(conf.int[, 2]), rev(as.vector(conf.int[, 1]))), 
            col = "lightyellow", border = FALSE)
    lines(x, avg.price, col = "blue")    
  }  
  invisible(list(avg.price = avg.price,
                 conf.int = conf.int))                  
}


#' Test for martingale property and market consistency
#' 
#' @description
#' Performs statistical tests to verify if a time series follows a martingale property
#' and maintains market consistency. Three different testing methods are available.
#' 
#' @param r Risk-free rate (scalar or time series)
#' @param X Time series of asset prices or values
#' @param p0 Initial price(s) for comparison
#' @param alpha Significance level for hypothesis tests (default: 0.05)
#' @param method Testing method:
#'   \itemize{
#'     \item "onevsone": One-to-one comparison test
#'     \item "trend": Trend analysis test
#'     \item "ratio": Ratio-based test
#'   }
#' 
#' @details
#' The function implements three different approaches to test for martingale properties:
#' \itemize{
#'   \item "onevsone": Compares each observation against its expected value
#'   \item "trend": Analyzes trends in the martingale differences
#'   \item "ratio": Tests the ratio between discounted values and initial prices
#' }
#' 
#' The test results help determine if the time series maintains the martingale property
#' and market consistency under the risk-neutral measure.
#' 
#' @return A list containing test results, which may include:
#' \itemize{
#'   \item t-statistic and p-value for ratio test
#'   \item trend analysis results
#'   \item mean and standard deviation of martingale differences
#' }
#' 
#' @examples
#' # Generate sample data
#' r <- 0.05
#' X <- ts(matrix(rnorm(1000), 100, 10), start = 0, deltat = 0.1)
#' p0 <- 100
#' 
#' # Perform martingale test using ratio method
#' result <- esgmartingaletest(r, X, p0, method = "ratio")
#' 
#' # Perform trend analysis
#' trend_result <- esgmartingaletest(r, X, p0, method = "trend")
#' 
#' @export
esgmartingaletest <- function(r, X, p0, alpha = 0.05, 
                              method=c("onevsone", 
                                       "trend",
                                       "ratio"))
{   
  method <- match.arg(method)
  delta_X <- deltat(X)
  
  if (length(r) == 1) 
  {
    r <- ts(data = matrix(data = r, nrow = nrow(X), ncol = ncol(X)), 
            start = start(X), deltat = delta_X)
  }  
  nrow.r <- nrow(r)
  ncol.r <- ncol(r)  
  delta_r <- deltat(r)  
  nb.p0 <- length(p0)
  
  if (nb.p0 == 1)
  {
    Y <- ts(matrix(data = rep(p0, nrow.r*ncol.r), 
                   nrow = nrow.r, ncol = ncol.r), 
            start = start(r), 
            deltat = delta_X)
  }
  else
  {
    Y <- ts(data = replicate(ncol.r, p0), 
            start = start(r), 
            deltat = delta_X)    
  }
  
  delta_Y <- delta_X  
  Dt <- esgdiscountfactor(r, X)
  if (method == "ratio")
  {
    return(t.test(x=-log(Dt/Y), conf.level = (1 - alpha)))
  }
  MartingaleDiff <- Dt - Y
  if (method == "trend")
  {
    return(suppressWarnings(martingale_test(MartingaleDiff, level = 100*(1 - alpha))))
  }
  
  n <- ncol(MartingaleDiff) 
  meanMartingaleDiff <- rowMeans(MartingaleDiff[-1, ])  
  sdMartingaleDiff <- apply(MartingaleDiff[-1, ], 1, sd)
  if (is.ts(MartingaleDiff))
  {
    meanMartingaleDiff <- ts(meanMartingaleDiff, 
                             deltat = deltat(MartingaleDiff),
                             end = end(MartingaleDiff))  
    sdMartingaleDiff <- ts(sdMartingaleDiff, 
                           deltat = deltat(MartingaleDiff),
                           end = end(MartingaleDiff))
  }
  qtStudent <- qt(p = 1 - alpha/2, df = (n - 1))  
  stat_t <- meanMartingaleDiff/(sdMartingaleDiff/sqrt(n))  
  p_value <- pt(q = abs(stat_t), df = (n - 1), lower.tail = F) + 
    pt(q = -abs(stat_t), df = (n - 1))  
  horizon <- end(r)[1]
  mc.ci <- list(t = stat_t, 
                p.value = p_value, 
                samplemean = meanMartingaleDiff, 
                conf.int = ts(cbind(c(0, meanMartingaleDiff - qtStudent * sdMartingaleDiff/sqrt(n)), 
                                    c(0, meanMartingaleDiff + qtStudent * sdMartingaleDiff/sqrt(n))),
                              deltat = deltat(MartingaleDiff),
                              end = end(MartingaleDiff)),
                truemean = rep.int(0, dim(MartingaleDiff)[1]), 
                true_prices = Y[1:nrow(MartingaleDiff), 1], mc.prices = rowMeans(Dt)) 
  
  start_Y <- start(Y)  
  
  mat.ci <- ts(mc.ci$conf.int, 
               deltat = deltat(MartingaleDiff),
               end = end(MartingaleDiff))  
  colnames(mat.ci) <- c("c.i lower bound", "c.i upper bound")  
  t_val_p_val <- ts(cbind(mc.ci$t, mc.ci$p.value), 
                    deltat = deltat(MartingaleDiff),
                    end = end(MartingaleDiff))  
  colnames(t_val_p_val) <- c("t", "p-value")
  
  cat("\n")
  cat(" martingale '1=1' one Sample t-test", "\n")
  cat("\n")
  cat(" alternative hypothesis: true mean of the martingale difference is not equal to 0", 
      "\n")
  cat("\n")
  cat("df = ", n - 1)
  cat("\n")
  print(t_val_p_val)
  cat("\n")
  cat((1 - alpha) * 100, "percent confidence intervals for the mean :", "\n")
  print(mat.ci)
  invisible(mc.ci)
}


# Correlation test for shocks --------------------------------------------------------


#' Correlation test for shocks
#' 
#' @description
#' Performs correlation tests between two sets of shocks at each time point.
#' 
#' @param x A list containing two time series objects to test for correlation
#' @param alternative A character string specifying the alternative hypothesis:
#'   "two.sided" (default), "less", or "greater"
#' @param method A character string indicating which correlation coefficient is to be used:
#'   "pearson" (default), "kendall", or "spearman"
#' @param conf.level Confidence level for the returned confidence interval
#' 
#' @details
#' The function performs correlation tests between two sets of shocks at each time point.
#' It returns both the correlation estimates and confidence intervals for the correlation
#' coefficient at each time point.
#' 
#' @return A list containing:
#' \itemize{
#'   \item cor.estimate: A time series of correlation estimates
#'   \item conf.int: A time series matrix containing the lower and upper bounds of the confidence intervals
#' }
#' 
#' @examples
#' # Generate sample shocks
#' x <- ts(matrix(rnorm(1000), 100, 10))
#' y <- ts(matrix(rnorm(1000), 100, 10))
#' 
#' # Test correlation between shocks
#' result <- esgcortest(list(x, y))
#' 
#' # Plot correlation estimates with confidence intervals
#' plot(result$cor.estimate)
#' lines(result$conf.int[,1], col="red")
#' lines(result$conf.int[,2], col="red")
#' 
#' @export
esgcortest <- function(x, alternative = c("two.sided", "less", "greater"),
                       method = c("pearson", "kendall", "spearman"),
                       conf.level = 0.95)
{
  y <- x[[2]]
  x <- x[[1]]
  delta_x <- deltat(x)
  delta_y <- deltat(y)
  
  if (prod(dim(x) == dim(y)) == 0) stop("We must have dim(x) == dim(y)")  
  if (delta_x != delta_y) stop("We must have deltat(x) == deltat(y)")  
  alternative <- match.arg(alternative)  
  method <- match.arg(method)  
  nbdates <- dim(x)[1]  
  return(list(cor.estimate = ts(sapply(1:nbdates, function(i) cor(x[i,], y[i,], 
                                                                  method = method)), 
                                start = delta_x, deltat = delta_x), 
              conf.int = ts(t(sapply(1:nbdates, function(i) cor.test(x[i,], y[i,],
                                                                     alternative = alternative, method = method, 
                                                                     conf.level = conf.level)$conf.int)), start = delta_x, 
                            deltat = delta_x)))
}


martingale_test <- function(X, level=95) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Compute Y as the difference between the last two rows
  Y <- X[n, ] - X[n-1, ]
  
  # Create predictors (past values of X up to t_n-1)
  X_past <- t(X[seq_len(n-1), ])
  
  # Fit the multiple regression model
  model <- lm(Y ~ X_past - 1)
  
  # Compute F-statistic
  regression_summary <- summary(model)
  f_stat <- regression_summary$fstatistic
  F_obs <- f_stat[1]
  df1 <- f_stat[2]
  df2 <- f_stat[3]
  
  # Compute p-value
  p_value <- pf(F_obs, df1, df2, lower.tail = FALSE)
  
  # Compute critical value at 5% significance level
  F_critical <- qf(level/100, df1, df2)
  
  # Extract residuals
  residuals <- model$residuals
  
  # Perform ADF test for stationarity on the residuals
  adf_result <- tseries::adf.test(residuals)
  
  # Perform Ljung-Box test for autocorrelation on the residuals
  lb_test <- stats::Box.test(residuals, lag = 1L, 
                             type = "Ljung-Box")
  lb_p_value <- lb_test$p.value
  
  return(list(
    model = model,
    confint = confint(model),
    regression_summary = regression_summary,
    F_statistic = F_obs,
    F_critical_value = F_critical,
    F_p_value = p_value,
    ADF_p_value = 1 - adf_result$p.value,
    Ljung_Box_p_value = lb_p_value
  ))
}

