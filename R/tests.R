# Martingale Tests and Monte Carlo convergence --------------------------------------------------

#'@title Stochastic discount factors or discounted values
#'
#'@description 
#'
#'This function provides calculation of stochastic discount factors 
#'or discounted values
#'
#'@details
#'
#'The function result is : 
#'
#'\deqn{X_t exp(-\int_0^t r_s ds)}
#'
#'where \eqn{X_t} is an asset value at a given maturity \eqn{t}, and 
#'\eqn{(r_s)_s} is the risk-free rate.
#'
#'@param r the short rate, a \code{numeric} (constant rate) or a time series object
#'
#'@param X the asset's price, a \code{numeric} (constant payoff or asset price) or a time series 
#'object
#'
#'@seealso \code{\link{esgmcprices}}, \code{\link{esgmccv}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@examples 
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma_v <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma_v
#'
#'# OU
#'r <- simdiff(n = 10, horizon = 5, 
#'                frequency = "quart",  
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'# Stochastic discount factors
#' esgdiscountfactor(r, 1)
esgdiscountfactor <- function(r, X)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  length.r <- length(r)
  length.X <- length(X)
  start.r <- start(r)
  deltat.r <- deltat(r)
  start.X <- start(X)
  deltat.X <- deltat(X)
  
  if(length.r == 1 && length.X == 1) 
  {
    return(X*exp(-r)) 
  }
  
  if(length.r == 1 && length.X != 1) 
  {
    r <- ts(matrix(r, nrow(X), ncol(X)), 
            start = start.X, deltat = deltat.X)
    
    if(tsp(X)[1] > 0)
    {
      return(ts(X*exp(-apply(r, 2, cumsum)*deltat.X), 
                start = 0, 
                deltat = deltat.X))      
    }
    else
    {
      Int_r <- exp(-apply(r, 2, cumsum)*deltat.X)
      return(ts(X*rbind(rep(1, ncol(X)), 
                        Int_r[1:(nrow(X)-1), ]),                                           , 
                start = 0, 
                deltat = deltat.X))
    }
  } 
  
  if(length.r != 1 && length.X == 1)
  {
    X <- ts(matrix(X, nrow(r), ncol(r)), 
            start = start.r, deltat = deltat.r)
    
    return(ts(X*exp(-apply(r, 2, cumsum)*deltat.r),                                          
              start = 0, 
              deltat = deltat.r))
  }
  
  if(length.r != 1 && length.X != 1)
  {
    if(tsp(X)[1] > 0)
    {
      return(suppressWarnings(ts(X*window(exp(-apply(r, 2, cumsum)*deltat.r), 
                                          start = start.X, 
                                          deltat = deltat.X), 
                                 start = 0, 
                                 deltat = deltat.X)))
    }
    else
    {
      Int_r <- ts(exp(-apply(r, 2, cumsum)*deltat.r), deltat = deltat.r)
      return(suppressWarnings(ts(X*rbind(rep(1, ncol(X)), 
                                         Int_r[1:(nrow(X)-1), ]),                                           , 
                                 start = 0, 
                                 deltat = deltat.X)))
    }
  }
}
esgdiscountfactor <- cmpfun(esgdiscountfactor)

#'@title
#'
#'Estimation of discounted asset prices
#'
#'@description
#'
#'This function computes estimators (sample mean) of 
#'
#'\deqn{E[X_T exp(-\int_0^T r_s ds)]}
#'
#'where \eqn{X_T} is an asset value at given maturities \eqn{T}, and 
#'\eqn{(r_s)_s} is the risk-free rate.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X asset prices obtained with \code{\link{simdiff}}
#'
#'@param maturity the corresponding maturity (optional). If missing, all the maturities 
#'available in \code{X} are used.
#'
#'@seealso \code{\link{esgdiscountfactor}}, \code{\link{esgmccv}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@examples
#'
#'# GBM
#'
#'r <- 0.03
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'
#'# monte carlo prices
#'esgmcprices(r, sim.GBM)
#'
#'# monte carlo price for a given maturity
#'esgmcprices(r, sim.GBM, 2)
#'
esgmcprices <- function(r, X, maturity = NULL)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  maturity.out <- maturity
  
  if(is.ts(X) && tsp(X)[1] > 0 && !is.null(maturity))
  {
    maturity.out <- maturity - deltat(X)
  }
  
  Y <- esgdiscountfactor(r, X)
  
  if(length(r) == 1 && length(X) == 1) 
  {
    return(Y) 
  }
  
  Z <- ts(rowMeans(Y), start = start(Y), deltat = deltat(Y))
  
  if(!is.null(maturity))
  {
    return(window(Z, start = maturity.out, end = maturity.out))
  }
  else  
  {
    return(Z)
  }
  
}
esgmcprices <- cmpfun(esgmcprices)

#'@title
#'
#'Convergence of Monte Carlo prices
#'
#'@description
#'
#'This function computes and plots confidence intervals around the estimated 
#'average price, as functions of the number of simulations. 
#'
#'@details
#'
#'Studying the convergence of the sample mean of :
#'
#'\deqn{E[X_T exp(-\int_0^T r_s ds)]}
#'
#'towards its true value.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X asset prices obtained with \code{\link{simdiff}}
#'
#'@param maturity the corresponding maturity (optional). If missing, all the maturities 
#'available in \code{X} are used.
#'
#'@param plot if \code{TRUE} (default), a plot of the convergence is displayed.
#'
#'@param ... additional parameters provided to \code{\link{matplot}}
#'
#'@return a list with estimated average prices and the confidence intervals around them.
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@examples
#'
#'r <- 0.03
#'
#'set.seed(1)
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'
#'# monte carlo prices
#'esgmcprices(r, sim.GBM)
#'
#'# convergence to a specific price
#'(esgmccv(r, sim.GBM, 2))
#'
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
esgmccv <- cmpfun(esgmccv)


#'@title Martingale and market consistency tests
#'
#'@description
#'
#'This function performs martingale and market consistency (t-)tests.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X a time series object, containing payoffs or projected asset values.
#' 
#'@param p0 a \code{numeric} or a vector or a univariate time series containing  
#'initial price(s) of an asset. 
#'
#'@param alpha 1 - confidence level for the test. Default value is 0.05.
#'
#'@return The function result can be just displayed. Otherwise, you can get a list 
#'by an assignation, containing (for each maturity) : 
#'\itemize{
#'\item the Student t values 
#'\item the p-values 
#'\item the estimated mean
#'of the martingale difference
#'\item  Monte Carlo prices
#'}
#'
#'@export
#'
#'@seealso \code{\link{esgplotbands}}
#'
#'@author Thierry Moudiki
#'
#'@examples
#'
#'r0 <- 0.03
#'S0 <- 100
#'
#'set.seed(10)
#'eps0 <- simshocks(n = 100, horizon = 3, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 3, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = S0, theta1 = r0, theta2 = 0.1, 
#'                eps = eps0)
#' 
#'mc.test <- esgmartingaletest(r = r0, X = sim.GBM, p0 = S0, 
#'alpha = 0.05)                               
#'esgplotbands(mc.test)                
#'
esgmartingaletest <- function(r, X, p0, alpha = 0.05)
{   
  delta_X <- deltat(X)
  
  if (length(r) == 1) 
  {
    r <- ts(data = matrix(data = r, nrow = nrow(X), ncol = ncol(X)), 
            start = 0, deltat = delta_X)
  }  
  nrow.r <- nrow(r)
  ncol.r <- ncol(r)  
  delta_r <- deltat(r)  
  nb.p0 <- length(p0)
  
  if (nb.p0 == 1)
  {
    Y <- ts(matrix(data = rep(p0, nrow.r*ncol.r), 
                   nrow = nrow.r, ncol = ncol.r), 
            start = 0, deltat = delta_X)
  }
  else
  {
    Y <- ts(data = replicate(ncol.r, p0), 
            start = 0, deltat = delta_X)    
  }
  
  delta_Y <- delta_X  
  Dt <- esgdiscountfactor(r, X)
  MartingaleDiff <- Dt - Y
  
  n <- ncol(MartingaleDiff) 
  meanMartingaleDiff <- rowMeans(MartingaleDiff[-1, ])  
  sdMartingaleDiff <- apply(MartingaleDiff[-1, ], 1, sd)  
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
                              start = 0, deltat = delta_Y),
                truemean = rep.int(0, dim(MartingaleDiff)[1]), 
                true_prices = Y[1:nrow(MartingaleDiff), 1], mc.prices = rowMeans(Dt))  
  start_Y <- start(Y)  
  
  mat.ci <- ts(mc.ci$conf.int, start = 0, deltat = delta_Y)  
  colnames(mat.ci) <- c("c.i lower bound", "c.i upper bound")  
  t_val_p_val <- ts(cbind(mc.ci$t, mc.ci$p.value), start = delta_Y, deltat = delta_Y)  
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

#'@title Correlation tests for the shocks
#'
#'@description
#'
#'This function performs correlation tests for the shocks generated by \code{\link{simshocks}}.
#'
#'@param x gaussian (bivariate) shocks, with correlation, generated by \code{\link{simshocks}}.
#' 
#'@param alternative indicates the alternative hypothesis and must be one of "two.sided", 
#'"greater" or "less". 
#'
#'@param method which correlation coefficient is to be used for the test : 
#'"pearson", "kendall", or "spearman".
#'
#'@param conf.level confidence level.
#'
#'@return a list with 2 components : estimated correlation coefficients, 
#'and confidence intervals for the estimated correlations. 
#'
#'@export
#'
#'@seealso \code{\link{esgplotbands}}
#'
#'@references
#'
#'D. J. Best & D. E. Roberts (1975), Algorithm AS 89: The Upper Tail 
#'Probabilities of Spearman's rho. Applied Statistics, 24, 377-379.
#'
#'Myles Hollander & Douglas A. Wolfe (1973), Nonparametric Statistical Methods.
#' New York: John Wiley & Sons. Pages 185-194 (Kendall and Spearman tests).
#'
#'@author Thierry Moudiki + stats package 
#'
#'@examples
#'
#'nb <- 500
#'
#'s0.par1 <- simshocks(n = nb, horizon = 3, frequency = "semi",
#'family = 1, par = 0.2)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 3, frequency = "semi", 
#'family = 1, par = 0.8)
#'
#'(test1 <- esgcortest(s0.par1))
#'(test2 <- esgcortest(s0.par2))
#'par(mfrow=c(2, 1))
#'esgplotbands(test1)
#'esgplotbands(test2)
#'
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
esgcortest <- cmpfun(esgcortest)

