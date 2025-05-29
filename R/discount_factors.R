#' Stochastic Discount Factors
#' 
#' @description
#' Computes stochastic discount factors or discounted values based on interest rates and cash flows.
#' 
#' @param r A numeric vector or time series of interest rates
#' @param X A numeric vector or time series of cash flows or values to be discounted
#' 
#' @details
#' This function handles various combinations of scalar and time series inputs for both interest rates and cash flows.
#' For scalar inputs, it performs simple exponential discounting.
#' For time series inputs, it computes cumulative discounting over time.
#' 
#' @return A time series of discounted values
#' 
#' @examples
#' # Simple scalar discounting
#' esgdiscountfactor(0.05, 100)
#' 
#' # Time series discounting
#' #r <- ts(rep(0.05, 10), start = 0, deltat = 1)
#' #X <- ts(rep(100, 10), start = 0, deltat = 1)
#' #esgdiscountfactor(r, X)
#' 
#' @export
esgdiscountfactor <- function(r, X)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  length_r <- length(r)
  length_X <- length(X)
  start_r <- start(r)
  deltat.r <- deltat(r)
  start_X <- start(X)
  deltat_X <- deltat(X)
  
  if(length_r == 1 && length_X == 1) 
  {
    return(X*exp(-r)) 
  }
  
  if(length_r == 1 && length_X != 1) 
  {
    r <- ts(matrix(r, nrow(X), ncol(X)), 
            start = start_X, deltat = deltat_X)
    
    if(tsp(X)[1] > 0)
    {
      return(ts(X*exp(-apply(r, 2, cumsum)*deltat_X), 
                start = 0, 
                deltat = deltat_X))      
    }
    else
    {
      Int_r <- exp(-apply(r, 2, cumsum)*deltat_X)
      return(ts(X*rbind(rep(1, ncol(X)), 
                        Int_r[1:(nrow(X)-1), ]),                                            
                start = 0, 
                deltat = deltat_X))
    }
  } 
  
  if(length_r != 1 && length_X == 1)
  {
    X <- ts(matrix(X, nrow(r), ncol(r)), 
            start = start_r, deltat = deltat.r)
    
    return(ts(X*exp(-apply(r, 2, cumsum)*deltat.r),                                          
              start = 0, 
              deltat = deltat.r))
  }
  
  if(length_r != 1 && length_X != 1)
  {
    if (all.equal(start(X), start(r)) && all.equal(frequency(X), frequency(r)))
    {
      return(X*exp(-apply(r, 2, cumsum)*deltat(r)))
    }
    
    if(tsp(X)[1] > 0)
    {
      return(suppressWarnings(ts(X*window(exp(-apply(r, 2, cumsum)*deltat.r), 
                                          start = start_X, 
                                          deltat = deltat_X), 
                                 start = 0, 
                                 deltat = deltat_X)))
    }
    else
    {
      Int_r <- ts(exp(-apply(r, 2, cumsum)*deltat.r), deltat = deltat.r)
      return(suppressWarnings(ts(X*rbind(rep(1, ncol(X)), 
                                         Int_r[1:(nrow(X)-1), ]),
                                 start = 0, 
                                 deltat = deltat_X)))
    }
  }
}


#' Estimation of discounted asset prices
#' 
#' @description
#' Calculates the discounted asset prices using simulated interest rates and asset values.
#' 
#' @param r Interest rate time series or constant value
#' @param X Asset price time series or constant value
#' @param maturity Optional maturity date for the calculation
#' 
#' @details
#' The function takes interest rates and asset prices as inputs and calculates the discounted
#' asset prices. It handles both time series and constant values for both inputs.
#' 
#' If a maturity date is specified, the function returns the discounted price at that maturity.
#' Otherwise, it returns the time series of discounted prices averaged across all scenarios.
#' 
#' @return A time series object containing the discounted asset prices, or a single value
#'         if a maturity date is specified
#' 
#' @examples
#' # Using constant values
#' r <- 0.05
#' X <- 100
#' price <- esgmcprices(r, X)
#' 
#' # Using time series
#' r_ts <- ts(matrix(rnorm(100), 10, 10), start = 0, deltat = 0.1)
#' X_ts <- ts(matrix(rnorm(100), 10, 10), start = 0, deltat = 0.1)
#' prices_ts <- esgmcprices(r_ts, X_ts)
#' 
#' # With maturity
#' price_at_maturity <- esgmcprices(r_ts, X_ts, maturity = 1)
#' 
#' @export
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
  
  Z <- ts(rowMeans(Y), 
          start = start(Y), 
          deltat = deltat(Y))
  
  if(!is.null(maturity))
  {
    return(window(Z, start = maturity.out, end = maturity.out))
  }
  else  
  {
    return(Z)
  }
  
}

