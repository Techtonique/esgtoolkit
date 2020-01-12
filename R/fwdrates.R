# Instantaneous forward rates ---------------------------------------------


#'@title 
#'
#'Instantaneous forward rates
#'
#'@description
#'
#'This function provides instantaneous forward rates. They can be used 
#'in no-arbitrage short rate models, to fit the yield curve exactly. 
#'
#'@param in.maturities input maturities
#'
#'@param in.zerorates input zero rates
#'
#'@param n number of independent observations
#'
#'@param horizon horizon of projection
#'
#'@param out.frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252)
#'
#'@param ... additional parameters provided to \code{\link{ycinter}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@references
#'
#' Thierry Moudiki (2013). ycinterextra: Yield curve or zero-coupon prices interpolation and extrapolation. R package version 0.1. URL 
#'\url{https://CRAN.R-project.org/package=ycinterextra}
#'
#'@examples
#'
#'# Yield to maturities
#'txZC <- c(0.01422,0.01309,0.01380,0.01549,0.01747,0.01940,0.02104,0.02236,0.02348,
#'          0.02446,0.02535,0.02614,0.02679,0.02727,0.02760,0.02779,0.02787,0.02786,0.02776
#'          ,0.02762,0.02745,0.02727,0.02707,0.02686,0.02663,0.02640,0.02618,0.02597,0.02578,0.02563)
#'
#'# Observed maturities
#'u <- 1:30
#'
#'\dontrun{
#'par(mfrow=c(2,2))
#'fwdNS <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "NS")
#'matplot(time(fwdNS), fwdNS, type = 'l')
#'
#'fwdSV <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "SV")
#'matplot(time(fwdSV), fwdSV, type = 'l')
#'
#'fwdSW <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "SW")
#'matplot(time(fwdSW), fwdSW, type = 'l')
#'
#'fwdHCSPL <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                        n = 10, horizon = 20, 
#'                        out.frequency = "semi-annual", method = "HCSPL")
#'matplot(time(fwdHCSPL), fwdHCSPL, type = 'l')
#'}
esgfwdrates <- function(in.maturities, in.zerorates,
                        n, horizon, 
                        out.frequency = c("annual", "semi-annual", 
                                          "quarterly", "monthly", 
                                          "weekly", "daily"),  
                        ...)
{
  if(is.null(in.maturities) || is.null(in.zerorates))
    stop("Zero rates and maturities must be provided")
  
  if (floor(n) != n) 
    stop("'n' must be an integer")
  
  if (floor(horizon) != horizon) 
    stop("'horizon' must be an integer")
  
  max.in.maturities <- max(in.maturities)
  
  if (horizon > max.in.maturities)
    stop("not enough maturities for the provided horizon")
  
  frequency <- match.arg(out.frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  
  # simply coumpounded zero-coupon price
  pricefromeuribor <- function(t, T, L)
  {
    # Brigo P. 7
    tau <- T-t
    return(as.vector(1/(1 + L*tau)))    
  } 
  
  p <- pricefromeuribor(t = 0, T = in.maturities, L = in.zerorates)
  tt <- seq(from = min(in.maturities), to = max.in.maturities, by = delta)
  nn <- length(tt)
  
  yc <- ycinter(matsin = in.maturities, matsout = tt, p = p, 
                typeres="prices", ...)  
  
  ZC.prices <- fitted(yc)
  
  ZC.prices <- c(1 + ((p[1] - 1)/(in.maturities[1] - 0))*(seq(0, 1, by = delta) - 0), 
                 ZC.prices[-1])
  
  tt <- seq(0, max(in.maturities), by = delta)
  
  nn <- length(tt)
  
  fwdrates <- ts(data = (1/delta)*(ZC.prices[-nn]/ZC.prices[-1]-1), 
                 start = 0,
                 deltat = delta)     
  
  return(ts(replicate(n, fwdrates), start = 0,
            deltat = delta))
}
esgfwdrates <- cmpfun(esgfwdrates)