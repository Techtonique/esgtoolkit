# Instantaneous forward rates ---------------------------------------------


#' Instantaneous forward rates
#' 
#' @description
#' Computes instantaneous forward rates from zero rates using various interpolation methods.
#' 
#' @param in.maturities Vector of input maturities
#' @param in.zerorates Vector of input zero rates
#' @param n Number of simulations
#' @param horizon Time horizon for forward rates
#' @param out.frequency Output frequency for forward rates. One of:
#'   \itemize{
#'     \item "annual" (default)
#'     \item "semi-annual"
#'     \item "quarterly" 
#'     \item "monthly"
#'     \item "weekly"
#'     \item "daily"
#'   }
#' @param method Interpolation method. One of:
#'   \itemize{
#'     \item "fmm" (default)
#'     \item "periodic"
#'     \item "natural"
#'     \item "monoH.FC"
#'     \item "hyman"
#'     \item "HCSPL"
#'     \item "SW"
#'   }
#' @param ... Additional arguments passed to interpolation functions
#' 
#' @details
#' The function computes instantaneous forward rates from zero rates using various interpolation methods.
#' It first converts zero rates to zero-coupon prices, then interpolates these prices using the specified method.
#' The forward rates are then computed from the interpolated prices.
#' 
#' The function supports different output frequencies and interpolation methods to suit various needs.
#' 
#' @return A time series object containing the instantaneous forward rates
#' 
#' @examples
#' # Generate sample data
#' maturities <- c(1, 2, 3, 5, 7, 10)
#' zero_rates <- c(0.01, 0.015, 0.02, 0.025, 0.03, 0.035)
#' 
#' # Compute forward rates with annual frequency
#' fwd_rates <- esgfwdrates(in.maturities = maturities,
#'                         in.zerorates = zero_rates,
#'                         n = 1000,
#'                         horizon = 10,
#'                         out.frequency = "annual",
#'                         method = "fmm")
#' 
#' @export
esgfwdrates <- function(in.maturities, in.zerorates,
                        n, horizon, 
                        out.frequency = c("annual", "semi-annual", 
                                          "quarterly", "monthly", 
                                          "weekly", "daily"),  
                        method = c("fmm", "periodic", "natural", 
                                   "monoH.FC", "hyman", "HCSPL", 
                                   "SW"),
                        ...)
{
  if(is.null(in.maturities) || is.null(in.zerorates))
    stop("Zero rates and maturities must be provided")
  
  if (floor(n) != n) 
    stop("'n' must be an integer")
  
  if (floor(horizon) != horizon) 
    stop("'horizon' must be an integer")
  
  max.in.maturities <- max(in.maturities)
  
  if (horizon > (max.in.maturities + .Machine$double.eps))
    stop("not enough maturities for the provided horizon")
  
  frequency <- match.arg(out.frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  method <- match.arg(method)
  
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
  
  # yc <- ycinter(matsin = in.maturities, matsout = tt, p = p, 
  #               typeres="prices", ...)  
  # ZC.prices <- fitted(yc)
  if (method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman"))
    ZC.prices <- stats::spline(x = in.maturities, y = p, xout = tt, 
                             method =  method, ...)$y
  
  if (base::identical(method, "HCSPL"))
    ZC.prices <- hermitecubicspline(p = p, 
                                    matsin = in.maturities, 
                                    matsout = tt, 
                                    typeres="prices")$values
  
  if (base::identical(method, "SW"))
    ZC.prices <- tZC_SW(p = p, 
                        u = in.maturities, 
                        t = tt, 
                        UFR = 0.0345, 
                        typeres="prices", 
                        ...)$values
    
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