#' Simulate SVJD model (with Feller condition check) with G2++ short rates model
#'
#' @param n Number of simulations.
#' @param horizon Time steps (e.g., 252 for daily data).
#' @param freq Frequency ("daily", "weekly", "monthly").
#' @param S0 Initial stock price (default=100).
#' @param rho_stock_vol Correlation stock price / stochastic volatility
#' @param rho_g2plus Correlation between G2++ factors
#' @param ... Additional parameters to be passed to \code{\link{rsvjd}} and
#' \code{\link{rg2plus}}
rsvjdg2plus <- function(n = 10,
                        horizon = 5,
                        freq = "quarterly",
                        S0 = 100,
                        rho_stock_vol = -0.7,
                        rho_g2plus = -0.99855687,
                        ...)
{
  short_rates <- esgtoolkit::rg2plus(
    n = n,
    freq = freq,
    horizon = horizon,
    rho = rho_g2plus,
    ...
  )
  stock_prices <- esgtoolkit::rsvjd(S0 = S0,
    r0 = short_rates,
    n = n,
    freq = freq,
    horizon = horizon, 
    rho = rho_stock_vol,
    ...
  )
  return(list(short_rates = short_rates, 
              stock_prices = stock_prices))
}
