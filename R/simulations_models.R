#' Simulate SVJD model using esgtoolkit (with Feller condition check)
#' 
#' @param n Number of simulations.
#' @param horizon Time steps (e.g., 252 for daily data).
#' @param freq Frequency ("daily", "weekly", "monthly").
#' @param S0 Initial price (default=100).
#' @param V0 Initial volatility (default=0.04).
#' @param r0 Risk-free rate (default=0.02).
#' @param kappa Mean reversion speed (must satisfy Feller: 2*kappa*theta >= volvol^2).
#' @param theta Long-run mean volatility.
#' @param volvol Volatility of volatility (must satisfy Feller).
#' @param rho Leverage effect (correlation between shocks).
#' @param lambda Jump intensity.
#' @param mu_J Mean jump size.
#' @param sigma_J Std dev of jump size.
#' @return List of simulated `price` and `vol`.
rsvjd <- function(
    n = 10, horizon = 1, freq = "daily",
    S0 = 100, V0 = 0.04, r0 = 0.02,
    kappa = 5, theta = 0.04, volvol = 0.5, rho = -0.7,
    lambda = 0.1, mu_J = 0.0, sigma_J = 0.1) { 
  
  # Check Feller condition
  if (2 * kappa * theta < volvol^2) {
    warning(
      "Feller condition violated (2*kappa*theta >= volvol^2 not satisfied). ",
      "Volatility may hit zero. Adjust parameters."
    )
  }  

  # Simulate shocks and paths
  shocks <- esgtoolkit::simshocks(
    n = n, horizon = horizon, frequency = freq,
    method = "anti", family = 1, par = rho
  )
  
  sim_vol <- esgtoolkit::simdiff(
    n = n, horizon = horizon, frequency = freq,
    model = "CIR", x0 = V0,
    theta1 = kappa * theta, theta2 = kappa, theta3 = volvol,
    eps = shocks[[1]]
  )
  
  sim_price <- esgtoolkit::simdiff(
    n = n, horizon = horizon, frequency = freq,
    model = "GBM", x0 = S0,
    theta1 = r0 - lambda * mu_J, theta2 = sim_vol,
    lambda = lambda, mu_z = mu_J, sigma_z = sigma_J,
    eps = shocks[[2]]
  )
  
  return(list(price = sim_price, vol = sim_vol))
}