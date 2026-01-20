#' Simulate Zero-Coupon Bond Prices using G2++ Model
#'
#' @param n Number of scenarios to simulate.
#' @param horizon Time steps for simulation (e.g., 5 for 5 years).
#' @param freq Frequency of simulation (default is "semi-annual").
#' @param u Observed maturities (vector).
#' @param txZC Yield to maturities (vector, same length as u).
#' @param a G2++ mean reversion for factor x.
#' @param b G2++ mean reversion for factor y.
#' @param sigma Volatility of factor x.
#' @param eta Volatility of factor y.
#' @param rho Correlation between factors.
#' @param maturities Bond maturity times.
#' @param start_ Starting time for ts object.
#' @param seed Random seed for reproducibility.
#' @param ... Additional parameters to be passed to simdiff or simshocks
#' @return A time series of simulated bond prices for each scenario.
rpriceg2plus <- function(n = 10L,
                         horizon = 5L,
                         freq = "semi-annual",
                         u = 1:30,
                         txZC = c(0.01422, 0.01309, 0.01380, 0.01549, 0.01747, 0.01940, 
                                  0.02104, 0.02236, 0.02348, 0.02446, 0.02535, 0.02614, 
                                  0.02679, 0.02727, 0.02760, 0.02779, 0.02787, 0.02786, 
                                  0.02776, 0.02762, 0.02745, 0.02727, 0.02707, 0.02686, 
                                  0.02663, 0.02640, 0.02618, 0.02597, 0.02578, 0.02563),
                         a = 0.5,
                         b = 0.3541203,
                         sigma = 0.09416266,
                         eta = 0.08439934,
                         rho = -0.99855687,
                         maturities = c(5, 7, 10),
                         start_ = c(2000,1),
                         seed = NULL, 
                         ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Parameter validation
  stopifnot(sigma > 0, eta > 0, a > 0, b > 0)
  stopifnot(rho >= -1, rho <= 1)
  stopifnot(length(u) == length(txZC))
  
  # Frequency handling
  freq_map <- list(
    annual = list(delta_t = 1, tsfreq = 1),
    'semi-annual' = list(delta_t = 1/2, tsfreq = 2),
    quarterly = list(delta_t = 1/4, tsfreq = 4),
    monthly = list(delta_t = 1/12, tsfreq = 12),
    weekly = list(delta_t = 1/52, tsfreq = 52),
    daily = list(delta_t = 1/252, tsfreq = 252)
  )
  
  if (is.character(freq)) {
    fkey <- tolower(freq)
    if (!fkey %in% names(freq_map)) stop("Unknown freq")
    delta_t <- freq_map[[fkey]]$delta_t
    ts_frequency <- freq_map[[fkey]]$tsfreq
  } else if (is.numeric(freq) && freq > 0) {
    delta_t <- 1 / freq
    ts_frequency <- as.integer(freq)
  } else stop("freq must be char or positive numeric")
  
  # Market curve setup
  u_extended <- c(0, u)
  txZC_extended <- c(0, txZC)
  P_M_0 <- exp(-u_extended * txZC_extended)
  spline_fit <- splinefun(u_extended, P_M_0, method = "hyman")
  
  # Time grid - CRITICAL FIX: Include final maturity time
  t_out <- seq(0, horizon, by = delta_t)
  n_time <- length(t_out)
  n_bonds <- length(maturities)
  
  # Simulate OU factors
  eps <- esgtoolkit::simshocks(n = n, horizon = horizon, frequency = freq, 
                               family = 1, par = rho, start_ = start_, 
                               ...)
  x <- esgtoolkit::simdiff(n = n, horizon = horizon, frequency = freq, model = "OU",
                           x0 = 0, theta1 = 0, theta2 = a, theta3 = sigma, 
                           eps = eps[[1]], start_ = start_, 
                           ...)
  y <- esgtoolkit::simdiff(n = n, horizon = horizon, frequency = freq, model = "OU",
                           x0 = 0, theta1 = 0, theta2 = b, theta3 = eta, 
                           eps = eps[[2]], start_ = start_, 
                           ...)
  
  # Ensure proper dimensions
  if (nrow(x) != n_time) { 
    x <- t(x) 
    y <- t(y) 
  }
  
  # CRITICAL FIX: Proper time-to-maturity calculation
  # Create time x scenarios x maturities arrays
  time_array <- array(rep(t_out, times = n * n_bonds), 
                      dim = c(n_time, n, n_bonds))
  maturity_array <- array(rep(maturities, each = n_time * n), 
                          dim = c(n_time, n, n_bonds))
  
  # Time to maturity: max(0, T - t) to avoid negative values
  tau_array <- pmax(maturity_array - time_array, 0)
  
  # M(t,T) - Lemma 4.2.1 Equation (4.9)
  coef_x <- (1 - exp(-a * tau_array)) / a
  coef_y <- (1 - exp(-b * tau_array)) / b
  x3d <- array(rep(x, times = n_bonds), dim = c(n_time, n, n_bonds))
  y3d <- array(rep(y, times = n_bonds), dim = c(n_time, n, n_bonds))
  M3d <- coef_x * x3d + coef_y * y3d
  
  # V(t,T) - Lemma 4.2.1 Equation (4.10)
  V3d <- calculate_V(tau_array, a, b, sigma, eta, rho)
  
  # Calculate V(0,t) and V(0,T) for integral term
  V_0_t <- calculate_V(t_out, a, b, sigma, eta, rho)
  V_0_maturity <- calculate_V(maturities, a, b, sigma, eta, rho)
  
  V_0_t_array <- array(rep(V_0_t, times = n * n_bonds), 
                       dim = c(n_time, n, n_bonds))
  V_0_maturity_array <- array(rep(V_0_maturity, each = n_time * n), 
                              dim = c(n_time, n, n_bonds))
  
  # Integral term - Theorem 4.2.1 + Corollary 4.2.1 Equation (4.13)
  P0_t <- spline_fit(t_out)
  P0_t_array <- array(rep(P0_t, times = n_bonds), 
                      dim = c(n_time, n, n_bonds))
  P0_mat_array <- array(rep(spline_fit(maturities), each = n_time), 
                        dim = c(n_time, n, n_bonds))
  
  integral_array <- -log(P0_mat_array / P0_t_array) - 0.5 * (V_0_maturity_array - V_0_t_array)
  
  # Bond prices - Theorem 4.2.1 Equation (4.11)
  bond_prices <- exp(-integral_array - M3d - 0.5 * V3d)
  
  # CRITICAL FIX: Only set to 1 when ACTUALLY at maturity (tau = 0)
  # Use a small tolerance for numerical stability
  bond_prices[tau_array <= 1e-10] <- 1.0
  
  # Additional arbitrage check: enforce P(t,T1) >= P(t,T2) for T1 < T2
  bond_prices <- enforce_monotonicity(bond_prices, maturities)
  
  # Return list of time series objects
  bond_ts_list <- lapply(seq_len(n_bonds), function(b) {
    ts(bond_prices[,,b], start = start_, frequency = ts_frequency)
  })
  names(bond_ts_list) <- paste0("Maturity_", maturities)
  
  return(bond_ts_list)
}

#' Enforce monotonicity: P(t,T1) >= P(t,T2) for T1 < T2
enforce_monotonicity <- function(bond_prices, maturities) {
  n_time <- dim(bond_prices)[1]
  n_scenarios <- dim(bond_prices)[2]
  n_bonds <- dim(bond_prices)[3]
  
  # Sort maturities and corresponding price arrays
  maturity_order <- order(maturities)
  sorted_maturities <- maturities[maturity_order]
  
  for (t in 1:n_time) {
    for (s in 1:n_scenarios) {
      # Get prices for this time/scenario across maturities
      prices <- bond_prices[t, s, maturity_order]
      
      # Enforce monotonicity: longer maturities should have lower prices
      for (i in 2:n_bonds) {
        if (prices[i] > prices[i-1]) {
          # If longer maturity has higher price, cap it at previous maturity's price
          prices[i] <- min(prices[i], prices[i-1])
        }
      }
      
      # Update the prices (maintain original order)
      bond_prices[t, s, maturity_order] <- prices
    }
  }
  
  return(bond_prices)
}

#' Calculate variance term V(t,T) for G2++ model
calculate_V <- function(tau, a, b, sigma, eta, rho) {
  # Handle both scalar and array inputs
  term1 <- (sigma^2 / a^2) * (
    tau + (2/a) * exp(-a * tau) - 
      (1/(2*a)) * exp(-2*a * tau) - (3/(2*a))
  )
  
  term2 <- (eta^2 / b^2) * (
    tau + (2/b) * exp(-b * tau) - 
      (1/(2*b)) * exp(-2*b * tau) - (3/(2*b))
  )
  
  term3 <- 2 * rho * (sigma * eta) / (a * b) * (
    tau + 
      (exp(-a * tau) - 1)/a +
      (exp(-b * tau) - 1)/b -
      (exp(-(a + b) * tau) - 1)/(a + b)
  )
  
  return(term1 + term2 + term3)
}