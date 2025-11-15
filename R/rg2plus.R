#' Simulate G2++ short rates model
#'
#' @param n Number of scenarios to simulate.
#' @param horizon Time steps (semi-annual, default = 20).
#' @param freq Frequency of simulation (default is "semi-annual").
#' @param u Observed maturities (vector).
#' @param txZC Yield to maturities (vector, same length as u).
#' @param a G2++ mean reversion for factor x.
#' @param b G2++ mean reversion for factor y.
#' @param sigma Volatility of factor x.
#' @param eta Volatility of factor y.
#' @param rho Correlation between factors.
#' @param methodyc Interpolation method for forward rates ("fmm", "hyman", "HCSPL", "SW").
#' @return A time series of simulated short rates for each scenario.
rg2plus <- function(n = 10L,
                    horizon = 5L,
                    freq = "semi-annual",
                    u = 1:30,
                    txZC = c(
                      0.01422,
                      0.01309,
                      0.01380,
                      0.01549,
                      0.01747,
                      0.01940,
                      0.02104,
                      0.02236,
                      0.02348,
                      0.02446,
                      0.02535,
                      0.02614,
                      0.02679,
                      0.02727,
                      0.02760,
                      0.02779,
                      0.02787,
                      0.02786,
                      0.02776,
                      0.02762,
                      0.02745,
                      0.02727,
                      0.02707,
                      0.02686,
                      0.02663,
                      0.02640,
                      0.02618,
                      0.02597,
                      0.02578,
                      0.02563
                    ),
                    a = 0.5,
                    b = 0.3541203,
                    sigma = 0.09416266,
                    eta = 0.08439934,
                    rho = -0.99855687,
                    methodyc = c("fmm", "hyman", "HCSPL", "SW")) {
  # Compute delta_t based on freq
  if (is.character(freq)) {
    delta_t <- switch(
      tolower(freq),
      "annual" = 1,
      "semi-annual" = 1 / 2,
      "quarterly" = 1 / 4,
      "monthly" = 1 / 12,
      "weekly" = 1 / 52,
      "daily" = 1 / 252,
      stop("Unknown frequency string")
    )
  } else if (is.numeric(freq) && freq > 0) {
    delta_t <- 1 / freq
  } else {
    stop("freq must be a character string or a positive numeric frequency")
  }
  
  # Generate correlated shocks
  eps <- esgtoolkit::simshocks(
    n = n,
    horizon = horizon,
    frequency = freq,
    family = 1,
    par = rho
  )
  
  # Simulate OU factors
  x <- esgtoolkit::simdiff(
    n = n,
    horizon = horizon,
    frequency = freq,
    model = "OU",
    x0 = 0,
    theta1 = 0,
    theta2 = a,
    theta3 = sigma,
    eps = eps[[1]]
  )
  y <- esgtoolkit::simdiff(
    n = n,
    horizon = horizon,
    frequency = freq,
    model = "OU",
    x0 = 0,
    theta1 = 0,
    theta2 = b,
    theta3 = eta,
    eps = eps[[2]]
  )
  
  # Forward rate curve
  methodyc <- match.arg(methodyc)
  fwdrates <- esgtoolkit::esgfwdrates(
    n = n,
    horizon = horizon,
    out.frequency = freq,
    in.maturities = u,
    in.zerorates = txZC,
    method = methodyc
  )
  fwdrates <- window(fwdrates, end = horizon)
  
  # Compute phi adjustment
  t.out <- seq(0, horizon, by = delta_t)
  param.phi <- 0.5 * (sigma^2) * (1 - exp(-a * t.out))^2 / a^2 +
    0.5 * (eta^2) * (1 - exp(-b * t.out))^2 / b^2 +
    rho * sigma * eta * (1 - exp(-a * t.out)) * (1 - exp(-b * t.out)) / (a * b)
  param.phi <- ts(replicate(n, param.phi),
                  start = start(x),
                  deltat = delta_t)
  
  # Short rates
  r <- x + y + fwdrates + param.phi
  colnames(r) <- paste0("Series ", 1:n)
  
  return(r)
}
