#' Generate random SVJD scenarios with Feller-compliant parameters
#' 
#' @param n_series Number of scenarios (default=1000).
#' @param horizon Time steps per scenario (default=252).
#' @param freq Frequency ("daily", "weekly", "monthly").
#' @param type_return Type of return: basic return ("basic") or log-return ("log")
#' @param seed Random seed (optional).
#' @return List of scenarios (each with `log_returns`, `vol`, and `params`).
generate_svjd <- function(
    n_series = 100, horizon = 252, frequency="daily", 
    type_return = c("log", "basic"), 
    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  type_return <- match.arg(type_return)
  scenarios <- list()
  for (i in 1:n_series) {
    # Sample parameters with Feller condition enforced
    repeat {
      kappa <- runif(1, 1, 10)          # Mean reversion speed
      theta <- runif(1, 0.02, 0.06)      # Long-run volatility
      volvol <- runif(1, 0.1, 0.8)       # Volatility of volatility
      if (2 * kappa * theta >= volvol^2) break  # Feller check
    }
    # Generate scenario
    scenarios[[i]] <- rsvjd(
      n = 1, horizon = horizon, freq = freq,
      S0 = runif(1, 50, 150),          # Initial price (50-150)
      V0 = runif(1, 0.01, 0.09),       # Initial volatility (1%-9%)
      r0 = runif(1, 0.00, 0.05),       # Risk-free rate (0%-5%)
      kappa = kappa, theta = theta, volvol = volvol,
      rho = runif(1, -0.9, -0.3),      # Leverage effect (-0.9 to -0.3)
      lambda = runif(1, 0.01, 0.2),    # Jump intensity (1-20% annual)
      mu_J = runif(1, -0.1, 0.1),      # Jump mean (-10% to +10%)
      sigma_J = runif(1, 0.05, 0.3)    # Jump volatility (5%-30%)
    )$sim_price   
    scenarios[[i]] <- esgtoolkit::calculatereturns(scenarios[[i]], type=type_return)
  }
  return(unlist(scenarios))
}

#' Generate a single synthetic stock return path using a stochastic volatility model with jumps and regime switching.
#'
#' This function simulates returns based on a Heston-style stochastic volatility process,
#' optionally incorporating jumps (with configurable distributions), microstructure noise,
#' and regime switching via a Markov chain.
#'
#' @param n_days Number of trading days to simulate (default: 252 * 10).
#' @param mu Drift (mean return) per time step.
#' @param kappa Mean-reversion speed of the variance process.
#' @param theta Long-run average variance.
#' @param sigma_v Volatility of volatility (vol of variance process).
#' @param rho Correlation between volatility and returns (leverage effect).
#' @param lambda_jump Jump intensity (expected jumps per day).
#' @param jump_size_dist Distribution of jump sizes ("normal", "log_normal", or "exponential").
#' @param sigma_jump Scale parameter for jump size distribution.
#' @param noise_dist Type of microstructure noise ("normal" or "student_t").
#' @param noise_scale Standard deviation (or scale) of added microstructure noise.
#' @param noise_df Degrees of freedom for Student-t noise (if used).
#' @param regime_params List specifying regime switching behavior. Includes a transition matrix, and multipliers for `theta` and `kappa` during high-volatility regimes.
#' @param random_seed Optional random seed for reproducibility.
#'
#' @return A `data.table` with columns: `date`, `returns`, `variance`, `volatility`, `regime`.
#' @export
generate_synthetic_returns <- function(
    n_days = 252 * 10,
    mu = 0.0002,
    kappa = 0.05,
    theta = 0.0001,
    sigma_v = 0.01,
    rho = -0.7,
    lambda_jump = 0.05,
    jump_size_dist = "normal",
    sigma_jump = 0.02,
    noise_dist = "normal",
    noise_scale = 0.0005,
    noise_df = 5.0,
    regime_params = NULL,
    random_seed = NULL
) {
  
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # Default regime switching parameters
  if (is.null(regime_params)) {
    regime_params <- list(
      transition_matrix = matrix(c(0.99, 0.01, 0.03, 0.97), nrow = 2, byrow = TRUE),
      theta_high_multiplier = 3.0,
      kappa_high_multiplier = 2.0
    )
  }
  
  # Validate transition matrix
  if (!is.null(regime_params$transition_matrix)) {
    row_sums <- rowSums(regime_params$transition_matrix)
    if (!all(abs(row_sums - 1) < 1e-6)) {
      stop("Transition matrix rows must sum to 1.")
    }
  }
  
  # Initialize arrays
  n <- n_days
  v <- numeric(n)
  r <- numeric(n)
  regime <- integer(n)
  
  # Initialize starting values
  v[1] <- theta
  regime[1] <- 0
  
  # Pre-generate all random numbers
  z_vol <- rnorm(n)
  z_return <- rnorm(n)
  jump_indicators <- rpois(n, lambda = lambda_jump)
  
  # 1. Simulate the Markov chain for regimes
  for (t in 2:n) {
    prev_regime <- regime[t-1] + 1
    probs <- regime_params$transition_matrix[prev_regime, ]
    regime[t] <- sample(0:1, size = 1, prob = probs)
  }
  
  # 2. Simulate the Heston process with jumps
  for (t in 2:n) {
    current_regime <- regime[t]
    if (current_regime == 1) {
      theta_t <- theta * regime_params$theta_high_multiplier
      kappa_t <- kappa * regime_params$kappa_high_multiplier
    } else {
      theta_t <- theta
      kappa_t <- kappa
    }
    
    # Robust volatility discretization
    v_prev <- v[t-1]
    eta <- z_vol[t]
    
    drift <- kappa_t * (theta_t - max(v_prev, 0))
    volvol_term <- sigma_v * sqrt(max(v_prev, 0)) * eta
    v_new <- v_prev + drift + volvol_term
    v[t] <- max(v_new, 0)
    
    # Return process
    epsilon_t <- rho * eta + sqrt(1 - rho^2) * z_return[t]
    diffusion_component <- sqrt(max(v_prev, 0)) * epsilon_t
    
    # Jump process (single jump per period)
    J <- 0
    if (jump_indicators[t] > 0) {
      if (jump_size_dist == "normal") {
        J <- rnorm(1, mean = 0, sd = sigma_jump)
      } else if (jump_size_dist == "log_normal") {
        log_J <- rnorm(1, mean = -0.5 * sigma_jump^2, sd = sigma_jump)
        J <- exp(log_J) - 1
      } else if (jump_size_dist == "exponential") {
        sign <- sample(c(-1, 1), 1)
        J <- sign * rexp(1, rate = 1/sigma_jump)
      } else {
        stop("Invalid jump_size_dist.")
      }
    }
    
    r[t] <- mu + diffusion_component + J
  }
  
  # 3. Add microstructure noise
  if (noise_dist == "normal") {
    noise <- rnorm(n, mean = 0, sd = noise_scale)
  } else if (noise_dist == "student_t") {
    if (noise_df <= 2) stop("noise_df must be > 2")
    scale_factor <- noise_scale / sqrt(noise_df / (noise_df - 2))
    noise <- rt(n, df = noise_df) * scale_factor
  } else {
    stop("Invalid noise_dist.")
  }
  r <- r + noise
  
  # 4. Create output data.table
  dt <- data.table(
    date = seq.Date(as.Date("1970-01-01"), by = "day", length.out = n),
    returns = r,
    variance = v,
    volatility = sqrt(v),
    regime = factor(regime, levels = c(0, 1), labels = c("Low Vol", "High Vol"))
  )
  
  return(dt)
}

#' Generate a list of synthetic return paths using stochastic volatility models with randomized parameters.
#'
#' This function produces a collection of synthetic return series, each simulated with randomized parameters
#' for drift, volatility, jump characteristics, and regime-switching behavior. Useful for pretraining models
#' or scenario testing.
#'
#' @param n_paths Number of return paths to generate.
#' @param horizon Length of each return path in trading days.
#' @param frequency Frequency of the data (e.g., "daily").
#' @param jump_type Type of jump size distribution: "normal", "log_normal", "exponential", or "mixed".
#' @param include_regime_switching Logical; whether to include Markov regime switching.
#' @param random_seed Optional seed for reproducibility.
#'
#' @return A list of class `diverse_sv_paths` containing:
#' - `paths`: A list of numeric vectors (returns)
#' - `parameters`: The parameter sets used to generate each path
#' - `horizon`: The path length
#' - `frequency`: The data frequency
#' - `n_paths`: Number of paths
#' - `generation_date`: Timestamp of generation
#'
#' @export
generate_diverse_sv_paths <- function(
    n_paths = 10000,
    horizon = 252 * 5,
    frequency = "daily",
    jump_type = "mixed",
    include_regime_switching = TRUE,
    random_seed = NULL
) {
  
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  all_paths <- vector("list", n_paths)
  all_params <- vector("list", n_paths)
  
  for (i in 1:n_paths) {
    
    params <- list()
    
    # Sample diverse parameters
    params$mu <- runif(1, -0.0005, 0.0005)
    params$kappa <- runif(1, 0.02, 0.15)
    params$theta <- runif(1, 5e-5, 3e-4)
    params$sigma_v <- runif(1, 0.005, 0.03)
    params$rho <- runif(1, -0.85, -0.4)
    
    # Jump parameters
    params$lambda_jump <- sample(c(
      runif(1, 0.005, 0.02),
      runif(1, 0.02, 0.08),
      runif(1, 0.08, 0.15)
    ), 1)
    
    if (jump_type == "mixed") {
      params$jump_size_dist <- sample(
        c("normal", "log_normal", "exponential"), 
        1,
        prob = c(0.4, 0.4, 0.2)
      )
    } else {
      params$jump_size_dist <- jump_type
    }
    
    params$sigma_jump <- runif(1, 0.01, 0.05)
    params$noise_dist <- sample(c("normal", "student_t"), 1, prob = c(0.7, 0.3))
    params$noise_scale <- runif(1, 1e-5, 2e-4)
    params$noise_df <- runif(1, 3, 8)
    
    # Regime switching parameters with FIXED transition matrices
    if (include_regime_switching) {
      regime_type <- sample(1:3, 1)
      
      if (regime_type == 1) {
        # Persistent regimes - FIXED: ensure rows sum to 1
        p11 <- runif(1, 0.97, 0.995)
        p12 <- 1 - p11
        p21 <- runif(1, 0.01, 0.05)
        p22 <- 1 - p21
        transition_matrix <- matrix(c(p11, p12, p21, p22), nrow = 2, byrow = TRUE)
      } else if (regime_type == 2) {
        # Mean-reverting regimes - FIXED: ensure rows sum to 1
        p11 <- runif(1, 0.85, 0.92)
        p12 <- 1 - p11
        p21 <- runif(1, 0.08, 0.15)
        p22 <- 1 - p21
        transition_matrix <- matrix(c(p11, p12, p21, p22), nrow = 2, byrow = TRUE)
      } else {
        # Rapid switching regimes - FIXED: ensure rows sum to 1
        p11 <- runif(1, 0.7, 0.8)
        p12 <- 1 - p11
        p21 <- runif(1, 0.2, 0.3)
        p22 <- 1 - p21
        transition_matrix <- matrix(c(p11, p12, p21, p22), nrow = 2, byrow = TRUE)
      }
      
      params$regime_params <- list(
        transition_matrix = transition_matrix,
        theta_high_multiplier = runif(1, 2.0, 5.0),
        kappa_high_multiplier = runif(1, 1.5, 3.0)
      )
    } else {
      params$regime_params <- NULL
    }
    
    # Generate the path
    path_data <- generate_synthetic_returns(
      n_days = horizon,
      mu = params$mu,
      kappa = params$kappa,
      theta = params$theta,
      sigma_v = params$sigma_v,
      rho = params$rho,
      lambda_jump = params$lambda_jump,
      jump_size_dist = params$jump_size_dist,
      sigma_jump = params$sigma_jump,
      noise_dist = params$noise_dist,
      noise_scale = params$noise_scale,
      noise_df = params$noise_df,
      regime_params = params$regime_params
    )
    
    all_paths[[i]] <- path_data$returns
    all_params[[i]] <- params
    
    if (i %% 1000 == 0) {
      message(sprintf("Generated %d/%d paths", i, n_paths))
    }
  }
  
  result <- list(
    paths = all_paths,
    parameters = all_params,
    horizon = horizon,
    frequency = frequency,
    n_paths = n_paths,
    generation_date = Sys.time()
  )
  
  class(result) <- "diverse_sv_paths"
  
  return(result)
}
