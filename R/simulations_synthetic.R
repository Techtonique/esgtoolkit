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


#' Enhanced R-vine copula simulation with proper distribution preservation
#'
#' Simulates multivariate time series data using R-vine copulas while preserving
#' the dependence structure and marginal distributions of the original data.
#'
#' @param data Matrix or data.frame of multivariate time series
#' @param n Number of simulations to generate
#' @param seed Random seed for reproducibility
#' @param verbose Whether to print fitting information
#' @param n_trials Number of trials to select best fit
#' @param oversample_factor Factor to oversample for better selection (default 1.5)
#' @param score_weights Vector of weights for scoring criteria. Must be length 5 and sum to 1.
#'                      Order: [Kendall_cor, Pearson_cor, KS_test, mean_error, sd_error]
#'
#' @return A list of class 'rvine_simulation' containing:
#' \itemize{
#'   \item \code{simulated_data} - The simulated multivariate time series
#'   \item \code{diagnostics} - Comprehensive diagnostic information including:
#'   \itemize{
#'     \item Correlation matrices and errors
#'     \item Distribution comparison metrics
#'     \item Quality scores and trial results
#'     \item Model information
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default weights
#' result <- simulate_rvine(uschange, n = 200)
#' 
#' # Custom weights emphasizing correlation preservation
#' result <- simulate_rvine(uschange, n = 200, 
#'                         score_weights = c(0.6, 0.2, 0.1, 0.05, 0.05))
#' 
#' # Plot results
#' plot(result)
#' }
#'
#' @export
simulate_rvine <- function(data, n = 100, seed = 123, 
                           verbose = FALSE, n_trials = 10,
                           oversample_factor = 1.5,
                           score_weights = c(0.4, 0.2, 0.2, 0.1, 0.1)) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data.frame")
  }
  
  # Validate score weights
  if (length(score_weights) != 5) {
    stop("score_weights must be a vector of length 5")
  }
  if (abs(sum(score_weights) - 1) > 1e-10) {
    stop("score_weights must sum to 1")
  }
  if (any(score_weights < 0)) {
    stop("score_weights must be non-negative")
  }
  
  data <- as.matrix(data)  # Ensure matrix format
  n_obs <- nrow(data)
  n_series <- ncol(data)
  series_names <- colnames(data)
  
  # Enhanced validation
  if (n_obs < 30) {
    stop("Sample size too small for reliable copula estimation (minimum 30 observations)")
  } else if (n_obs < 50) {
    warning("Small sample size may affect copula estimation quality")
  }
  
  if (any(is.na(data))) {
    stop("Data contains missing values. Please handle missing data before simulation.")
  }
  
  if (n_series < 2) {
    stop("At least 2 variables required for copula modeling")
  }
  
  # Store original statistics for comparison
  original_cor <- cor(data, method = "kendall")  # Use Kendall's tau for copulas
  original_cor_pearson <- cor(data)
  original_means <- colMeans(data)
  original_sds <- apply(data, 2, sd)
  
  # Improved uniformization with better boundary handling
  if (verbose) message("Transforming data to uniform margins with improved boundary handling...")
  U <- matrix(NA, nrow = n_obs, ncol = n_series)
  for (i in 1:n_series) {
    # Use empirical CDF with interpolation for smoother transformation
    ecdf_func <- ecdf(data[, i])
    u_vals <- ecdf_func(data[, i])
    # Adjust for boundary issues
    u_vals[u_vals == 0] <- 1 / (2 * n_obs)
    u_vals[u_vals == 1] <- 1 - 1 / (2 * n_obs)
    U[, i] <- u_vals
  }
  
  # Fit the R-vine copula model with enhanced family selection
  if (verbose) message("Fitting R-vine copula model...")
  # Extended set of valid copula families
  valid_families <- c(1, 2, 3, 4, 5, 6,           # Basic families
                      13, 14, 16, 17, 18, 19, 20,  # Rotated Clayton
                      23, 24, 26, 27, 28, 29, 30,  # Rotated Gumbel  
                      33, 34, 36, 37, 38, 39, 40)  # Rotated Frank
  
  # Try different selection criteria and choose best
  tryCatch({
    RVM_U <- VineCopula::RVineStructureSelect(
      data = U, 
      familyset = valid_families, 
      type = 0, 
      selectioncrit = "BIC",
      progress = verbose,
      treecrit = "tau",
      trunclevel = min(n_series - 1, 3)  # Limit complexity for small datasets
    )
  }, error = function(e) {
    if (verbose) message("Full model failed, trying simplified approach...")
    # Fallback to basic families
    RVM_U <- VineCopula::RVineStructureSelect(
      data = U, 
      familyset = c(1:6), 
      type = 0, 
      selectioncrit = "AIC",
      progress = verbose
    )
  })
  
  if (verbose) {
    message("R-vine copula model fitted successfully")
    print(summary(RVM_U))
  }
  
  # Multiple simulation trials with improved back-transformation
  best_sim <- NULL
  best_score <- Inf
  trial_scores <- numeric(n_trials)
  
  if (verbose) message(sprintf("Running %d simulation trials...", n_trials))
  
  for (trial in 1:n_trials) {
    # Better seed management
    set.seed(seed * 1000 + trial * 17)  # More distributed seeds
    # Simulate with oversampling for better selection
    n_sim <- max(n, ceiling(n * oversample_factor))
    tryCatch({
      rvine_simulation <- VineCopula::RVineSim(N = n_sim, RVM = RVM_U)
      # Select best n observations if oversampling
      if (n_sim > n) {
        # Calculate quality for each potential subsample
        subsample_scores <- numeric(n_sim - n + 1)
        for (start_idx in 1:(n_sim - n + 1)) {
          end_idx <- start_idx + n - 1
          temp_sim <- rvine_simulation[start_idx:end_idx, , drop = FALSE]
          temp_cor <- cor(temp_sim, method = "kendall")
          subsample_scores[start_idx] <- mean(abs(temp_cor - original_cor))
        }
        best_start <- which.min(subsample_scores)
        rvine_simulation <- rvine_simulation[best_start:(best_start + n - 1), , drop = FALSE]
      }
      
      # Improved back-transformation using empirical quantiles
      simulated_data <- matrix(NA, nrow = n, ncol = n_series)
      
      for (i in 1:n_series) {
        # Use empirical quantile function with interpolation
        simulated_data[, i] <- quantile(data[, i], 
                                        probs = pmax(pmin(rvine_simulation[, i], 
                                                          1 - 1e-10), 1e-10), 
                                        type = 8)
      }
      
      # Enhanced quality scoring using provided weights
      sim_cor <- cor(simulated_data, method = "kendall")
      sim_cor_pearson <- cor(simulated_data)
      
      # Correlation preservation (Kendall's tau)
      cor_error_tau <- mean(abs(sim_cor - original_cor))
      cor_error_pearson <- mean(abs(sim_cor_pearson - original_cor_pearson))
      
      # Distribution similarity using multiple metrics
      ks_stats <- numeric(n_series)
      for (i in 1:n_series) {
        # Suppress KS test warnings about ties (expected with continuous data)
        ks_result <- suppressWarnings(ks.test(data[, i], simulated_data[, i]))
        ks_stats[i] <- ks_result$statistic
      }
      mean_ks <- mean(ks_stats)
      
      # Moment preservation
      sim_means <- colMeans(simulated_data)
      sim_sds <- apply(simulated_data, 2, sd)
      mean_error <- mean(abs(sim_means - original_means) / abs(original_means))
      sd_error <- mean(abs(sim_sds - original_sds) / original_sds)
      
      # Combined score with user-defined weights
      score <- score_weights[1] * cor_error_tau + 
        score_weights[2] * cor_error_pearson + 
        score_weights[3] * mean_ks + 
        score_weights[4] * mean_error + 
        score_weights[5] * sd_error
      
      trial_scores[trial] <- score
      if (score < best_score) {
        best_score <- score
        best_sim <- simulated_data
      }
    }, error = function(e) {
      if (verbose) message(sprintf("Trial %d failed: %s", trial, e$message))
      trial_scores[trial] <- Inf
    })
  }
  
  if (is.null(best_sim)) {
    stop("All simulation trials failed. Check data quality and model specification.")
  }
  
  colnames(best_sim) <- series_names
  
  # Calculate comprehensive diagnostics
  sim_cor <- cor(best_sim, method = "kendall")
  sim_cor_pearson <- cor(best_sim)
  cor_error <- sim_cor - original_cor
  cor_error_pearson <- sim_cor_pearson - original_cor_pearson
  
  # Additional diagnostic statistics
  diagnostics <- list(
    # Correlation analysis
    original_correlation_tau = original_cor,
    simulated_correlation_tau = sim_cor,
    correlation_error_tau = cor_error,
    original_correlation_pearson = original_cor_pearson,
    simulated_correlation_pearson = sim_cor_pearson,
    correlation_error_pearson = cor_error_pearson,
    
    # Error metrics
    mean_absolute_error_tau = mean(abs(cor_error)),
    max_absolute_error_tau = max(abs(cor_error)),
    mean_absolute_error_pearson = mean(abs(cor_error_pearson)),
    max_absolute_error_pearson = max(abs(cor_error_pearson)),
    
    # Quality assessment
    quality_score = best_score,
    score_weights_used = score_weights,
    trial_scores = trial_scores,
    successful_trials = sum(is.finite(trial_scores)),
    
    # Model information
    RVM_model = RVM_U,
    n_observations = n_obs,
    n_variables = n_series,
    n_simulations = n,
    
    # Distribution comparison
    original_means = original_means,
    simulated_means = colMeans(best_sim),
    original_sds = original_sds,
    simulated_sds = apply(best_sim, 2, sd),
    
    # KS test results for distribution similarity
    ks_test_statistics = sapply(1:n_series, function(i) {
      suppressWarnings(ks.test(data[, i], best_sim[, i])$statistic)
    }),
    ks_test_pvalues = sapply(1:n_series, function(i) {
      suppressWarnings(ks.test(data[, i], best_sim[, i])$p.value)
    })
  )
  
  if (verbose) {
    message(sprintf("Best simulation achieved quality score: %.4f", best_score))
    message(sprintf("Score weights used: [%.1f, %.1f, %.1f, %.1f, %.1f]", 
                    score_weights[1], score_weights[2], score_weights[3], 
                    score_weights[4], score_weights[5]))
    message(sprintf("Mean absolute correlation error (Kendall): %.4f", 
                    diagnostics$mean_absolute_error_tau))
    message(sprintf("Mean absolute correlation error (Pearson): %.4f", 
                    diagnostics$mean_absolute_error_pearson))
  }
  
  # Create result object with class
  result <- list(
    simulated_data = best_sim,
    diagnostics = diagnostics,
    original_data = data,
    call = match.call()
  )
  
  class(result) <- "rvine_simulation"
  return(result)
}

#' Plot method for rvine_simulation objects
#'
#' Creates diagnostic plots for R-vine copula simulation results.
#'
#' @param x An object of class 'rvine_simulation'
#' @param type Type of plot: "distribution" (default), "correlation", or "both"
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot2 object or grid of plots
#'
#' @method plot rvine_simulation
#' @export
plot.rvine_simulation <- function(x, type = "distribution", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("reshape2 package required for plotting")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for plotting")
  }
  
  # Prepare data for plotting
  original_df <- data.frame(x$original_data, type = "Original")
  simulated_df <- data.frame(x$simulated_data, type = "Simulated")
  combined_df <- rbind(original_df, simulated_df)
  melted_df <- reshape2::melt(combined_df, id.vars = "type")
  
  if (type == "distribution") {
    # Distribution comparison plot
    p <- ggplot2::ggplot(melted_df, ggplot2::aes(x = value, fill = type)) +
      ggplot2::geom_density(alpha = 0.6) +
      ggplot2::facet_wrap(~variable, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Distribution Comparison: Original vs Simulated",
        x = "Value", 
        y = "Density",
        fill = "Data Type"
      ) +
      ggplot2::scale_fill_manual(values = c("Original" = "blue", "Simulated" = "red"))
    
    return(p)
    
  } else if (type == "correlation") {
    # Correlation comparison plot
    orig_cor <- x$diagnostics$original_correlation_pearson
    sim_cor <- x$diagnostics$simulated_correlation_pearson
    
    # Melt correlation matrices
    orig_melt <- reshape2::melt(orig_cor, varnames = c("Var1", "Var2"))
    sim_melt <- reshape2::melt(sim_cor, varnames = c("Var1", "Var2"))
    orig_melt$type <- "Original"
    sim_melt$type <- "Simulated"
    cor_df <- rbind(orig_melt, sim_melt)
    
    p1 <- ggplot2::ggplot(orig_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                    midpoint = 0, limit = c(-1, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Original Correlation Matrix", fill = "Correlation") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot2::ggplot(sim_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                    midpoint = 0, limit = c(-1, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Simulated Correlation Matrix", fill = "Correlation") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    return(gridExtra::grid.arrange(p1, p2, ncol = 2))
    
  } else if (type == "both") {
    # Both distribution and correlation plots
    p_dist <- plot.rvine_simulation(x, type = "distribution")
    p_corr <- plot.rvine_simulation(x, type = "correlation")
    
    # Create a simple correlation error plot for the grid
    cor_error <- x$diagnostics$correlation_error_pearson
    error_melt <- reshape2::melt(cor_error, varnames = c("Var1", "Var2"))
    p_error <- ggplot2::ggplot(error_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                    midpoint = 0, limit = c(-1, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Correlation Error (Original - Simulated)", fill = "Error") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    return(gridExtra::grid.arrange(p_dist, p_corr, p_error, ncol = 1, heights = c(2, 1, 1)))
  }
}

#' Print method for rvine_simulation objects
#'
#' @param x An object of class 'rvine_simulation'
#' @param ... Additional arguments passed to print
#'
#' @method print rvine_simulation
#' @export
print.rvine_simulation <- function(x, ...) {
  cat("R-vine Copula Simulation Results\n")
  cat("================================\n\n")
  cat(sprintf("Original observations: %d\n", x$diagnostics$n_observations))
  cat(sprintf("Variables: %d\n", x$diagnostics$n_variables))
  cat(sprintf("Simulated observations: %d\n", x$diagnostics$n_simulations))
  cat(sprintf("Quality score: %.4f\n", x$diagnostics$quality_score))
  cat(sprintf("Successful trials: %d/%d\n", 
              x$diagnostics$successful_trials, 
              length(x$diagnostics$trial_scores)))
  cat(sprintf("Mean absolute correlation error (Kendall): %.4f\n", 
              x$diagnostics$mean_absolute_error_tau))
  cat(sprintf("Mean absolute correlation error (Pearson): %.4f\n", 
              x$diagnostics$mean_absolute_error_pearson))
  cat("\nUse plot() to visualize results and $diagnostics for detailed metrics.\n")
}
