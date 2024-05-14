# simulation of risk factors ----------------------------------------------


# Simulation of diffusion processes. 
simdiff <- function(n, horizon, 
                    frequency = c("annual", "semi-annual", 
                                  "quarterly", "monthly", 
                                  "weekly", "daily"), 
                    model = c("GBM", "CIR", "OU"), 
                    x0, theta1 = NULL, theta2 = NULL, theta3 = NULL,
                    lambda = NULL, mu_z = NULL, sigma_z = NULL, 
                    p = NULL, eta_up = NULL, eta_down = NULL,
                    eps = NULL, 
                    start = NULL,
                    seed = 123)
{
  if (floor(n) != n) stop("'n' must be an integer")
  if (n <= 1) stop("'n' must be > 1")
  if (floor(horizon) != horizon) stop("'horizon' must be an integer")  
  frequency <- match.arg(frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  m <- horizon/delta 
  nbdates <- m + 1
  model <- match.arg(model)
  
  set.seed(seed)
  
  if (is.null(eps))
  { 
    eps <- rnormESGcpp(N = n, M = m)
  } 
  else
  {
    if(!all.equal(dim(eps), c(m, n))) 
      stop(paste("dimensions required for 'eps' : nrow (observation dates) = "
                 , m, "\n and ncol (number of scenarios) = ", n, ".\n Use 'simshocks' function
                 with same number of observations, horizon and frequency to obtain 'eps'"))
  }
  
  if(model == "GBM")
  {
    if (!is.null(theta3)) 
      warning("unused parameter : 'theta3'")
    
    if(length(theta1) == 1) 
    {
      theta1 <- matrix(data = theta1, nrow = nbdates, ncol = n)
    }
    else
    {
      if(prod(dim(theta1) == c(nbdates, n)) == 0) 
        stop(paste("dimensions required for 'theta1' : a numeric or nrow (observation dates) = "
                   , nbdates, "\n and ncol (number of scenarios) = ", n))
    }
    
    if(length(theta2) == 1) 
    {
      theta2 <- matrix(data = theta2, nrow = nbdates, ncol = n)
    }
    else
    {
      if(prod(dim(theta2) == c(nbdates, n)) == 0) 
        stop(paste("dimensions required for 'theta2' : a numeric or nrow (observation dates) = "
                   , nbdates, "\n and ncol (number of scenarios) = ", n))
    }   
  }
  
  if (model == "CIR")
  {
    if (is.null(theta1) || is.null(theta2) || is.null(theta3))
      stop("'theta1', 'theta2' and 'theta3' must be provided")
    
    if(length(theta1) != 1 || length(theta2) != 1 || length(theta3) != 1) 
      stop(paste("length required for 'theta1', 'theta2' and 'theta3' : 1"))
    
    if (!is.null(start))
      {# rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, 
      # NumericVector theta, NumericMatrix eps)
      return(ts(rCIRESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                theta = c(theta1, theta2, theta3), eps = eps), 
                start = start, deltat = delta))}
    
    # rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)
    return(ts(rCIRESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                              theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, deltat = delta))
  }
  
  if (model == "OU")
  {
    if (is.null(theta1) || is.null(theta2) || is.null(theta3))
      stop("'theta1', 'theta2' and 'theta3' must be provided")
    
    if(length(theta1) != 1 || length(theta2) != 1 || length(theta3) != 1) 
      stop(paste("length required for 'theta1', 'theta2' and 'theta3' : 1"))
    
    if (!is.null(start))
    {
      # rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, 
      # NumericVector theta, NumericMatrix eps)     
      return(ts(rOUESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                               theta = c(theta1, theta2, theta3), eps = eps), 
                start = start, deltat = delta))
    }
    
    # rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)     
    return(ts(rOUESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                             theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, deltat = delta))
  }
  
  if (model == "GBM")
  {
    if (is.null(mu_z) && is.null(sigma_z) && is.null(p) && is.null(eta_up) && is.null(eta_down))
    {
      if (!is.null(start))
      {
        # rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, 
        # NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps)       
        return(ts(rGBMESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                             theta1 = theta1, theta2 = theta2, eps = eps), 
                  start = start, deltat = delta))
      }
        
      # rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, 
      # NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps)       
      return(ts(rGBMESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                           theta1 = theta1, theta2 = theta2, eps = eps), 
                start = 0, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(mu_z) && !is.null(sigma_z))
    {
      if (!is.null(start))
      {
        # rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta,
        # const double x0, NumericMatrix theta1, NumericMatrix theta2,
        # const double lambda, const double mu, const double sigma, 
        # NumericMatrix eps)       
        return(ts(rGBMjumpsnormESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                      theta1 = theta1, theta2 = theta2,
                                      lambda = lambda, mu = mu_z, sigma = sigma_z, 
                                      eps = eps), start = start, deltat = delta))
      }
        
      # rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double mu, const double sigma, 
      # NumericMatrix eps)       
      return(ts(rGBMjumpsnormESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                    theta1 = theta1, theta2 = theta2,
                                    lambda = lambda, mu = mu_z, sigma = sigma_z, 
                                    eps = eps), start = 0, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(p) && !is.null(eta_up) && !is.null(eta_down))
    {
      if (!is.null(start))
      {
        # rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta,
        # const double x0, NumericMatrix theta1, NumericMatrix theta2,
        # const double lambda, const double eta_up, const double eta_down,
        # const double p, NumericMatrix eps)
        return(ts(rGBMjumpskouESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                     theta1 = theta1, theta2 = theta2,
                                     lambda = lambda, eta_up = eta_up, eta_down = eta_down,
                                     p = p, eps = eps), 
                  start = start, deltat = delta))
      }
        
      # rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double eta_up, const double eta_down,
      # const double p, NumericMatrix eps)
      return(ts(rGBMjumpskouESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                   theta1 = theta1, theta2 = theta2,
                                   lambda = lambda, eta_up = eta_up, eta_down = eta_down,
                                   p = p, eps = eps), 
                start = 0, deltat = delta))
    }
  }
}                                                                                                                                  

