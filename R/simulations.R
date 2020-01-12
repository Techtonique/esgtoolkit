
# simulation of gaussian shocks --------------------------------------------

#'@title 
#'
#'Underlying gaussian shocks for risk factors' simulation.
#'
#'@description 
#'
#'This function makes simulations of correlated or dependent gaussian shocks for risk factors.
#'
#'@details The function shall be used along with \code{\link{simdiff}}, in order to embed  
#'correlated or dependent random gaussian shocks into simulated diffusions. 
#'\code{\link{esgplotshocks}} can help in visualizing the type of dependence 
#'between the shocks. 
#'
#'@param n number of independent observations for each risk factor.
#'
#'@param horizon horizon of projection.
#'
#'@param frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252).
#'
#'@param method either classic monte carlo, antithetic variates, moment matching, 
#'hybrid antithetic variates + moment matching or "TAG" (see the 4th reference for 
#'the latter). 
#'
#'@param family the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 integer vector of C-/D-vine pair-copula families with values 
#' 0 = independence copula, 
#' 1 = Gaussian copula, 
#' 2 = Student t copula (t-copula), 
#' 3 = Clayton copula, 
#' 4 = Gumbel copula, 
#' 5 = Frank copula, 
#' 6 = Joe copula, 
#' 7 = BB1 copula, 
#' 8 = BB6 copula, 
#' 9 = BB7 copula, 
#' 10 = BB8 copula, 
#' 13 = rotated Clayton copula (180 degrees; "survival Clayton"), 
#' 14 = rotated Gumbel copula (180 degrees; "survival Gumbel"), 
#' 16 = rotated Joe copula (180 degrees; "survival Joe"), 
#' 17 = rotated BB1 copula (180 degrees; "survival BB1"),
#' 18 = rotated BB6 copula (180 degrees; "survival BB6"),
#' 19 = rotated BB7 copula (180 degrees; "survival BB7"),
#' 20 = rotated BB8 copula (180 degrees; "survival BB8"),
#' 23 = rotated Clayton copula (90 degrees), 
#' 24 = rotated Gumbel copula (90 degrees),
#' 26 = rotated Joe copula (90 degrees), 
#' 27 = rotated BB1 copula (90 degrees), 
#' 28 = rotated BB6 copula (90 degrees), 
#' 29 = rotated BB7 copula (90 degrees), 
#' 30 = rotated BB8 copula (90 degrees), 
#' 33 = rotated Clayton copula (270 degrees), 
#' 34 = rotated Gumbel copula (270 degrees), 
#' 36 = rotated Joe copula (270 degrees), 
#' 37 = rotated BB1 copula (270 degrees), 
#' 38 = rotated BB6 copula (270 degrees), 
#' 39 = rotated BB7 copula (270 degrees), 
#' 40 = rotated BB8 copula (270 degrees)  
#'
#'@param par the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 vector of pair-copula parameters.
#'
#'@param par2 the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 vector of second parameters for pair-copula families with two 
#'parameters (t, BB1, BB6, BB7, BB8; no default).
#'
#'@param type type of the vine model:
#' 1 : C-vine
#' 2 : D-vine
#'
#'@param seed reproducibility seed
#'
#'@return
#'If \code{family} and \code{par} are not provided, a univariate time 
#'series object with simulated gaussian shocks for one risk factor. Otherwise, 
#'a list of time series objects, containing gaussian shocks for each risk factor. 
#'
#'@author Thierry Moudiki
#'
#'@references
#'
#'Brechmann, E., Schepsmeier, U. (2013). Modeling Dependence with C-
#'and D-Vine Copulas: The R Package CDVine. Journal of Statistical Software,
#'52(3), 1-27. URL \url{http://www.jstatsoft.org/v52/i03/}.
#'
#'Genz, A. Bretz, F., Miwa, T. Mi, X., Leisch, F., Scheipl, F., Hothorn, T. (2013).
#' mvtnorm: Multivariate Normal and t Distributions. R package version 0.9-9996.
#'
#'Genz, A. Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. 
#'Lecture Notes in Statistics, Vol. 195., Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.
#'
#'Nteukam T, O., & Planchet, F. (2012). Stochastic evaluation of life insurance 
#'contracts: Model point on asset trajectories and measurement of the error 
#'related to aggregation. Insurance: Mathematics and Economics, 51(3), 624-631.
#'URL \url{http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/ab539dcebcc4e77ac12576c6004afa67/$FILE/Article_US_v1.5.pdf}
#'
#'@export
#'
#'@seealso \code{\link{simdiff}}, \code{\link{esgplotshocks}}
#'
#'@examples
#'
#'# Number of risk factors
#'d <- 6
#'
#'# Number of possible combinations of the risk factors
#'dd <- d*(d-1)/2
#'
#'# Family : Gaussian copula for all
#'fam1 <- rep(1,dd)
#'
#'# Correlation coefficients between the risk factors (d*(d-1)/2)
#'par1 <- c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,
#'         -0.2,-0.32,-0.19,-0.17,-0.06)
#'
#'                  
#'# Simulation of shocks for the 6 risk factors
#'simshocks(n = 10, horizon = 5, family = fam1, par = par1)
#'
#'
#'# Simulation of shocks for the 6 risk factors
#'# on a quarterly basis
#'simshocks(n = 10, frequency = "quarterly", horizon = 2, family = fam1, 
#'par = par1)
#'
#'
#'# Simulation of shocks for the 6 risk factors simulation
#'# on a quarterly basis, with antithetic variates and moment matching. 
#'s0 <- simshocks(n = 10, method = "hyb", horizon = 4, 
#'family = fam1, par = par1)
#'
#'  
#' s0[[2]]
#' colMeans(s0[[1]])
#' colMeans(s0[[5]])
#' apply(s0[[3]], 2, sd)
#' apply(s0[[4]], 2, sd)
#'
simshocks <- function(n, horizon, 
                      frequency = c("annual", "semi-annual", 
                                    "quarterly", "monthly", 
                                    "weekly", "daily"), 
                      method = c("classic", "antithetic", 
                                 "mm", "hybridantimm", "TAG"), 
                      family = NULL, par = NULL, par2 = NULL, type = c("CVine", "DVine"),
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
  method <- match.arg(method)
  m <- horizon/delta  
  type <- match.arg(type)
  
  set.seed(seed)
  
  if (is.null(family) && is.null(par))
  {
    if (method == "classic")
    {
      return(ts(data = rnormESGcpp(N = n, M = m), 
                start = delta, end = horizon, deltat = delta))  
    }
    
    if (method == "antithetic")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")        
      temp <- rnormESGcpp(N = half.n, M = m)        
      return(ts(data = cbind(temp, -temp), 
                start = delta, end = horizon, deltat = delta))
    }
    
    if (method == "mm") 
    {
      return(ts(data = scaleESG(rnormESGcpp(N = n, M = m)), 
                start = delta, end = horizon, deltat = delta))  
    }
    
    if (method == "hybridantimm")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")
      temp <- rnormESGcpp(N = half.n, M = m)
      return(ts(data = scaleESG(cbind(temp, -temp)), 
                start = delta, end = horizon, deltat = delta))
    }
    
    if (method == "TAG")
    {
      return(ts(data = TAG(n, m), 
                start = delta, end = horizon, deltat = delta))        
    }    
  } 
  else # !is.null(family) && !is.null(par)
  { 
    nb.sim <- n*m
    
    if (method == "classic")
    {
      shocks.sim <- qnorm(CDVineSim(N = nb.sim, 
                                    family, par, par2, type)) 
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n), 
           start = delta, end = horizon, deltat = delta)))        
    }
    
    if (method == "antithetic")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")      
      temp <- qnorm(CDVineSim(N = half.nb.sim, 
                              family, par, par2, type))         
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n), 
           start = delta, end = horizon, deltat = delta)))
    }
    
    if (method == "mm") 
    {
      shocks.sim <- qnorm(CDVineSim(N = nb.sim, 
                                    family, par, par2, type)) 
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)), 
           start = delta, end = horizon, deltat = delta))) 
    }
    
    if (method == "hybridantimm")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")
      temp <- qnorm(CDVineSim(N = half.nb.sim, 
                              family, par, par2, type)) 
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)), 
           start = delta, end = horizon, deltat = delta)))  
    }
    
    if (method == "TAG")
    {
      if (!is.null(family) && family > 0) warning("for method == 'TAG', only the independence copula is implemented")   
      return(lapply(1:d, function(i) 
        ts(data = TAG(n = n, m = m), 
           start = delta, end = horizon, deltat = delta)))            
    }    
  }
}
simshocks <- cmpfun(simshocks)


# simulation of risk factors ----------------------------------------------

#'@title 
#'
#'Simulation of diffusion processes. 
#'
#'@description 
#'
#'This function makes simulations of diffusion processes, that are building 
#'blocks for various risk factors' models. 
#'
#'@param n number of independent observations.
#'
#'@param horizon horizon of projection.
#'
#'@param frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252).
#'
#'@param model either Geometric Brownian motion-like (\code{"GBM"}), 
#'Cox-Ingersoll-Ross (\code{"CIR"}), or Ornstein-Uhlenbeck 
#'(\code{"OU"}).
#'
#' GBM-like (GBM, Merton, Kou, Heston, Bates)
#'
#'\deqn{dX_t = \theta_1(t) X_t dt + \theta_2(t) X_t dW_t + X_t JdN_t}
#'  
#'  CIR
#'
#'\deqn{dX_t = (\theta_1 - \theta_2 X_t) dt + \theta_3\sqrt(X_t) dW_t}
#'  
#'  Ornstein-Uhlenbeck
#'\deqn{dX_t = (\theta_1 - \theta_2 X_t)dt + \theta_3 dW_t}
#'  
#'Where \eqn{(W_t)_t} is a standard brownian motion :
#'  
#'  \deqn{dW_t ~~ \epsilon \sqrt(dt)}
#'  
#'and \deqn{\epsilon ~~ N(0, 1)}
#'
#'The \eqn{\epsilon} is a gaussian increment that 
#'can be an output from \code{\link{simshocks}}. 
#'
#'For 'GBM-like', \eqn{\theta_1} and \eqn{\theta_2} can be held constant, and the
#' jumps part \eqn{JdN_t} is optional. In case the jumps are used, they arise 
#' following a Poisson process \eqn{(N_t)}, with intensities \eqn{J} drawn either 
#'from lognormal or asymmetric double-exponential distribution. 
#'
#'@param x0 starting value of the process. 
#'
#'@param theta1 a \code{numeric} for \code{model = "GBM"}, \code{model = "CIR"},
#'\code{model = "OU"}. Can also be a time series object (an output from 
#'\code{simdiff} with the same number of scenarios, horizon and frequency) for 
#'\code{model = "GBM"}, and time-varying parameters. 
#'
#'@param theta2 a \code{numeric} for \code{model = "GBM"}, \code{model = "CIR"},
#'\code{model = "OU"}. Can also be a time series object (an output from 
#'\code{simdiff} with the same number of scenarios, horizon and frequency) 
#'for \code{model = "GBM"}, and time-varying parameters.
#'
#'@param theta3 a \code{numeric}, volatility for \code{model = "CIR"} and 
#'\code{model = "OU"}.
#'
#'@param lambda intensity of the Poisson process counting the jumps. Optional.
#'
#'@param mu_z mean parameter for the lognormal jumps size. Optional.
#'
#'@param sigma_z standard deviation parameter for the lognormal jumps size. 
#'Optional.
#'
#'@param p probability of positive jumps. Must belong to ]0, 1[. Optional.
#'
#'@param eta_up mean of positive jumps in Kou's model. Must belong to 
#']0, 1[. Optional.
#'
#'@param eta_down mean of negative jumps. Must belong to ]0, 1[. Optional.
#'
#'@param eps gaussian shocks. If not provided, independent shocks are 
#'generated internally by the function. Otherwise, for custom shocks, 
#'must be an output from \code{\link{simshocks}}. 
#'
#'@param seed reproducibility seed
#'
#'@return a time series object. 
#'
#'@seealso \code{\link{simshocks}}, \code{\link{esgplotts}}
#'
#'@author Thierry Moudiki
#'
#'@references
#'
#'Black, F., Scholes, M.S. (1973) The pricing of options and corporate liabilities, 
#'Journal of Political Economy, 81, 637-654.
#'
#'Cox, J.C., Ingersoll, J.E., Ross, S.A. (1985) A theory of the term structure of
#' interest rates, Econometrica, 53, 385-408.
#'
#'Iacus, S. M. (2009). Simulation and inference for stochastic differential 
#'equations: with R examples (Vol. 1). Springer.
#'
#'Glasserman, P. (2004). Monte Carlo methods in financial engineering (Vol. 53). 
#'Springer.
#'
#'Kou S, (2002), A jump diffusion model for option pricing, Management Sci-
#'ence Vol. 48, 1086-1101.
#'
#'Merton, R. C. (1976). Option pricing when underlying stock returns are 
#'discontinuous. Journal of financial economics, 3(1), 125-144.
#'
#'Uhlenbeck, G. E., Ornstein, L. S. (1930) On the theory of Brownian motion, 
#'Phys. Rev., 36, 823-841.
#'
#'Vasicek, O. (1977) An Equilibrium Characterization of the Term Structure, 
#'Journal of Financial Economics, 5, 177-188.
#'
#'@export
#'
#'@examples
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma_v <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma_v
#'
#'# OU
#'
#'sim.OU <- simdiff(n = 10, horizon = 5, 
#'                frequency = "quart",  
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'head(sim.OU)
#'par(mfrow=c(2,1))
#'esgplotbands(sim.OU, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.OU), sim.OU, type = 'l', main = "with matplot")                
#'
#'
#'# OU with simulated shocks (check the dimensions)
#'
#'eps0 <- simshocks(n = 50, horizon = 5, frequency = "quart", method = "anti")
#'sim.OU <- simdiff(n = 50, horizon = 5, frequency = "quart",   
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3, 
#'                eps = eps0)
#'par(mfrow=c(2,1))
#'esgplotbands(sim.OU, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.OU), sim.OU, type = 'l', main = "with matplot")
#'# a different plot
#'esgplotts(sim.OU)
#'
#'
#'# CIR
#'
#'sim.CIR <- simdiff(n = 50, horizon = 5, 
#'                frequency = "quart",  
#'                model = "CIR", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = 0.05)
#'esgplotbands(sim.CIR, xlab = "time", ylab = "values", main = "with esgplotbands")                  
#'matplot(time(sim.CIR), sim.CIR, type = 'l', main = "with matplot")
#'
#'
#'
#'# GBM
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'esgplotbands(sim.GBM, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.GBM), sim.GBM, type = 'l', main = "with matplot")
#'
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)                
#'esgplotbands(sim.GBM, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.GBM), sim.GBM, type = 'l', main = "with matplot")
simdiff <- function(n, horizon, 
                    frequency = c("annual", "semi-annual", 
                                  "quarterly", "monthly", 
                                  "weekly", "daily"), 
                    model = c("GBM", "CIR", "OU"), 
                    x0, theta1 = NULL, theta2 = NULL, theta3 = NULL,
                    lambda = NULL, mu_z = NULL, sigma_z = NULL, 
                    p = NULL, eta_up = NULL, eta_down = NULL,
                    eps = NULL, seed = 123)
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
    
    # rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)
    return(ts(rCIRESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                              theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, end = horizon, deltat = delta))
  }
  
  if (model == "OU")
  {
    if (is.null(theta1) || is.null(theta2) || is.null(theta3))
      stop("'theta1', 'theta2' and 'theta3' must be provided")
    
    if(length(theta1) != 1 || length(theta2) != 1 || length(theta3) != 1) 
      stop(paste("length required for 'theta1', 'theta2' and 'theta3' : 1"))
    
    # rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)     
    return(ts(rOUESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                             theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, end = horizon, deltat = delta))
  }
  
  if (model == "GBM")
  {
    if (is.null(mu_z) && is.null(sigma_z) && is.null(p) && is.null(eta_up) && is.null(eta_down))
    {
      # rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, 
      # NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps)       
      return(ts(rGBMESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                           theta1 = theta1, theta2 = theta2, eps = eps), 
                start = 0, end = horizon, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(mu_z) && !is.null(sigma_z))
    {
      # rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double mu, const double sigma, 
      # NumericMatrix eps)       
      return(ts(rGBMjumpsnormESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                    theta1 = theta1, theta2 = theta2,
                                    lambda = lambda, mu = mu_z, sigma = sigma_z, 
                                    eps = eps), start = 0, end = horizon, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(p) && !is.null(eta_up) && !is.null(eta_down))
    {
      # rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double eta_up, const double eta_down,
      # const double p, NumericMatrix eps)
      return(ts(rGBMjumpskouESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                   theta1 = theta1, theta2 = theta2,
                                   lambda = lambda, eta_up = eta_up, eta_down = eta_down,
                                   p = p, eps = eps), 
                start = 0, end = horizon, deltat = delta))
    }
  }
  }                                                                                                                                  
simdiff <- cmpfun(simdiff)


