
# simulation of gaussian shocks --------------------------------------------


# util. function to reproduce results from package CDVine
CDRVineSim <- function(N, family, par=NULL, 
                       par2 = rep(0, length(par)), 
                       type=c("CVine", "DVine", "RVine"), 
                       RVM=NULL, U=NULL) 
{
  
  type <- match.arg(type)
  
  if (type == "RVine")
  {
    
    stopifnot(!is.null(RVM))
    
  } else {
    
    stopifnot(!(is.null(par)))
    
    # solve length(family) = d * (d - 1) / 2 for d
    d <- (1 + sqrt(1 + 8 * length(family))) / 2
    
    order <- 1:d
    
    RVM <- switch(type, 
                  "CVine" = VineCopula::C2RVine(order, family, par, par2), 
                  "DVine" = VineCopula::D2RVine(order, family, par, par2))
    
  }
    
  return(VineCopula::RVineSim(N, RVM, U))
}


# Underlying gaussian shocks for risk factors' simulation.
simshocks <- function(n, horizon, 
                      frequency = c("annual", "semi-annual", 
                                    "quarterly", "monthly", 
                                    "weekly", "daily"), 
                      method = c("classic", "antithetic", 
                                 "mm", "hybridantimm", "TAG"), 
                      family = NULL, 
                      par = NULL, par2 = rep(0, length(par)), 
                      RVM = NULL, 
                      type = c("CVine", "DVine", "RVine"),
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
                  "monthly" = 1/12L,
                  "weekly" = 1/52L,
                  "daily" = 1/252L)
  type <- match.arg(type)
  
  method <- match.arg(method)
  m <- horizon/delta  
  type <- match.arg(type)
  
  set.seed(seed)
  
  if (is.null(family) && is.null(par)) # independent
  {
    if (method == "classic")
    {
      if (!is.null(start))
        return(ts(data = rnormESGcpp(N = n, M = m), 
           start = start, deltat = delta))
        
      return(ts(data = rnormESGcpp(N = n, M = m), 
                start = delta, deltat = delta))  
    }
    
    if (method == "antithetic")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")        
      temp <- rnormESGcpp(N = half.n, M = m)        
      
      if (!is.null(start))
        return(ts(data = cbind(temp, -temp), 
                  start = start, deltat = delta))
      
      return(ts(data = cbind(temp, -temp), 
                start = delta, deltat = delta))
    }
    
    if (method == "mm") 
    {
      if (!is.null(start))
        return(ts(data = scaleESG(rnormESGcpp(N = n, M = m)), 
                  start = start, deltat = delta))  
        
      return(ts(data = scaleESG(rnormESGcpp(N = n, M = m)), 
                start = delta, deltat = delta))  
    }
    
    if (method == "hybridantimm")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")
      temp <- rnormESGcpp(N = half.n, M = m)
      
      if (!is.null(start))
        return(ts(data = scaleESG(cbind(temp, -temp)), 
                  start = start, deltat = delta))
      
      return(ts(data = scaleESG(cbind(temp, -temp)), 
                start = delta, deltat = delta))
    }
    
    if (method == "TAG")
    {
      if (!is.null(start))
        return(ts(data = TAG(n, m), 
                  start = start, deltat = delta))        
        
      return(ts(data = TAG(n, m), 
                start = delta, deltat = delta))        
    }    
  } 
  else # !is.null(family) && !is.null(par) # C, D, R-Vine
  {
    nb.sim <- n*m

    if (method == "classic")
    {
      shocks.sim <- qnorm(CDRVineSim(N = nb.sim, RVM=RVM,
                                    family, par, par2, type))
      d <- dim(shocks.sim)[2]
      
      if (!is.null(start))
        return(lapply(1:d, function(i)
          ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n),
             start = start, deltat = delta)))
      
      return(lapply(1:d, function(i)
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n),
           start = delta, deltat = delta)))
    }
    
    if (method == "antithetic")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")
      temp <- qnorm(CDRVineSim(N = half.nb.sim, RVM=RVM,
                              family, par, par2, type))
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]
      
      if (!is.null(start))
        return(lapply(1:d, function(i)
          ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n),
             start = start, deltat = delta)))
      
      return(lapply(1:d, function(i)
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n),
           start = delta, deltat = delta)))
    }

    if (method == "mm")
    {
      shocks.sim <- qnorm(CDRVineSim(N = nb.sim, RVM=RVM,
                                    family, par, par2, type))
      d <- dim(shocks.sim)[2]
      
      if (!is.null(start))
        return(lapply(1:d, function(i)
          ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)),
             start = start, deltat = delta)))
        
      return(lapply(1:d, function(i)
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)),
           start = delta, deltat = delta)))
    }

    if (method == "hybridantimm")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")
      temp <- qnorm(CDRVineSim(N = half.nb.sim, RVM=RVM,
                              family, par, par2, type))
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]
      
      if (!is.null(start))
        return(lapply(1:d, function(i)
          ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)),
             start = start, deltat = delta)))
      
      return(lapply(1:d, function(i)
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)),
           start = delta, deltat = delta)))
    }

    if (method == "TAG")
    {
      if (!is.null(family) && family > 0) warning("for method == 'TAG', only the independence copula is implemented")
      
      if (!is.null(start))
        return(lapply(1:d, function(i)
          ts(data = TAG(n = n, m = m),
             start = start, deltat = delta)))
        
      return(lapply(1:d, function(i)
        ts(data = TAG(n = n, m = m),
           start = delta, deltat = delta)))
    }
  }
}