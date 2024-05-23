# Stochastic discount factors or discounted values -----
esgdiscountfactor <- function(r, X)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  length_r <- length(r)
  length_X <- length(X)
  start_r <- start(r)
  deltat.r <- deltat(r)
  start_X <- start(X)
  deltat_X <- deltat(X)
  
  if(length_r == 1 && length_X == 1) 
  {
    return(X*exp(-r)) 
  }
  
  if(length_r == 1 && length_X != 1) 
  {
    r <- ts(matrix(r, nrow(X), ncol(X)), 
            start = start_X, deltat = deltat_X)
    
    if(tsp(X)[1] > 0)
    {
      return(ts(X*exp(-apply(r, 2, cumsum)*deltat_X), 
                start = 0, 
                deltat = deltat_X))      
    }
    else
    {
      Int_r <- exp(-apply(r, 2, cumsum)*deltat_X)
      return(ts(X*rbind(rep(1, ncol(X)), 
                        Int_r[1:(nrow(X)-1), ]),                                            
                start = 0, 
                deltat = deltat_X))
    }
  } 
  
  if(length_r != 1 && length_X == 1)
  {
    X <- ts(matrix(X, nrow(r), ncol(r)), 
            start = start_r, deltat = deltat.r)
    
    return(ts(X*exp(-apply(r, 2, cumsum)*deltat.r),                                          
              start = 0, 
              deltat = deltat.r))
  }
  
  if(length_r != 1 && length_X != 1)
  {
    if (all.equal(start(X), start(r)) && all.equal(frequency(X), frequency(r)))
    {
      return(X*exp(-apply(r, 2, cumsum)*deltat(r)))
    }
    
    if(tsp(X)[1] > 0)
    {
      return(suppressWarnings(ts(X*window(exp(-apply(r, 2, cumsum)*deltat.r), 
                                          start = start_X, 
                                          deltat = deltat_X), 
                                 start = 0, 
                                 deltat = deltat_X)))
    }
    else
    {
      Int_r <- ts(exp(-apply(r, 2, cumsum)*deltat.r), deltat = deltat.r)
      return(suppressWarnings(ts(X*rbind(rep(1, ncol(X)), 
                                         Int_r[1:(nrow(X)-1), ]),
                                 start = 0, 
                                 deltat = deltat_X)))
    }
  }
}


# Estimation of discounted asset prices -----
esgmcprices <- function(r, X, maturity = NULL)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  maturity.out <- maturity
  
  if(is.ts(X) && tsp(X)[1] > 0 && !is.null(maturity))
  {
    maturity.out <- maturity - deltat(X)
  }
  
  Y <- esgdiscountfactor(r, X)
  
  if(length(r) == 1 && length(X) == 1) 
  {
    return(Y) 
  }
  
  Z <- ts(rowMeans(Y), 
          start = start(Y), 
          deltat = deltat(Y))
  
  if(!is.null(maturity))
  {
    return(window(Z, start = maturity.out, end = maturity.out))
  }
  else  
  {
    return(Z)
  }
  
}

