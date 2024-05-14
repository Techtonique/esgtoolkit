########## Tools
# simply coumpounded euribor rate
euriborfromprice <- function(t, T, ZC)
{
  # Brigo P. 7
  tau <- T-t
  return(as.vector((1/ZC - 1)/tau))
}
euriborfromprice <- compiler::cmpfun(euriborfromprice)

#' @export
debug_print <- function(x) {
  cat("\n")
  print(paste0(deparse(substitute(x)), "'s value:"))
  print(x)
  cat("\n")
}

# simply coumpounded zero-coupon price
pricefromeuribor <- function(t, T, L)
{
  # Brigo P. 7
  tau <- T-t
  return(as.vector(1/(1 + L*tau)))
}
pricefromeuribor <- compiler::cmpfun(pricefromeuribor)

# simply coumpounded forward rate
fwdrate <- function(T_, S, ZC_T, ZC_S)
{
  # Brigo P. 12
  tau <- S-T_
  return((ZC_T/ZC_S - 1)/tau)
}
fwdrate <- compiler::cmpfun(fwdrate)

# Swap curve bootstrap
bootstrapswapcurve <- function(swaprate, maturity, typeres=c("rates", "prices"))
{  
  swaprate <- swaprate/100
  accrual <- diff(maturity)
  n <- length(swaprate)  
  ZC <- numeric(n)  
  ZC[1] <- 1/(1+maturity[1]*swaprate[1])
  
  for (i in seq_len(n)[-1])
  {
    ind <- seq_len(i-1)
    ZC[i] <- (1 - swaprate[i]*sum(accrual[ind]*ZC[ind]))/(1+swaprate[i]*accrual[i-1])
  }
  
  typeres <- match.arg(typeres)
  
  if (missing(typeres) || typeres == "prices")
  {return(ts(ZC))}
  
  if (typeres == "rates")
  {return((1-ZC)/(maturity*ZC))}  
}
bootstrapswapcurve <- compiler::cmpfun(bootstrapswapcurve)