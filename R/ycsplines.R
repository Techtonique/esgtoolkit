

# Polynomial Hermite cubic spline interpolation
hermitecubicspline <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
  u <- matsin
  t <- matsout
  
  if (max(t) > max(u)) stop("Only interpolation can be performed with cubic splines, check the output maturities")
  
  if (!is.null(yM) && !is.null(p)) stop("either yields OR prices must be provided")
  
  if (is.null(yM) && is.null(p)) stop("either yields or prices must be provided")
  
  if (!is.null(yM))
  {Sw <- yM}
  else{Sw <- p}
  
  n <- length(u)
  i <- seq_len(n)
  iup <- i[-1]
  idown <- i[-n]  
  Delta <- diff(u)
  delta <- diff(Sw)
  
  A <- delta/Delta    
  A_up <- A[iup]
  A_down <- A[idown] 
  
  Delta_up <- Delta[iup]
  Delta_down <- Delta[idown] 
  
  P <- (A_up*Delta_down + A_down*Delta_up)/(Delta_down + Delta_up)
  P <- c(2*A[1] - P[2], P[!is.na(P)], 0)
  
  m <- length(t)
  z <- rep.int(0, m)
  SWout <- rep.int(0, m)
  
  for (i in  seq_len(m))
  {          
    u_down <-  pmatch(floor(t[i]),u)
    u_up <-  pmatch(ceiling(t[i]),u)
    Sw_down <- Sw[u_down]
    Sw_up <- Sw[u_up]
    P_down <- P[u_down]
    P_up <- P[u_up]
    Delta_down <- Delta[u_down]
    Delta_up <- Delta[u_up]
    
    if (Sw_down != Sw_up) 
    {
      z <- (t[i] - u_down)/(u_up - u_down)                    
      SWout[i] <- Sw_down*((1-z)^3) + (3*Sw_down + Delta_down*P_down)*
        ((1-z)^2)*z + (3*Sw_up - Delta_down*P_up)*(1-z)*(z^2) + 
        Sw_up*(z^3)      
    }
    else 
    {
      SWout[i] <- Sw_down
    }    
  }
  
  if (typeres=="rates")
  {P <- pricefromeuribor(0, t, SWout)}
  else {P <- SWout}
  
  Fwd <- fwdrate(t[-m], t[-1], P[-m], P[-1])
  
  if ((!is.null(yM) && typeres == "rates") || (!is.null(p) && typeres == "prices"))
  {return (list(coefficients = NA, values = SWout, fwd = Fwd))}
  
  if (!is.null(yM) && typeres == "prices")
  {
    return (list(coefficients = NA, values = pricefromeuribor(0, t, SWout), fwd = Fwd))
  }
  
  if (!is.null(p) && typeres == "rates")
  {
    return (list(coefficients = NA, values = euriborfromprice(0, t, SWout), fwd = Fwd))
  }
}
hermitecubicspline <- compiler::cmpfun(hermitecubicspline)
