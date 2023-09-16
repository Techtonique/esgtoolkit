
# Smith-Wilson curve -----
tZC_SW <- function(yM = NULL, p = NULL, u, t, UFR, 
                   typeres=c("rates", "prices"), T_UFR = NULL)
{
  N <- length(u)  
  J <- length(t)
  
  fonctionsWilson <- function(t,u,alpha,UFR)
  {    
    N <- length(u)  
    J <- length(t)    
    u_mat <- matrix(rep.int(u,J), nrow=J, byrow=T)
    t_mat <- t(matrix(rep.int(t,N), nrow=N, byrow=T))
    min_u <- u_mat*(u_mat<=t_mat) + t_mat*(u_mat>t_mat)
    max_u <- u_mat + t_mat - min_u    
    return(exp(-UFR*(u_mat+t_mat))*(alpha*min_u - 0.5*exp(-alpha*max_u)*(exp(alpha*min_u)-exp(-alpha*min_u))))
  }
  
  if ((is.null(p) || missing(p)) && (!is.null(yM) || !missing(yM)))
  {
    p <- pricefromeuribor(t=0, T=u, L=yM)
  }
  
  alpha <- 0.1
  mu <- exp(-UFR*u)
  W <- fonctionsWilson(u,u,alpha,UFR)    
  Xhi <- solve(W)%*%(p-mu)
  W_interp <- fonctionsWilson(t,u,alpha,UFR)  
  P <- exp(-UFR*t) + W_interp %*%Xhi
  Fwd <- fwdrate(t[-J], t[-1], P[-J], P[-1])
  
  typeres <- match.arg(typeres)  
  if(max(t) > max(u))
  {
    if (missing(T_UFR)) stop("For the extrapolation T_UFR must be provided, check the output maturities")        
    alpha <- 0.1
    
    if (max(u)+T_UFR > max(t)) stop("Not enough extrapolation dates. Try a lower T_UFR")
    
    while(abs(Fwd[pmatch(max(u)+T_UFR, t)]- UFR) >= 0.0003)
    {
      alpha <- alpha+0.001
      mu <- exp(-UFR*u)
      W <- fonctionsWilson(u,u,alpha,UFR)    
      Xhi <- solve(W)%*%(p-mu)
      W_interp <- fonctionsWilson(t,u,alpha,UFR)  
      P <- exp(-UFR*t) + W_interp %*%Xhi      
      Fwd <- fwdrate(t[-J], t[-1], P[-J], P[-1])
    }    
  }
  
  if (typeres == "prices")
  {return (list(coefficients = list(alpha=alpha, Xhi=as.vector(Xhi)), 
                values = as.vector(P),
                fwd = Fwd))}
  
  if (typeres == "rates")
  {
    return (list(coefficients = list(alpha=alpha, Xhi=as.vector(Xhi)), values = euriborfromprice(t=0, T=t, ZC=P), fwd = Fwd))
  }
}
tZC_SW <- compiler::cmpfun(tZC_SW)


# Polynomial Hermite cubic spline interpolation -----
hermitecubicspline <- function(yM = NULL, p = NULL, matsin, matsout, 
                               typeres=c("rates", "prices"))
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

# utils -----

########## Tools
# simply coumpounded euribor rate
euriborfromprice <- function(t, T, ZC)
{
  # Brigo P. 7
  tau <- T-t
  return(as.vector((1/ZC - 1)/tau))
}
euriborfromprice <- compiler::cmpfun(euriborfromprice)


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