
# Smith-Wilson
tZC_SW <- function(yM = NULL, p = NULL, u, t, UFR, typeres=c("rates", "prices"), T_UFR = NULL)
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
