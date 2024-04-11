

# Svensson interpolation
tZC_SV <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
  #   require(randtoolbox)
  #   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=0, T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {   	
    gam1 <- matsin/betaV[5]
    gam2 <- matsin/betaV[6]
    aux1 <- 1-exp(-gam1)
    aux2 <- 1-exp(-gam2)      
    y <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
      betaV[4]*(aux2/gam2+aux2-1)
    
    return(crossprod(y - yM))
  }
  OF <- compiler::cmpfun(OF)
  set.seed(123)    
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(0, -15, -30, -30 ,0 ,3),
                           upper = c (15, 30, 30, 30 ,3 ,6),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[5]
  gam2 <- matsout/betaV[6]
  aux1 <- 1-exp(-gam1)
  aux2 <- 1-exp(-gam2)  
  L <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
    betaV[4]*(aux2/gam2+aux2-1)
  aux3 <- 1 - aux1
  aux4 <- 1 - aux2
  Fwd <- betaV[1] + betaV[2]*aux3 + betaV[3]*gam1*aux3 + betaV[4]*gam2*aux4
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=rep.int(0, J), T=matsout, L=L), fwd = Fwd))
  }
}
tZC_SV <- compiler::cmpfun(tZC_SV)


# Svensson extrapolation
tZC_SV_extra <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"), UFR)
{  
  #   require(randtoolbox)
  #   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=0, T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[5]
    gam2 <- matsin/betaV[6]
    aux1 <- 1-exp(-gam1)
    aux2 <- 1-exp(-gam2)      
    y <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
      betaV[4]*(aux2/gam2+aux2-1)
    return(crossprod(y - yM))
  }
  OF <- compiler::cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(UFR, -15, -30, -30 ,0 ,3),
                           upper = c (UFR, 30, 30, 30 ,3 ,6),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[5]
  gam2 <- matsout/betaV[6]
  aux1 <- 1-exp(-gam1)
  aux2 <- 1-exp(-gam2)  
  L <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
    betaV[4]*(aux2/gam2+aux2-1)
  aux3 <- 1 - aux1
  aux4 <- 1 - aux2
  Fwd <- UFR + betaV[2]*aux3 + betaV[3]*gam1*aux3 + betaV[4]*gam2*aux4
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=rep.int(0, J), T=matsout, L=L), fwd = Fwd))
  }
}
tZC_SV_extra <- compiler::cmpfun(tZC_SV_extra)
