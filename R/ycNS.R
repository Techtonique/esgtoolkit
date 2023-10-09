

# Nelson-Siegel interpolation
tZC_NS <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
  #   require(randtoolbox)
  #   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=rep.int(0, N), T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[4]
    aux1 <- 1-exp(-gam1)     
    y <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
    return(crossprod(y - yM))
  }
  OF <- compiler::cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(0, -15, -30, 0),
                           upper = c (15, 30, 30, 3),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[4]
  aux1 <- 1-exp(-gam1)
  L <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
  aux2 <- 1 - aux1
  Fwd <- betaV[1] + betaV[2]*aux2 + betaV[3]*gam1*aux2
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=0, T=matsout, L=L), fwd = Fwd))}
}
tZC_NS <- compiler::cmpfun(tZC_NS)

# Nelson-Siegel extrapolation
tZC_NS_extra <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"), UFR)
{
  #   require(randtoolbox)
  #   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=rep.int(0, N), T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[4]
    aux1 <- 1-exp(-gam1)     
    y <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
    return(crossprod(y - yM))
  }
  OF <- compiler::cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(UFR, -15, -30, 0),
                           upper = c (UFR, 30, 30, 3),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[4]
  aux1 <- 1-exp(-gam1)
  L <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
  aux2 <- 1 - aux1
  Fwd <- UFR + betaV[2]*aux2 + betaV[3]*gam1*aux2
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=0, T=matsout, L=L), fwd = Fwd))
  }
}
tZC_NS_extra <- compiler::cmpfun(tZC_NS_extra)
