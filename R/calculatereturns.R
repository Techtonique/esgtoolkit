
#' Calculate returns or log-returns for multivariate time series
#'
#' @param x Multivariate time series 
#' @param type Type of return: basic return ("basic") or log-return ("log")
#'
#' @return
#' @export
#'
#' @examples
#' 
#' kappa <- 1.5
#' V0 <- theta <- 0.04
#' sigma_v <- 0.2
#' theta1 <- kappa*theta
#' theta2 <- kappa
#' n_simulations <- 10
#' 
#' eps0 <- simshocks(n = n_simulations, horizon = 5, frequency = "quart")
#' sim_GBM <- simdiff(n = n_simulations, horizon = 5, frequency = "quart",   
#'                    model = "GBM", 
#'                 x0 = 100, theta1 = 0.03, theta2 = 0.1, eps = eps0) 
#' 
#' returns <- calculatereturns(sim_GBM)
#' log_returns <- diff(log(sim_GBM))
#' 
#' print(identical(returns, log_returns))
#' 
#' par(mfrow=c(1, 2))
#' matplot(returns, type = 'l')
#' matplot(log_returns, type = 'l')
#' 
calculatereturns <- function(x, type = c("basic", "log"))
{
  # ?rlang::abort #?
  stopifnot(is.ts(x))
  names_x <- colnames(x)
  type <- match.arg(type)
  
  if (identical(type, "basic"))
  {
    res <- diff(x)/lag(x)
    colnames(res) <- names_x
    return (res)
  } 
  res <- diff(log(x))
  colnames(res) <- names_x
  return (res)
}