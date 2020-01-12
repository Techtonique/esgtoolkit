
# Misc plots --------------------------------------------------------------


#'@title Plot time series percentiles and confidence intervals
#'
#'@description
#'
#'This function plots colored bands for time series percentiles and confidence 
#'intervals. You can use it for outputs from \code{link{simdiff}}, 
#'\code{link{esgmartingaletest}}, \code{link{esgcortest}}.
#'
#'@param x a times series object 
#'
#'@param ... additionnal (optional) parameters provided to \code{plot}
#'
#'@export
#'
#'@author Thierry Moudiki
#'
#'@seealso \code{\link{esgplotts}}
#'
#'@examples
#'
#'# Times series
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma
#'x <- simdiff(n = 100, horizon = 5, 
#'frequency = "quart",  
#'model = "OU", 
#'x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'par(mfrow=c(2,1))
#'esgplotbands(x, xlab = "time", ylab = "values")
#'matplot(time(x), x, type = 'l', xlab = "time", ylab = "series values")
#'
#'# Martingale test
#'
#'r0 <- 0.03
#'S0 <- 100
#'sigma0 <- 0.1
#'nbScenarios <- 100
#'horizon0 <- 10
#'eps0 <- simshocks(n = nbScenarios, horizon = horizon0, frequency = "quart",
#' method = "anti")
#'sim.GBM <- simdiff(n = nbScenarios, horizon = horizon0, frequency = "quart",   
#'                  model = "GBM", 
#'                  x0 = S0, theta1 = r0, theta2 = sigma0, 
#'                  eps = eps0)
#'
#'mc.test <- esgmartingaletest(r = r0, X = sim.GBM, p0 = S0, alpha = 0.05)   
#'esgplotbands(mc.test)
#'
#'# Correlation test
#'
#'nb <- 500
#'
#'s0.par1 <- simshocks(n = nb, horizon = 3, frequency = "semi",
#'family = 1, par = 0.2)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 3, frequency = "semi", 
#'family = 1, par = 0.8)
#'
#'(test1 <- esgcortest(s0.par1))
#'(test2 <- esgcortest(s0.par2))
#'par(mfrow=c(2, 1))
#'esgplotbands(test1)
#'esgplotbands(test2)
esgplotbands <- function(x, ...)
{
  if (is.ts(x))
  {    
    nrow.x <- nrow(x)      
    x0 <- x[1, 1]
    x1 <- (as.numeric(x[nrow.x, 1]))
    x2 <- (as.numeric(x[nrow.x, 2]))
    x3 <- (as.numeric(x[nrow.x, 3]))
    
    cond1 <- (x1 == x2)
    cond2 <- (x2 == x3)
    
    if(cond1 && cond2) 
    {
      x_up <- x[-c(1, nrow.x), ]
    }
    else
    {
      x_up <- x[-1, ]
    }
    
    qt.95 <- c(x0, apply(x_up, 1, function(x) quantile(x, 0.95)))
    qt.05 <- c(x0, apply(x_up, 1, function(x) quantile(x, 0.05)))
    
    x.summary <- cbind(rep(x0, 5), 
                       apply(x_up, 1, function(x) summary(x))[-3, ])
    x.ci <- cbind(rep(x0, 3), 
                  rbind(apply(x_up, 1, function(x) t.test(x)$conf.int)[1, ],
                        rowMeans(x_up),
                        apply(x_up, 1, function(x) t.test(x)$conf.int)[2, ]))
    jet.colors <- colorRampPalette( c("lightyellow", "lightgreen") )
    nbcol <- 3
    color <- jet.colors(nbcol)
    
    if(cond1 && cond2) 
    {
      abs <- as.numeric(time(x))[-(nrow.x - 1)]
    }
    else{
      abs <- as.numeric(time(x))
    }
    
    y.mean <- x.summary[3, ]
    
    bands.plot(abs, y.mean, ci.upper = x.summary[1, ], ci.lower = x.summary[5, ], 
               col = color[1], ylim = c(min(x.summary[1,]), max(x.summary[5,])), ...)
    bands.add(abs, y.mean, col = color[2], ci.upper = qt.95, 
              ci.lower = qt.05)
    bands.add(abs, y.mean, col = color[3], ci.upper = x.summary[2, ], 
              ci.lower = x.summary[4, ])  
  }
  
  if (is.list(x))
  {
    is.cor.test <- is.numeric(try(x$cor.estimate, silent = TRUE))
    
    if (is.cor.test)
    {
      conf.int <- x$conf.int
      abs <- as.numeric(time(conf.int))
      bands.plot(abs, x$cor.estimate, ci.upper = conf.int[ , 2], ci.lower = conf.int[ , 1], 
                 col = "#C7F6B8", ylim = c(min(conf.int[ , 1]), max(conf.int[ , 2])), 
                 main = "conf. int for the correlations", xlab = "time", ylab = "conf. int.", ...)
      points(abs, x$cor.estimate, pch = 16)
    }
    
    if (!is.cor.test)
    {
      abs <- as.numeric(time(x$conf.int))
      y.mean <- rep(0, length(abs))
      
      par(mfrow = c(2, 1))
      bands.plot(abs, y.mean, ci.upper = x$conf.int[, 2], ci.lower =  x$conf.int[, 1], 
                 col = "gray80", ylim = c(min(x$conf.int[, 1]), max(x$conf.int[, 2])), 
                 xlab = "time", ylab = "conf. int.", 
                 main = "conf. int. \n for the martingale difference", ...)
      lines(abs, y.mean, col = "blue", lty = 2)
      
      plot(abs, x$true_prices, type = 'l', col = "black", 
           ylim = x$true_prices[1]*c(1 - 0.03, 1 + 0.03), 
           main = "true (black) vs \n monte carlo (blue) prices", 
           xlab = "time", ylab = "prices")
      lines(abs, x$mc.prices, col = "blue")
      points(abs, x$true_prices, col = "black", pch = 16)
      points(abs, x$mc.prices,  col = "blue",  pch = 16)
    }
  }
}



#'@title Plot time series objects
#'
#'@description This function plots outputs from \code{\link{simdiff}}.
#'
#'@details For a large number of simulations, it's preferable to use 
#'\code{\link{esgplotbands}} for a synthetic view by percentiles.
#'
#'@param x a time series object, an output from \code{\link{simdiff}}.
#'
#'@seealso \code{\link{simdiff}}, \code{\link{esgplotbands}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@references
#'
#'H. Wickham (2009), ggplot2: elegant graphics for data analysis. Springer 
#'New York.
#'
#'@examples
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma
#'x <- simdiff(n = 10, horizon = 5, frequency = "quart",  
#'model = "OU", 
#'x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'esgplotts(x)
#'
esgplotts <- function(x)
{ 
  x0 <- start(x)[1]
  dt_x <- deltat(x)
  tt <- cumsum(c(x0, rep.int(dt_x, dim(x)[1]-1)))
  meltdf <- melt(cbind(tt, as.data.frame(x)), id="tt")
  group <- rep(LETTERS[1:dim(x)[2]], each = dim(x)[1])
  qplot(tt, meltdf$value, data = meltdf, geom = "line", colour=group) + 
    xlab("Maturity") + ylab("Values") + 
    ggplot2::theme(legend.position = "none")
}
esgplotts <- cmpfun(esgplotts)



#'@title Visualize the dependence between 2 gaussian shocks
#'
#'@description
#'
#'This function helps you in visualizing the dependence between 2 gaussian
#' shocks. 
#'
#'@param x an output from \code{\link{simshocks}}, a list with 2 components.
#'
#'@param y an output from \code{\link{simshocks}}, a list with 2 components 
#'(Optional). 
#'
#'@export
#'
#'@seealso \code{\link{simshocks}}
#'
#'@author Thierry Moudiki + some nice blogs :) 
#'
#'@references
#'
#'H. Wickham (2009), ggplot2: elegant graphics for data analysis. Springer 
#'New York.
#'
#'@examples
#'
#'# Number of risk factors
#'d <- 2
#'
#'# Number of possible combinations of the risk factors
#'dd <- d*(d-1)/2
#'
#'# Family : Gaussian copula 
#'fam1 <- rep(1,dd)
#'# Correlation coefficients between the risk factors (d*(d-1)/2)
#'par0.1 <- 0.1
#'par0.2 <- -0.9
#'
#'# Family : Rotated Clayton (180 degrees)
#'fam2 <- 13
#'par0.3 <- 2
#'
#'# Family : Rotated Clayton (90 degrees)
#'fam3 <- 23
#'par0.4 <- -2
#'
#'# number of simulations
#'nb <- 500
#'
#'# Simulation of shocks for the d risk factors
#'s0.par1 <- simshocks(n = nb, horizon = 4, 
#'family = fam1, par = par0.1)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 4, 
#'family = fam1, par = par0.2)
#'
#'s0.par3 <- simshocks(n = nb, horizon = 4, 
#'family = fam2, par = par0.3)
#'
#'s0.par4 <- simshocks(n = nb, horizon = 4, 
#'family = fam3, par = par0.4)
#'
#'\dontrun{
#'esgplotshocks(s0.par1, s0.par2)
#'esgplotshocks(s0.par2, s0.par3)
#'esgplotshocks(s0.par2, s0.par4)
#'esgplotshocks(s0.par1, s0.par4)} 
#'
esgplotshocks <-  function(x, y = NULL)
{
  x <- matrix(unlist(x), ncol = 2)
  if (!is.null(y))
  {
    y <- matrix(unlist(y), ncol = 2)
    if (length(x) != length(y)) stop("'x' and 'y' must have the same dimensions")
    xvar <- c(x[,1], y[,1])
    yvar <- c(x[,2], y[,2])
    nb  <- dim(x)[1]
    zvar <- as.factor(c(rep("x", nb), rep("y", nb)))
    xy <- data.frame(xvar, yvar, zvar)
  }
  else 
  {
    xvar <- x[,1]
    yvar <- x[,2]
    nb  <- dim(x)[1]
    zvar <- as.factor(rep("x", nb))
    xy <- data.frame(xvar, yvar, zvar)
  }
  
  #placeholder plot - prints nothing at all
  empty <- ggplot()+ggplot2::geom_point(ggplot2::aes(1,1), colour="white") +
    ggplot2::theme(                              
      plot.background = ggplot2::element_blank(), 
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(), 
      panel.border = ggplot2::element_blank(), 
      panel.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  #scatterplot of x and y variables
  scatter <- ggplot(xy, ggplot2::aes(xvar, yvar)) + 
    ggplot2::geom_point(ggplot2::aes(color=zvar)) + 
    ggplot2::scale_color_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position=c(1,1),legend.justification=c(1,1)) 
  
  #marginal density of x - plot on top
  plot_top <- ggplot(xy, ggplot2::aes(xvar, fill=zvar)) + 
    ggplot2::geom_density(alpha=.5) + 
    ggplot2::scale_fill_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position = "none")
  
  #marginal density of y - plot on the right
  plot_right <- ggplot(xy, ggplot2::aes(yvar, fill=zvar)) + 
    ggplot2::geom_density(alpha=.5) + 
    ggplot2::coord_flip() + 
    ggplot2::scale_fill_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position = "none") 
  
  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}
esgplotshocks  <- cmpfun(esgplotshocks)

