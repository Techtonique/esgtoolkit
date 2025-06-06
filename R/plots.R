
# Misc plots --------------------------------------------------------------


#' Plot time series percentiles and confidence intervals
#' 
#' @description
#' Creates a plot showing time series data with confidence intervals, percentiles, and mean values.
#' 
#' @param x A time series object containing the data to plot
#' @param ... Additional arguments passed to the plotting function
#' 
#' @details
#' The function creates a plot with:
#' \itemize{
#'   \item Mean values as the central line
#'   \item 95th and 5th percentiles as outer bands
#'   \item Upper and lower quartiles as inner bands
#'   \item Confidence intervals based on t-tests
#' }
#' 
#' The plot uses a color gradient from light yellow to light green for the different bands.
#' 
#' @return A plot object
#' 
#' @examples
#' # Create sample time series
#' x <- ts(matrix(rnorm(1000), ncol=10), start=0, deltat=1)
#' 
#' # Plot with default settings
#' esgplotbands(x)
#' 
#' # Plot with custom title
#' esgplotbands(x, main="Custom Title")
#' 
#' @export
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
      
      par(mfrow = c(1, 2))
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


#' Plot time series objects
#' 
#' @description
#' Creates a line plot of time series data using ggplot2, with each series represented by a different color.
#' 
#' @param x A time series object (ts) containing the data to plot
#' 
#' @details
#' The function converts the time series data into a format suitable for ggplot2 plotting.
#' Each column of the time series is plotted as a separate line with a different color.
#' The x-axis represents time/maturity and the y-axis shows the values.
#' 
#' @return A ggplot2 object containing the line plot
#' 
#' @examples
#' # Create sample time series data
#' x <- ts(matrix(rnorm(100), 20, 5), start = 0, deltat = 0.1)
#' 
#' # Plot the time series
#' esgplotts(x)
#' 
#' @export
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


#' Visualize the dependence between Gaussian shocks
#' 
#' @description
#' Creates a scatter plot with marginal density plots to visualize the relationship between two sets of Gaussian shocks.
#' 
#' @param x A matrix or list containing two columns of Gaussian shocks
#' @param y Optional second matrix or list of Gaussian shocks to compare with x
#' 
#' @details
#' The function creates a scatter plot of the shocks with marginal density plots on the top and right sides.
#' If y is provided, both sets of shocks are plotted with different colors (blue for x, red for y).
#' The marginal density plots show the distribution of shocks along each dimension.
#' 
#' @return A ggplot2 object containing the scatter plot with marginal densities
#' 
#' @examples
#' # Generate sample Gaussian shocks
#' x <- matrix(rnorm(1000), 500, 2)
#' 
#' # Plot single set of shocks
#' esgplotshocks(x)
#' 
#' # Plot two sets of shocks for comparison
#' y <- matrix(rnorm(1000), 500, 2)
#' esgplotshocks(x, y)
#' 
#' @export
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
  empty <- ggplot2::ggplot()+ggplot2::geom_point(ggplot2::aes(1,1), colour="white") +
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
  scatter <- ggplot2::ggplot(xy, ggplot2::aes(xvar, yvar)) + 
    ggplot2::geom_point(ggplot2::aes(color=zvar)) + 
    ggplot2::scale_color_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position=c(1,1),legend.justification=c(1,1)) 
  
  #marginal density of x - plot on top
  plot_top <- ggplot2::ggplot(xy, ggplot2::aes(xvar, fill=zvar)) + 
    ggplot2::geom_density(alpha=.5) + 
    ggplot2::scale_fill_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position = "none")
  
  #marginal density of y - plot on the right
  plot_right <- ggplot2::ggplot(xy, ggplot2::aes(yvar, fill=zvar)) + 
    ggplot2::geom_density(alpha=.5) + 
    ggplot2::coord_flip() + 
    ggplot2::scale_fill_manual(values = c("blue", "red")) + 
    ggplot2::theme(legend.position = "none") 
  
  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(plot_top, empty, scatter, plot_right, 
               ncol=2, nrow=2, widths=c(4, 1), 
               heights=c(1, 4))
}