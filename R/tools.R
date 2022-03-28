
# Tools -------------------------------------------------------------------

# quantiles
quantileESG <- function (x, probs) 
{
  eps <- 100 * .Machine$double.eps
  if (any((p.ok <- !is.na(probs)) & (probs < -eps | probs > 
                                     1 + eps))) 
    stop("'probs' outside [0,1]")
  n <- length(x)
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
    probs <- pmax(0, pmin(1, probs))
  }
  np <- length(probs)
  if (n > 0 && np > 0) {
    index <- 1 + (n - 1) * probs
    lo <- floor(index)
    hi <- ceiling(index)
    x <- sort(x, partial = unique(c(lo, hi)))
    qs <- x[lo]
    i <- which(index > lo)
    h <- (index - lo)[i]
    qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
  }
  else {
    qs <- rep(NA_real_, np)
  }
  qs
}
quantileESG <- compiler::cmpfun(quantileESG)

# normal distrib. with mean = 0, sd = 1 
rnormESG <- function(n, m = NULL)
{
  if (m == 1 || is.null(m))
  {
    return(as.vector(rnormESGcpp(N = n, M = 1)))
  }
  else 
  {
    return(rnormESGcpp(N = n, M = m))
  }
}
rnormESG <- compiler::cmpfun(rnormESG)

# simulation with TAG
TAG <- function(n, m) 
{
  TAGbase <- function(n)
  {
    n2 <- 10000
    sim <- rnormESG(n = n2)
    sj <- quantileESG(sim, (0:n)/n)
    sj_up <- sj[-1]
    sj_down <- sj[-(n+1)]
    out <- TAGcorecpp(sim = sim, sj_down = sj_down, 
                      sj_up = sj_up, n = n2, p = n) 
    sample(out)
  }
  
  if (m == 1 || is.null(m))
  {
    return(TAGbase(n))
  }
  else
  {
    return(t(replicate(m, TAGbase(n))))
  }
}
TAG <- compiler::cmpfun(TAG)

# scaling a matrix
scaleESG <- function (x, center = TRUE, scale = TRUE) 
{
  #  x <- as.matrix(x)
  nc <- ncol(x)
  if (is.logical(center)) {
    if (center) {
      center <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, center, check.margin = FALSE)
    }
  }
  else if (is.numeric(center) && (length(center) == nc)) 
    x <- sweep(x, 2L, center, check.margin = FALSE)
  else stop("length of 'center' must equal the number of columns of 'x'")
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/max(1, length(v) - 1L))
      }
      scale <- apply(x, 2L, f)
      x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
    }
  }
  else if (is.numeric(scale) && length(scale) == nc) 
    x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  else stop("length of 'scale' must equal the number of columns of 'x'")
  return(x)
}
scaleESG <- compiler::cmpfun(scaleESG)

# plot bands
bands.plot <- function(x, y.mean, ci.upper, ci.lower, col, y.goal = NULL, goal.col = "blue", ...)
{
  if (missing(x)) stop("'x' must be provided")
  if (missing(y.mean)) stop("'x' must be provided")
  if (missing(ci.upper) & missing(ci.lower)) stop("'ci.upper' and 'ci.lower' must be provided")
  
  plot(x = x, y = y.mean, type = 'l', ...)
  polygon(c(x, rev(x)), 
          c(ci.upper, rev(ci.lower)), 
          col = col, border = FALSE)
  lines(x, y.mean, lwd = 2)  
  if (!is.null(y.goal))
  {
    abline(h = y.goal, lty = 2, lwd = 2, col = goal.col)
  }
}

# add bands on a plot
bands.add <- function(x, y.mean, col, ci.upper, ci.lower)
{
  if (missing(x)) stop("'x' must be provided")
  if (missing(col)) stop("'col' must be provided")
  
  polygon(c(x, rev(x)), 
          c(ci.upper, rev(ci.lower)), 
          col = col, border = FALSE)
  lines(x, y.mean, lwd = 2)
}
