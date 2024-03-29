\name{ycinter}
\alias{ycinter}
\title{Yield curve or zero-coupon prices interpolation}
\usage{
  ycinter(yM = NULL, p = NULL, matsin, matsout,
    method = c("NS", "SV", "SW", "HCSPL"),
    typeres = c("rates", "prices"))
}
\arguments{
  \item{yM}{A vector of non-negative numerical quantities,
  containing the yield to maturities.}

  \item{p}{A vector of non-negative numerical quantities,
  containing the zero-coupon prices.}

  \item{matsin}{A vector containing the observed
  maturities.}

  \item{matsout}{the output maturities needed.}

  \item{method}{A character string giving the type of
  method used fo intepolation. \code{method} can be either
  "NS" for Nelson-Siegel, "SV" for Svensson, "HCSPL" for
  Hermite cubic spline, or "SW" Smith-Wilson.}

  \item{typeres}{A character string, giving the type of
  return. Either "prices" or "rates".}
}
\value{
  An S4 Object
}
\description{
  Yield curve or zero-coupon bonds prices curve
  interpolation using the Nelson-Siegel , Svensson,
  Smith-Wilson models and an Hermite cubic spline.
}
\details{
  This function interpolates between observed points of a
  yield curve, or zero-coupon prices, using the
  Nelson-Siegel, Svensson, Smith-Wilson models and an
  Hermite cubic spline. The result can be either prices or
  zero rates.
}
\examples{
## Interpolation of yields to matuities with prices as outputs

 # Yield to maturities
 txZC <- c(0.01422,0.01309,0.01380,0.01549,0.01747,0.01940,0.02104,0.02236,0.02348,
 0.02446,0.02535,0.02614,0.02679,0.02727,0.02760,0.02779,0.02787,0.02786,0.02776
 ,0.02762,0.02745,0.02727,0.02707,0.02686,0.02663,0.02640,0.02618,0.02597,0.02578,0.02563)

 # Zero-coupon prices
 p <- c(0.9859794,0.9744879,0.9602458,0.9416551,0.9196671,0.8957363,0.8716268,0.8482628,
 0.8255457,0.8034710,0.7819525,0.7612204,0.7416912,0.7237042,0.7072136
 ,0.6922140,0.6785227,0.6660095,0.6546902,0.6441639,0.6343366,0.6250234,0.6162910,0.6080358,
 0.6003302,0.5929791,0.5858711,0.5789852,0.5722068,0.5653231)

 # Observed maturities
 u <- 1:30

 # Output maturities
 t <- seq(from = 1, to = 30, by = 0.5)

 # Cubic splines interpolation
 (yc <- ycinter(yM = txZC, matsin = u, matsout = t,
 method="HCSPL", typeres="rates"))

 # Nelson-Siegel interpolation
 (yc <- ycinter(yM = txZC, matsin = u, matsout = t,
 method="NS", typeres="prices"))

 # Svensson interpolation
 (yc <- ycinter(p = p, matsin = u, matsout = t,
 method="SV", typeres="prices"))

 #Smith-Wilson interpolation
 (yc <- ycinter(p = p, matsin = u, matsout = t,
 method="SW", typeres="rates"))
}
\author{
  Thierry Moudiki
}
\keyword{curve}
\keyword{interpolation,}
\keyword{yield}

