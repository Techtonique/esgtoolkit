\name{esgfwdrates}
\alias{esgfwdrates}
\title{Instantaneous forward rates}
\usage{
  esgfwdrates(in.maturities, in.zerorates, n, horizon,
    out.frequency = c("annual", "semi-annual", "quarterly", "monthly", "weekly", "daily"),
    ...)
}
\arguments{
  \item{in.maturities}{input maturities}

  \item{in.zerorates}{input zero rates}

  \item{n}{number of independent observations}

  \item{horizon}{horizon of projection}

  \item{out.frequency}{either "annual", "semi-annual",
  "quarterly", "monthly", "weekly", or "daily" (1, 1/2,
  1/4, 1/12, 1/52, 1/252)}

  \item{...}{additional parameters provided to
  \code{\link{ycinter}}}
}
\description{
  This function provides instantaneous forward rates. They
  can be used in no-arbitrage short rate models, to fit the
  yield curve exactly.
}
\examples{
# Yield to maturities
txZC <- c(0.01422,0.01309,0.01380,0.01549,0.01747,0.01940,0.02104,0.02236,0.02348,
         0.02446,0.02535,0.02614,0.02679,0.02727,0.02760,0.02779,0.02787,0.02786,0.02776
         ,0.02762,0.02745,0.02727,0.02707,0.02686,0.02663,0.02640,0.02618,0.02597,0.02578,0.02563)

# Observed maturities
u <- 1:30

\dontrun{
par(mfrow=c(2,2))
fwdNS <- esgfwdrates(in.maturities = u, in.zerorates = txZC,
                    n = 10, horizon = 20,
                    out.frequency = "semi-annual", method = "NS")
matplot(time(fwdNS), fwdNS, type = 'l')

fwdSV <- esgfwdrates(in.maturities = u, in.zerorates = txZC,
                    n = 10, horizon = 20,
                    out.frequency = "semi-annual", method = "SV")
matplot(time(fwdSV), fwdSV, type = 'l')

fwdSW <- esgfwdrates(in.maturities = u, in.zerorates = txZC,
                    n = 10, horizon = 20,
                    out.frequency = "semi-annual", method = "SW")
matplot(time(fwdSW), fwdSW, type = 'l')

fwdHCSPL <- esgfwdrates(in.maturities = u, in.zerorates = txZC,
                       n = 10, horizon = 20,
                       out.frequency = "semi-annual", method = "HCSPL")
matplot(time(fwdHCSPL), fwdHCSPL, type = 'l')
}
}
\author{
  Thierry Moudiki
}
\references{
  Thierry Moudiki (2013). ycinterextra: Yield curve or
  zero-coupon prices interpolation and extrapolation. R
  package version 0.1. URL
  \url{http://CRAN.R-project.org/package=ycinterextra}
}

