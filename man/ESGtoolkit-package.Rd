\name{ESGtoolkit-package}
\alias{ESGtoolkit-package}
\alias{ESGtoolkit}
\docType{package}
\title{
Toolkit for financial assets and interest rates simulation. 
}
\description{
Toolkit for Monte Carlo simulation of financial assets and interest rates, involved in an Economic Scenario Generator (ESG). 
}
\details{
\tabular{ll}{
Package: \tab ESGtoolkit\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-06-13\cr
License: \tab GPL-2 | GPL-3\cr
}

The main functions of the package are : \cr

- \code{\link{simdiff}} for the simulation of diffusion processes.\cr
- \code{\link{simshocks}} for the custom simulation of the gaussian shocks embedded into the diffusion processes. 

There are also several functions for statistical tests on the simulations, and 
for visualization. 
}
\author{
Jean-Charles Croix, Thierry Moudiki, Frederic Planchet, Wassim Youssef

Maintainer: Thierry Moudiki <thierry.moudiki@gmail.com>
}
\references{
Auguie, B. (2012). gridExtra: functions in Grid graphics. R package version
  0.9.1. URL \url{http://CRAN.R-project.org/package=gridExtra}

D. J. Best & D. E. Roberts (1975), Algorithm AS 89: The Upper Tail 
Probabilities of Spearman's rho. Applied Statistics, 24, 377-379.

Black, F., Scholes, M.S. (1973) The pricing of options and corporate liabilities, 
Journal of Political Economy, 81, 637-654.

Brechmann, E., Schepsmeier, U. (2013). Modeling Dependence with C-
and D-Vine Copulas: The R Package CDVine. Journal of Statistical Software,
52(3), 1-27. URL \url{http://www.jstatsoft.org/v52/i03/}.

Cox, J.C., Ingersoll, J.E., Ross, S.A. (1985) A theory of the term structure of
 interest rates, Econometrica, 53, 385-408.

Genz, A. Bretz, F., Miwa, T. Mi, X., Leisch, F., Scheipl, F., Hothorn, T. (2013).
 mvtnorm: Multivariate Normal and t Distributions. R package version 0.9-9996. 
 URL \url{http://CRAN.R-project.org/package=mvtnorm}

Genz, A. Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. 
Lecture Notes in Statistics, Vol. 195., Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.

Glasserman, P. (2004). Monte Carlo methods in financial engineering (Vol. 53). 
Springer.

Hollander, M. & D. A. Wolfe (1973), Nonparametric Statistical Methods.
 New York: John Wiley & Sons. Pages 185-194 (Kendall and Spearman tests).

Iacus, S. M. (2009). Simulation and inference for stochastic differential 
equations: with R examples (Vol. 1). Springer.

Kou S, (2002), A jump diffusion model for option pricing, Management 
Science Vol. 48, 1086-1101.

Merton, R. C. (1976). Option pricing when underlying stock returns are 
discontinuous. Journal of financial economics, 3(1), 125-144.

Moudiki, T. (2013). ycinterextra: Yield curve or zero-coupon prices
interpolation and extrapolation. R package version 0.1. URL \url{http://CRAN.R-project.org/package=ycinterextra}

Nteukam T, O., & Planchet, F. (2012). Stochastic evaluation of life insurance 
contracts: Model point on asset trajectories and measurement of the error 
related to aggregation. Insurance: Mathematics and Economics, 51(3), 624-631.
URL \url{http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/ab539dcebcc4e77ac12576c6004afa67/$FILE/Article_US_v1.5.pdf}

Uhlenbeck, G. E., Ornstein, L. S. (1930) On the theory of Brownian motion, 
Phys. Rev., 36, 823-841.

Vasicek, O. (1977) An Equilibrium Characterization of the Term Structure, 
Journal of Financial Economics, 5, 177-188.

Wickham, H. (2009), ggplot2: elegant graphics for data analysis. Springer 
New York. Package URL \url{http://CRAN.R-project.org/package=ggplot2}
}

\keyword{ESG, Economic Scenario Generator, Finance, Insurance, Risk Management}
