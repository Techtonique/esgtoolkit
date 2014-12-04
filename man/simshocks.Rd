\name{simshocks}
\alias{simshocks}
\title{Underlying gaussian shocks for risk factors' simulation.}
\usage{
  simshocks(n, horizon,
    frequency = c("annual", "semi-annual", "quarterly", "monthly", "weekly", "daily"),
    method = c("classic", "antithetic", "mm", "hybridantimm", "TAG"),
    family = NULL, par = NULL, par2 = NULL,
    type = c("CVine", "DVine"))
}
\arguments{
  \item{n}{number of independent observations for each risk
  factor.}

  \item{horizon}{horizon of projection.}

  \item{frequency}{either "annual", "semi-annual",
  "quarterly", "monthly", "weekly", or "daily" (1, 1/2,
  1/4, 1/12, 1/52, 1/252).}

  \item{method}{either classic monte carlo, antithetic
  variates, moment matching, hybrid antithetic variates +
  moment matching or "TAG" (see the 4th reference for the
  latter).}

  \item{family}{the same as \code{\link{CDVineSim}} from
  package \code{CDVine}. A d*(d-1)/2 integer vector of
  C-/D-vine pair-copula families with values 0 =
  independence copula, 1 = Gaussian copula, 2 = Student t
  copula (t-copula), 3 = Clayton copula, 4 = Gumbel copula,
  5 = Frank copula, 6 = Joe copula, 7 = BB1 copula, 8 = BB6
  copula, 9 = BB7 copula, 10 = BB8 copula, 13 = rotated
  Clayton copula (180 degrees; "survival Clayton"), 14 =
  rotated Gumbel copula (180 degrees; "survival Gumbel"),
  16 = rotated Joe copula (180 degrees; "survival Joe"), 17
  = rotated BB1 copula (180 degrees; "survival BB1"), 18 =
  rotated BB6 copula (180 degrees; "survival BB6"), 19 =
  rotated BB7 copula (180 degrees; "survival BB7"), 20 =
  rotated BB8 copula (180 degrees; "survival BB8"), 23 =
  rotated Clayton copula (90 degrees), 24 = rotated Gumbel
  copula (90 degrees), 26 = rotated Joe copula (90
  degrees), 27 = rotated BB1 copula (90 degrees), 28 =
  rotated BB6 copula (90 degrees), 29 = rotated BB7 copula
  (90 degrees), 30 = rotated BB8 copula (90 degrees), 33 =
  rotated Clayton copula (270 degrees), 34 = rotated Gumbel
  copula (270 degrees), 36 = rotated Joe copula (270
  degrees), 37 = rotated BB1 copula (270 degrees), 38 =
  rotated BB6 copula (270 degrees), 39 = rotated BB7 copula
  (270 degrees), 40 = rotated BB8 copula (270 degrees)}

  \item{par}{the same as \code{\link{CDVineSim}} from
  package \code{CDVine}. A d*(d-1)/2 vector of pair-copula
  parameters.}

  \item{par2}{the same as \code{\link{CDVineSim}} from
  package \code{CDVine}. A d*(d-1)/2 vector of second
  parameters for pair-copula families with two parameters
  (t, BB1, BB6, BB7, BB8; no default).}

  \item{type}{type of the vine model: 1 : C-vine 2 :
  D-vine}
}
\value{
  If \code{family} and \code{par} are not provided, a
  univariate time series object with simulated gaussian
  shocks for one risk factor. Otherwise, a list of time
  series objects, containing gaussian shocks for each risk
  factor.
}
\description{
  This function makes simulations of correlated or
  dependent gaussian shocks for risk factors.
}
\details{
  The function shall be used along with
  \code{\link{simdiff}}, in order to embed correlated or
  dependent random gaussian shocks into simulated
  diffusions. \code{\link{esgplotshocks}} can help in
  visualizing the type of dependence between the shocks.
}
\examples{
# Number of risk factors
d <- 6

# Number of possible combinations of the risk factors
dd <- d*(d-1)/2

# Family : Gaussian copula for all
fam1 <- rep(1,dd)

# Correlation coefficients between the risk factors (d*(d-1)/2)
par1 <- c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,
        -0.2,-0.32,-0.19,-0.17,-0.06)


# Simulation of shocks for the 6 risk factors
simshocks(n = 10, horizon = 5, family = fam1, par = par1)


# Simulation of shocks for the 6 risk factors
# on a quarterly basis
simshocks(n = 10, frequency = "quarterly", horizon = 2, family = fam1,
par = par1)


# Simulation of shocks for the 6 risk factors simulation
# on a quarterly basis, with antithetic variates and moment matching.
s0 <- simshocks(n = 10, method = "hyb", horizon = 4,
family = fam1, par = par1)


s0[[2]]
colMeans(s0[[1]])
colMeans(s0[[5]])
apply(s0[[3]], 2, sd)
apply(s0[[4]], 2, sd)
}
\author{
  Thierry Moudiki
}
\references{
  Brechmann, E., Schepsmeier, U. (2013). Modeling
  Dependence with C- and D-Vine Copulas: The R Package
  CDVine. Journal of Statistical Software, 52(3), 1-27. URL
  \url{http://www.jstatsoft.org/v52/i03/}.

  Genz, A. Bretz, F., Miwa, T. Mi, X., Leisch, F., Scheipl,
  F., Hothorn, T. (2013). mvtnorm: Multivariate Normal and
  t Distributions. R package version 0.9-9996.

  Genz, A. Bretz, F. (2009), Computation of Multivariate
  Normal and t Probabilities. Lecture Notes in Statistics,
  Vol. 195., Springer-Verlag, Heidelberg. ISBN
  978-3-642-01688-2.

  Nteukam T, O., & Planchet, F. (2012). Stochastic
  evaluation of life insurance contracts: Model point on
  asset trajectories and measurement of the error related
  to aggregation. Insurance: Mathematics and Economics,
  51(3), 624-631. URL
  \url{http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/ab539dcebcc4e77ac12576c6004afa67/$FILE/Article_US_v1.5.pdf}
}
\seealso{
  \code{\link{simdiff}}, \code{\link{esgplotshocks}}
}

