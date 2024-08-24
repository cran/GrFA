\name{CCA}
\alias{CCA}
\title{
Canonical Correlation Estimation
}
\description{
Canonical Correlation Estimation
}
\usage{
CCA(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, method = "CCD", type = "IC3")
}
\arguments{
  \item{y}{a list of the observation data, each element is a data matrix of each group  with dimension \eqn{T * N_m}.}
  \item{rmax}{the maximum factor numbers of all groups.}
  \item{r0}{the number of global factors, default is \code{NULL}, the algorithm will automatically estimate the number of global factors. If you have the prior information about the true number of global factors, you can set it by your own.}
  \item{r}{the number of local factors in each group, default is \code{NULL}, the algorithm will automatically estimate the number of local factors. If you have the prior information about the true number of local factors, you can set it by your own, notice it should be an integer vector of length \eqn{M} (the number of groups).}
  \item{localfactor}{if \code{localfactor = FALSE}, then we would not estimate the local factors; if \code{localfactor = TRUE}, then we will further estimate the local factors.}
  \item{method}{the method used in the algorithm, default is \code{CCD}, it can also be \code{MCC}.}
  \item{type}{the method used in estimating the factor numbers in each group initially, default is \code{IC3}.}
}
\value{
  \item{r0hat}{the estimated number of the global factors.}
  \item{rho}{the estimated number of the local factors.}
  \item{Ghat}{the estimated global factors.}
  \item{Fhat}{the estimated local factors.}
  \item{loading_G}{a list consisting of the estimated global factor loadings.}
  \item{loading_F}{a list consisting of the estimated local factor loadings.}
  \item{e}{a list consisting of the residuals.}
  \item{threshold}{the threshold used in determining the number of global factors, only for \code{method = "MCC"}.}
}
\references{
Choi, I., Lin, R., & Shin, Y. (2021). Canonical correlation-based model selection for the multilevel factors. Journal of Econometrics.
}
\examples{
dat = gendata()
dat
CCA(dat$y, rmax = 8, localfactor = TRUE, method = "CCD")
CCA(dat$y, rmax = 8, localfactor = TRUE, method = "MCC")
}

