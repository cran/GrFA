\name{WPE}
\alias{WPE}
\title{
Weighted Projection Estimation
}
\description{
Weighted Projection Estimation
}
\usage{
WPE(y, rmax = 8, r0 = NULL, r = NULL, weight = TRUE, method = "ic1", type = "BIC3")
}
\arguments{
  \item{y}{a list of the observation data, each element is a data matrix of each group  with dimension \eqn{T * N_m}.}
  \item{rmax}{the maximum factor numbers of all groups.}
  \item{r0}{the number of global factors, default is \code{NULL}, the algorithm will automatically estimate the number of global factors. If you have the prior information about the true number of global factors, you can set it by your own.}
  \item{r}{the number of local factors in each group, default is \code{NULL}, the algorithm will automatically estimate the number of local factors. If you have the prior information about the true number of local factors, you can set it by your own, notice it should be an integer vector of length \eqn{M} (the number of groups).}
  \item{weight}{the weight of each projection matrix, default is \code{TRUE}, means \eqn{w_m = N_m/N}, if \code{weight = FALSE}, then simply calculate the mean of all projection matrices.}
  \item{method}{the method used in the algorithm, default is \code{ic1}, it can also be \code{ic2} and \code{gap}.}
  \item{type}{the method used in estimating the factor numbers in each group initially, default is \code{BIC3}, it can also be \code{IC3}}
}
\value{
  \item{r0hat}{the estimated number of the global factors.}
  \item{rho}{the estimated number of the local factors.}
  \item{Ghat}{the estimated global factors.}
  \item{Fhat}{the estimated local factors.}
  \item{loading_G}{a list consisting of the estimated global factor loadings.}
  \item{loading_F}{a list consisting of the estimated local factor loadings.}
  \item{e}{a list consisting of the residuals.}
  \item{threshold}{the threshold used in determining the number of global factors, only for \code{method = ic1} and \code{method = ic2}.}
}
\author{
Jiaqi Hu
}
\examples{
dat = gendata()
dat
WPE(dat$y, method = "ic1")
WPE(dat$y, method = "ic2")
WPE(dat$y, method = "gap")
}
