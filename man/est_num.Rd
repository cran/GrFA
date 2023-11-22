\name{est_num}
\alias{est_num}
\title{
Estimate factor numbers
}
\description{
Estimate factor numbers.
}
\usage{
est_num(X, kmax = 8, type = "BIC3")
}
\arguments{
  \item{X}{the observation data matrix of dimension \eqn{T\times N}.}
  \item{kmax}{the maximum number of factors.}
  \item{type}{the criterion used in determining the number of factors, default is \code{type = "BIC3"}, it can also be \code{"PC1", "PC2", "PC3", "IC1", "IC2","IC3", "AIC3", "BIC3", "ER", "GR"}.}
}
\value{
  \item{rhat}{the estimated number of factors.}
}
\references{
Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.

Ahn, S. C., & Horenstein, A. R. (2013). Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203-1227.
}

