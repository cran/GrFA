\name{FA}
\alias{FA}
\title{
Classical factor analysis
}
\description{
Classical factor analysis.
}
\usage{
FA(X, r)
}
\arguments{
  \item{X}{the observation data matrix of dimension \eqn{T\times N}.}
  \item{r}{the factor numbers need to estimated.}
}

\value{
  \item{F}{the estimated factors.}
  \item{L}{the estimated factor loadings.}
}
\references{
Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
}
\author{
Jiaqi Hu
}

