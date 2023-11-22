\name{gendata}
\alias{gendata}
\title{
Generate the grouped data.
}
\description{
Generate the grouped data.
}
\usage{
gendata(seed = 1, T = 50, N = rep(30, 5), r0 = 2, r = rep(2, 5),
        Phi_G = 0.5, Phi_F = 0.5, Phi_e = 0.5, W_F = 0.5, beta = 0.1,
        kappa = 1, case = 1)
}

\arguments{
  \item{seed}{the seed used in \code{set.seed}.}
  \item{T}{the number of time points.}
  \item{N}{a vector representing the number of variables in each group.}
  \item{r0}{the number of global factors.}
  \item{r}{a vector representing the number of the local factors. Notice, the length of \eqn{r} is the same as \eqn{N}.}
  \item{Phi_G}{hyperparameter of the global factors, default is 0.5, the value should between 0 and 1.}
  \item{Phi_F}{hyperparameter of the local factors, default is 0.5, the value should between 0 and 1.}
  \item{Phi_e}{hyperparameter of the errors, default is 0.5, the value should between 0 and 1.}
  \item{W_F}{hyperparameter of the correlation of local factors, only applicable in \code{case = 3}, the value should between 0 and 1.}
  \item{beta}{hyperparameter of the errors, default is 0.1.}
  \item{kappa}{hyperparameter of signal to noise ratio, default is 1.}
  \item{case}{the case of the data-generating process, default is 1, it can also be 2 and 3.}
}

\value{
  \item{y}{a list of the data.}
  \item{G}{the global factors.}
  \item{F}{a list of the local factors.}
  \item{loading_G}{the global factor loadings.}
  \item{loading_F}{the local factor loadings.}
  \item{T}{the number of time points.}
  \item{N}{a vector representing the number of variables in each group.}
  \item{M}{the number of groups.}
  \item{r0}{the number of global factors.}
  \item{r}{a vector representing the number of the local factors.}
  \item{case}{the case of the data-generating process.}
}
\examples{
dat = gendata()
dat
}
