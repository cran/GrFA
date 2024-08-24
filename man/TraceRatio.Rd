\name{TraceRatio}
\alias{TraceRatio}
\title{
Trace ratio
}
\description{
Evaluation of the estimated factors by trace ratios, the values is between 0 and 1, higher values means better estimation.
}
\usage{
TraceRatio(G, Ghat)
}
\arguments{
  \item{G}{the true factors.}
  \item{Ghat}{the estimated factors.}
}
\value{
  \item{trace ratio}{defined as \eqn{\mathrm{TR} = \mathrm{tr} ( \mathbf{G}' \widehat{\mathbf{G}} (\widehat{\mathbf{G}}'\widehat{\mathbf{G}})^{-1} \widehat{\mathbf{G}}'\mathbf{G})/\mathrm{tr}(\mathbf{G'G})}.}
}
