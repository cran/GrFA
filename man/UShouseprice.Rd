\name{UShouseprice}
\alias{UShouseprice}
\docType{data}
\title{
Housing price data for 16 states in the U.S.
}
\description{
Housing price data for 16 states in the U.S over the period Jan 2000 to April 2023.
}
\usage{data("UShouseprice")}
\format{
A list with a length of 16. Each element is a matrix of dimension \eqn{T*N_m}.
}
\source{
The original data is downloaded from the website of Zillow at \url{www.zillow.com/research/data}.
}

\examples{
data(UShouseprice)
log_diff = function(x){
  T = nrow(x)
  res = log(x[2:T,]/x[1:(T-1),])*100
  scale(res, center = TRUE, scale = TRUE)
}
UShouseprice1 = lapply(UShouseprice, log_diff)
}
\keyword{datasets}
