\name{Y}
\alias{Y}
\docType{data}
\title{Matrix of data}
\description{A data frame [(n x M) x 3] containing 5 Gaussian series with size n=50 each simulated as in the paper arXiv:1505.05660 (with rho=0.6 and sigma=0.2). The total number of segments is K=11 or 6 breakpoints (at position 39 for series 1; 35 for series 2; no breaks for series 3; 11 for series 4 and 2, 3 and 12 for series 5).}
\usage{data(Y)}
\format{
  A data frame with 250 observations on the following 3 variables.
  \describe{
    \item{\code{series}}{a numeric vector}
    \item{\code{position}}{a numeric vector}
    \item{\code{signal}}{a numeric vector}
  }
}
\details{series: the number of the series; position: the grid \{1:n\}; signal: the values of the observed signal
}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(FASeg)
data(Y)
}
\keyword{datasets}
