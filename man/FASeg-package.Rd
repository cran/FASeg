\name{FASeg-package}
\alias{FASeg-package}
\alias{FASeg}
\docType{package}
\title{Joint Segmentation of Set of Correlated Time-Series}
\description{FASeg contains a function designed to the joint segmentation (the segmentation is series-specific) in the mean of several correlated series. The form of the correlation is assumed to be arbitrary and we propose to model it with a factor model. A EM algorithm is used to estimate the parameters and a model selection strategy is proposed to determine both the number of breakpoints and the number of factors}

\author{
\packageAuthor{FASeg}

Maintainer: \packageMaintainer{FASeg}
}
%\author{Xavier Collilieux, Emilie Lebarbier and Stephane Robin; 
%Maintainer: E. Lebarbier <emilie.lebarbier@agroparistech.fr>
%}
\references{A factor model approach for the joint segmentation with between-series correlation (arXiv:1505.05660)}

\keyword{Dynamic programming; EM algorithm; Factor model; Segmentation; Model selection; Multivariate time-series}

\examples{
library(FASeg)
data(Y)
M=max(Y$series)
uniKmax=3
multiKmax=11
qmax=M-1
selection=FALSE
WithoutCorr=FALSE
seg=F_FASeg(Y,uniKmax,multiKmax,qmax,selection,WithoutCorr)
}
