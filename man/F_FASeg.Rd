\name{F_FASeg}
\alias{F_FASeg}
\title{Joint Segmentation of Set of Correlated Time-Series}
\description{This function is dedicated to the joint segmentation (the segmentation is series-specific) in the mean of several correlated series. The form of the correlation is assumed to be arbitrary and we propose to model it with a factor model. A EM algorithm is used to estimate the parameters. A model selection procedure is also proposed to determine both the number of breakpoints and the number of factors.}
\usage{
F_FASeg(Y, uniKmax, multiKmax, qmax, selection, WithoutCorr)
}
\arguments{
  \item{Y}{Data frame, with size [(n*M) x 3], which contains the data and other informations, n is the length of each series and M is the number of series}
  \item{uniKmax}{Maximal number of segments per series (uniKmax will be lower or equal to n)}
  \item{multiKmax}{Maximal number of segments for all the series (multiKmax will be greater or equal to M)}
  \item{qmax}{Maximal number of factors (qmax will be lower or equal to M-1) (default qmax=M-1). If qmax=0 then a joint segmentation with multiKmax segments and without taking into account the correlation between series is performed}
  \item{selection}{A logical value indicating if the selection of the number of segments K and the number of factors Q is performed (default=TRUE). If it is TRUE, K and Q are selected; if it is FALSE, K is fixed to multiKmax and Q is fixed to qmax}
  \item{WithoutCorr}{A logical value indicating if, when K and Q are selected, the joint segmentation without taking into account the correlation between series is also a possible solution in the selection (default=FALSE)}
}


\value{Contains the following attributes:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{SelectedK}{Selected number of segments for all the series if selection=TRUE, the number of segments fixed by the user otherwise (K=multiKmax)}
  \item{Selectedq}{Selected number of factors if selection=TRUE, the number of factors fixed by the user otherwise (Q=qmax)}
 \item{SelectedSigma}{Estimation of the covariance matrix Sigma}
 \item{SelectedPsi}{Estimation of the matrix Psi}
 \item{SelectedB}{Estimation of the matrix of coefficients B}
 \item{SelectedZ}{Estimation of the latent vectors Z}
 \item{SelectedSeg}{Optimal segmentation with a selected or fixed value of the number of segments and the number of factors}
}
\references{A factor model approach for the joint segmentation with between-series correlation (arXiv:1505.05660)}


\author{Xavier Collilieux, Emilie Lebarbier and Stephane Robin}

%% ~Make other sections like Warning with \section{Warning }{....} ~
