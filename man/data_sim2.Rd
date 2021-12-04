\name{data_sim2}

\alias{data_sim2}

\docType{data}

\title{Simulated data 2}

\description{Simulated capture-recapture data from a latent class model with 2 classes, one covariate affecting the class weights, and bivariate loglinear interactions between consecutive lists.}

\usage{data(data_sim2)}

\format{
A list containing capture configurations at individual level and corresponding covariates.
  \describe{
    \item{\code{Y}}{capture configurations}
    \item{\code{W}}{values of the covariate affecting the class weights}
  }
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\examples{
\donttest{
data(data_sim2)
est = estLCCR(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
}
}

\keyword{datasets}
