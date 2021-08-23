\name{data_sim4}

\alias{data_sim4}

\docType{data}

\title{Simulated data 4}

\description{Simulated stratified capture-recapture data from a latent class model with 2 classes and covariates affecting both the class weights and the conditional capture probabilities given the latent class.}

\usage{data(data_sim4)}

\format{
A list containing capture configurations at individual level and corresponding covariates.
  \describe{
    \item{\code{Y}}{capture configurations}
    \item{\code{W}}{values of the covariate affecting the class weights}
    \item{\code{X}}{values of the covariate affecting the conditional capture probabilities given the latent class}
  }
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\examples{
data(data_sim4)
est = estLCCR(Y=data_sim4$Y,H=2,X=data_sim4$X,W=data_sim4$W)
}

\keyword{datasets}
