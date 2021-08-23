\name{data_sim3}

\alias{data_sim3}

\docType{data}

\title{Simulated data 3}

\description{Simulated capture-recapture data from a latent class model with 3 classes, one covariate affecting the logits of each response, and lag dependence.}

\usage{data(data_sim3)}

\format{
A list containing capture configurations at individual level and corresponding covariates.
  \describe{
    \item{\code{Y}}{capture configurations}
    \item{\code{X}}{values of the covariate affecting the conditional capture probabilities given the latent class}
  }
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\examples{
  \donttest{
  data(data_sim3)
  est = estLCCR(Y=data_sim3$Y,H=3,model="logit",X=data_sim3$X,flag="atleast")
  }
}

\keyword{datasets}
