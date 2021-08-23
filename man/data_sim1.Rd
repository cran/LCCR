\name{data_sim1}

\alias{data_sim1}

\docType{data}

\title{Simulated data 1}

\description{Simulated capture-recapture data from a latent class model with 2 classes having the same weight on 5 lists.}

\usage{data(data_sim1)}

\format{
  A matrix containing the capture configurations at individual level.
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\examples{
\donttest{
data(data_sim1)
est = estLCCR(Y=data_sim1,H=2)
}
}

\keyword{datasets}
