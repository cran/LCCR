\name{matrix_logit}

\alias{matrix_logit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Design matrices of saturated logit parametrization}

\description{For J binary responses, it provides design matrices A and B used to express the saturated recursive logit parametrization of type log(p) = [A*eta-B*log(1+exp(eta))].}

\usage{matrix_logit(J)}

\arguments{
  \item{J}{number of capture occasions}
}

\value{
  \item{A }{first matrix involved in the recursive logit parametrization}
  \item{B }{second matrix involved in the recursive logit parametrization}
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\author{Francesco Bartolucci, Antonio Foricna}

\seealso{
\code{\link{matrix_logit}}, \code{\link{estLCCR}}
}

