\name{confint.estLCCR}

\alias{confint.estLCCR}

\title{Confidence interval}

\description{Profile confidence interval for a latent class model with covariates.}

\usage{
\method{confint}{estLCCR}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{object}{output from \code{\link{estLCCR}}}
  \item{parm}{step of the grid of values of the population size}
  \item{level}{the confidence level required}
  \item{...}{additional argument(s) for methods}
}

\value{No return value}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.

Liu, Y., Li, P., and Qin, J. (2017). Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data. \emph{Biometrika}, \bold{104}, 527-543.
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{estLCCR}}, \code{\link{plot.confLCCR}}
}

\examples{
\donttest{
data(data_sim2)
est = estLCCR(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
conf = confint(est)
plot(conf)
}
}
