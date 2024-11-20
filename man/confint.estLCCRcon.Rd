\name{confint.estLCCRcon}

\alias{confint.estLCCRcon}

\title{Confidence interval}

\description{Profile confidence interval for a latent class model with covariates based on the CML method.}

\usage{
\method{confint}{estLCCRcon}(object, parm = list(), level = 0.95, ...)
}

\arguments{
  \item{object}{output from \code{\link{estLCCR}}}
  \item{parm}{a list containing control arguments for the step length of the N values(\code{step}, default value depending on the estimate of N), range of N values in terms of distance of the log-likelhood from the its maximum (\code{mult}, default value 1.5), and maximum value of this grid as a multiple of the estimate of this parameter (\code{max}, default value 5)}
  \item{level}{the confidence level required}
  \item{...}{additional argument(s) for methods}
}

\value{
  \item{conf }{confidence interval}
  \item{Nv }{vector of values of N tried}
  \item{lkv }{corresponding vector of log-likelihood values}
  \item{level }{level of confidence}
  \item{Nh }{point estimate of N}
  \item{devh }{minimum deviance}
  \item{dev1 }{deviance at the lower bound}
  \item{dev2 }{deviance at the upper bound}
  \item{step }{step used in the grid search}
  \item{err }{error in the grid search}
}

\references{
Bartolucci, F. and Forcina, A. (2006). Estimating the size of a closed population by modeling latent and observed heterogeneity. \emph{Journal of the American Statistical Association}, \bold{101}, 786-794.
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{estLCCR}}, \code{\link{plot.confLCCR}}
}

\examples{
\donttest{
data(data_sim2)
est = estLCCRcon(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
conf = confint(est)
plot(conf)
}
}
