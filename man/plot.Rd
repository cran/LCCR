\name{plot}
\docType{methods}
\alias{plot.confLCCR}

\title{Plot profile confidence interval}

\description{Plot for output from \code{confint.estLCCR}.}

\usage{
  \method{plot}{confLCCR}(x,...)
}

\arguments{
  \item{x}{an object of class \code{condLCCR}}
  \item{...}{further possible arguments}
}

\value{
  \item{M }{matrix of binary configurations}
}

\author{Francesco Bartolucci, Antonio Forcina}

\examples{
\donttest{
data(data_sim2)
est = estLCCR(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
conf = confint(est)
plot(conf)
}
}
