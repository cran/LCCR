\name{freq_data}
\alias{freq_data}

\title{Convert capture-recapture configurations into frequency data}

\description{Given a matrix of observed capture-recapture configurations, convert the data to a matrix of frequencies having the same number of rows and a number of columns equal to the number of possible configurations.}

\usage{
freq_data(R,count=rep(1,nrow(R)))
}


\arguments{
  \item{R}{Matrix of size n (sample size) x J (number of capture occastions) of observed capture-recapture configurations}
  \item{count}{Vector of counts for each configuration}
}


\value{
  \item{Y }{Matrix of frequency of each capture configuration of size n x (2^J)}
}

\author{
Francesco Bartolucci, Antonio Forcina
}

\seealso{
\code{\link{aggr_data}}
}

\examples{
R = rbind(c(0,0,0,1),c(0,1,0,0),c(0,0,1,1),c(1,1,1,1),c(1,1,1,1))
Y = freq_data(R)
}
