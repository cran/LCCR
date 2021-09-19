\name{aggr_data}

\alias{aggr_data}

\title{Aggregate capture-recapture data}

\description{Given a matrix of individual capture configurations, and the corresponding matrix and array of covariates (affecting the class weights and the conditional capture probabilities given the latent class), the function aggregates data, providing a stratified capture-recapture format.}

\usage{
aggr_data(Y, W = NULL, X = NULL)
}

\arguments{
  \item{Y}{matrix of frequencies of capture configurations}
  \item{W}{matrix of covariates affecting the class weights}
  \item{X}{array of covariates affecting the conditional capture probabilities}
}

\value{
  \item{Ya}{matrix of aggregated frequencies of capture configurations}
  \item{Wa}{matrix of aggregated covariates affecting the class weights}
  \item{Xa}{array of aggregated covariates affecting the conditional capture probabilities}
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{freq_data}}
}

\examples{
# case without covariates
Y = rbind(c(0,1,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0),c(0,0,0,0,0,0,1,0))
aggr_data(Y)

# case with covariates
W = c(1,1,2); X = as.matrix(c(0,0,2))
aggr_data(Y,W,X)
}
