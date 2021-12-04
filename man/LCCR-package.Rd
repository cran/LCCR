\name{LCCR-package}

\alias{LCCR-package}

\alias{LCCR}

\docType{package}

\title{
\packageTitle{LCCR}
}

\description{
\packageDescription{LCCR}
}

\details{
The package provides functions to specify latent class models for capture-recapture data with individual covariates, simulate data from the specified model and estimate it on the basis of observed data, including the possibility to obtain a profile confidence interval for the population size.

The latent class models may be specified by a log-linear or a recursive logit parametrization of the conditional distribution of the capture configurations given the latent class. The log-linear parametrization is based on main effects that may depend on some covariates and may include bivariate interaction terms. The recursive logit parametrization allows us to model the dependence of capture probabilities on previous capture history. Linear constraints on the model parameters may be assumed to make the model more parsimonious. Covariates may also affect the latent class weights by a Multinomial logit parametrization.

Estimation of the specified model is based on the unconditional likelihood method following the approach formalized in Forcina and Bartolucci (2021), which may be seen as a refinement of the one in Liu et al. (2017). The approach uses as additional parameters the weights associated to the different strata that are estimated together with the model parameters and the population size. Functions for constructing a profile confidence interval for the population size are provided; this proceduire is recommended as an alternative to using the normal approximation.

Functions for conditional maximum likelihood estimation as described in Bartolucci and Forcina (2006) are also included in the package.
}

\author{
\packageAuthor{LCCR}

Maintainer: \packageMaintainer{LCCR}
}

\references{
Bartolucci, F. and Forcina, A. (2006). Estimating the size of a closed population by modeling latent and observed heterogeneity. \emph{Journal of the American Statistical Association}, \bold{101}, 786-794.

Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.

Liu, Y., Li, P., and Qin, J. (2017). Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data. \emph{Biometrika}, \bold{104}, 527-543.
}

\examples{
\donttest{
data(data_sim2)
# estimate latent class model with 2 classes, same main log-linear effect across lists,
# one covariate affecting the weights, and bivariate interaction between consecutive lists
est = estLCCR(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
est
# compute profile confidence interval
conf = confint(est)
conf
plot(conf)
}
}
