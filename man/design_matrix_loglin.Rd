\name{design_matrix_loglin}

\alias{design_matrix_loglin}

\title{Design matrix for loglinear parametrization}

\description{For a latent class model for capture-recapture data, it provides design matrices for the loglinear parametrization used to formulate the conditional distribution of the capture configurations given the latent class. It is possible to include the effect of covariates and bivariate interactions.
}

\usage{
design_matrix_loglin(J, H = 1, main = c("LC", "same", "Rasch"), X = NULL,
                     free_cov = c("no", "class", "resp", "both"),
                     biv = NULL, free_biv = c("no", "class", "int", "both"))
}

\arguments{
  \item{J}{number of capture occasions}
  \item{H}{number of latent classes}
  \item{main}{"LC" for the latent class model in which there is a separate main effect for each capture occasion and latent class; "same" for the constrained model in which the main effect is the same across capture occasions but different across latent classes; "Rasch" for the constrained model in which each main effect is expressed in an additive form with a parameter related to the latent class and another parameter related to the capture occasion}
  \item{X}{array of covariates (n. strata x n. covariates x n. responses)}
  \item{free_cov}{"no" for constant effect of the covariates across capture occasions and latent classes; "class" for effect of covariates varying only with the latent class; "resp" for effect of covariates varying only with the capture occasion; "both" for effect of covariates varying with the capture occasion and the latent class}
  \item{biv}{matrix with two columns containing the list of bivariate interactions}
  \item{free_biv}{"no" for constant bivariate interation effect with respect to the occasion and the latent class; "class" for free interaction with respect to the latent class; "int" for free effect only with respect to the interation; "both" for free effect with respect to interation and latent class}
}

\value{
  \item{M }{design matrices}
  \item{par_list }{list of parameter names}
  \item{Main }{list of capture configurations}
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\author{Francesco Bartolucci, Antonio Foricna}

\seealso{
\code{\link{design_matrix_logit}}, \code{\link{estLCCR}}
}