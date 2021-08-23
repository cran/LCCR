\name{design_matrix_logit}

\alias{design_matrix_logit}

\title{Design matrix for recursive logit parametrization}

\description{For a latent class model for capture-recapture data, it provides design matrices for the recursive logit parametrization used to formulate the conditional distribution of capture configurations given the latent class. It is possible to include the effect of covariates and of previous captures.}

\usage{
design_matrix_logit(J, H = 1, main = c("LC", "same", "Rasch"), X = NULL,
                    free_cov = c("no", "class", "resp", "both"),
                    flag = c("no", "prev", "sum", "atleast"),
                    free_flag = c("no", "class", "resp", "both"))
}

\arguments{
  \item{J}{number of capture occasions}
  \item{H}{number of latent classes}
  \item{main}{"LC" for the latent class model in which there is a separate main effect for each capture occasion and latent class; "same" for the constrained model in which the main effect is the same across capture occasions but different across latent classes; "Rasch" for the constrained model in which each main effect is expressed in an additive form with a parameter related to the latent class and another parameter related to the capture occasion}
  \item{X}{array of covariates (n. strata x n. covariates x n. responses)}
  \item{free_cov}{"no" for constant effect of the covariates across capture occasions and latent classes; "class" for effect of covariates varying only with the latent class; "resp" for effect of covariates varying only with the capture occasion; "both" for effect of covariates varying with the capture occasion and the latent class}
  \item{flag}{"no" for no lag effect; "prev" for effect of the previous capture occasion only; "sum" for the effect of the sum of the previous capture occasions; "atleast" for the effect of at least one capture at the previous occasions}
  \item{free_flag}{"no" for constant effect of the previous capture occasions with respect to the occasion and the latent class; "class" for free effect only with respect to the latent class; "int" for free effect only with respect to the occasion; "both" for free effect with respect to capture occasion and latent class}
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
\code{\link{design_matrix_loglin}}, \code{\link{matrix_logit}}, \code{\link{estLCCR}}
}
