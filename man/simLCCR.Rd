\name{simLCCR}

\alias{simLCCR}

\title{Simulate capture-recapture data from a latent class model with individual covariates}

\description{
The function simulates capture-recapture data from a latent class model with individual covariates that may affect the class weights and/or the conditional distribution of capture configurations given the latent class. The data may be in disaggregated form (with each stratum having unitary dimension) or aggregated form (with strata having generic dimension).
}

\usage{
simLCCR(H, J, be, la, N, model = c("loglin", "logit"), Wc = NULL, Xc = NULL, biv = NULL,
        flag = c("no", "prev", "sum", "atleast"),
        main = c("LC", "same", "Rasch"), 
        free_cov = c("no", "class", "resp", "both"),
        free_biv = c("no", "class", "int", "both"),
        free_flag = c("no", "class", "resp", "both"))
}

\arguments{
  \item{H}{number of latent classes}
  \item{J}{number of capture occasions}
  \item{be}{value of the parameters affecting the class weights}
  \item{la}{value of the parameters affecting the conditional distribution of capture configurations given the latent}
  \item{N}{population size (with individual data); vector containing the size of any stratum (with aggregated data)}
  \item{model}{population size (with individual data); vector containing the size of any stratum (with aggregated data)}
  \item{Wc}{matrix of covariates affecting the class weights at population level}
  \item{Xc}{array of covariates at population level (n. strata x n. covariates x n. responses)}
  \item{biv}{matrix with two columns containing the list of bivariate interactions (for loglinear parametrization)}
  \item{flag}{"no" for no lag effect; "prev" for effect of the previous capture occasion only; "sum" for the effect of the sum of the previous capture occasions; "atleast" for the effect of at least one capture at the previous occasions (for recursive logit parametrization)}
  \item{main}{"LC" for the latent class model in which there is a separate main effect for each capture occasion and latent class; "same" for the constrained model in which the main effect is the same across capture occasions but different across latent classes; "Rasch" for the constrained model in which each main effect is expressed in an additive form with a parameter related to the latent class and another parameter related to the capture occasion}
  \item{free_cov}{"no" for constant effect of the covariates across capture occasions and latent classes; "class" for effect of covariates varying only with the latent class; "resp" for effect of covariates varying only with the capture occasion; "both" for effect of covariates varying with the capture occasion and the latent class}
  \item{free_biv}{"no" for constant bivariate interation effect with respect to the occasion and the latent class; "class" for free interaction with respect to the latent class; "int" for free effect only with respect to the interation; "both" for free effect with respect to interation and latent class (for loglinear parametrization)}
  \item{free_flag}{"no" for constant effect of the previous capture occasions with respect to the occasion and the latent class; "class" for free effect only with respect to the latent class; "int" for free effect only with respect to the occasion; "both" for free effect with respect to capture occasion and latent class (for recursive logit parametrization)}
}
\value{
  \item{Y}{matrix of frequencies for each stratum (by row), only for captured indivdiuals}
  \item{Yc}{matrix of frequencies for each stratum (by row), for all indivdiuals}
  \item{Piv}{matrix of class weights for each stratum}
  \item{Q}{matrix of the conditional distribution of capture configurations given the latent class}
  \item{Pm}{matrix of the distribution of the capture configurations (with aggregated data)}
  \item{R}{matrix of single capture occasions (with individual data), only for captured individuals}
  \item{U}{vector of latent classes (with individual data)}
  \item{Rc}{matrix of single capture occasions (with individual data), for all individuals}
  \item{W}{matrix of covariates affecting the class weights, only for captured individuals}
  \item{X}{array of covariates affecting the conditional distribution of capture configurations given the latent class, only for captured individuals}
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{design_matrix_logit}}, \code{\link{design_matrix_loglin}}, \code{\link{estLCCR}}
}
\examples{
# simulate data from latent class model with 2 classes having the same weight on 5 lists
out = simLCCR(2,5,be=0,la=c(rep(-1,5),rep(1,5)),N=200)

# simulate data from a latent class model with 2 classes, one covariate affecting the weights and 
# bivariate loglinear interactions between consecutive lists
Wc = rnorm(200)
out = simLCCR(2,6,be=c(0,1),la=c(-1,1,1),N=200,Wc=Wc,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")

# simulate data from a latent class model with 3 classes, one covariate affecting the logits of 
# each response, and lag dependence
Xc = rnorm(200)
out = simLCCR(3,6,model="logit",be=c(0,0),la=c(rep(-1,6),rep(0,6),rep(1,6),1,1),N=200,Xc=Xc,
              flag="atleast")

# simulate data from latent class model with 2 classes and covariates affecting both class weights
# and conditional probabilities of capture configurations given the latent class
Wc = c(-1,0,1); Xc = rnorm(3)
out = simLCCR(2,5,be=c(0,0),la=c(rep(-1,5),rep(1,5),1),N=c(100,100,100),Wc=Wc,Xc=Xc)
}