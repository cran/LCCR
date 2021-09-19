\name{estLCCR}

\alias{estLCCR}

\title{Estimate latent class models for capture-recapture data with individual covariates}

\description{For a latent class model for stratified capture-recapture data with individual covariates, it estimates the model on the basis of observed data by the unconditional likelihood method, exploiting weights associated to the different strata. Estimation of the model parameters, included the population size, is based on an EM algorithm.
}

\usage{
estLCCR(Y, H, model = c("loglin", "logit"), W = NULL, X = NULL, N = NULL, biv = NULL,
        flag = c("no", "prev", "sum", "atleast"),
        main = c("LC", "same", "Rasch"),
        free_cov = c("no", "class", "resp", "both"),
        free_biv = c("no", "class", "int", "both"),
        free_flag = c("no", "class", "resp", "both"),
        be = NULL, la = NULL, maxit = 5000,
        verb = TRUE, init_rand = FALSE, se_out = FALSE)
}

\arguments{
  \item{Y}{matrix of frequencies for each stratum (by row)}
  \item{H}{number of latent classes}
  \item{model}{"loglin" for loglinear parametrization; "logit" for recursive logit parametrization}
  \item{W}{matrix of covariates on the class weights}
  \item{X}{array of covariates (n. strata x n. covariates x n. responses)}
  \item{N}{fixed population size}
  \item{biv}{matrix with two columns containing the list of bivariate interactions (for loglinear parametrization)}
  \item{flag}{"no" for no lag effect; "prev" for effect of the previous capture occasion only; "sum" for the effect of the sum of the previous capture occasions; "atleast" for the effect of at least one capture at the previous occasions (for recursive logit parametrization)}
  \item{main}{"LC" for the latent class model in which there is a separate main effect for each capture occasion and latent class; "same" for the constrained model in which the main effect is the same across capture occasions but different across latent classes; "Rasch" for the constrained model in which each main effect is expressed in an additive form with a parameter related to the latent class and another parameter related to the capture occasion}
  \item{free_cov}{"no" for constant effect of the covariates across capture occasions and latent classes; "class" for effect of covariates varying only with the latent class; "resp" for effect of covariates varying only with the capture occasion; "both" for effect of covariates varying with the capture occasion and the latent class}
  \item{free_biv}{"no" for constant bivariate interation effect with respect to the occasion and the latent class; "class" for free interaction with respect to the latent class; "int" for free effect only with respect to the interation; "both" for free effect with respect to interation and latent class (for loglinear parametrization)}
  \item{free_flag}{"no" for constant effect of the previous capture occasions with respect to the occasion and the latent class; "class" for free effect only with respect to the latent class; "int" for free effect only with respect to the occasion; "both" for free effect with respect to capture occasion and latent class (for recursive logit parametrization)}
  \item{be}{initial value of the parameters affecting the class weights}
  \item{la}{initial value of the parameters affecting the conditional distribution of capture configurations given the latent class}
  \item{maxit}{maximum number of iterations of the EM algorithm}
  \item{verb}{to have partial output during the model fitting}
  \item{init_rand}{to use a random initialization of the parameters}
  \item{se_out}{to require computation of the standard errors}
}

\value{
  \item{be}{estimate of the parameters affecting the class weights}
  \item{la}{estimate of the parameters affecting the conditional distribution of capture configurations given the latent class}
  \item{lk}{final log-likelihood value}
  \item{N}{estimate of the population size}
  \item{np}{number of free parameters}
  \item{AIC}{value of AIC for model selection}
  \item{BIC}{value of BIC for model selection}
  \item{M}{design matrices used for the recursive logit or loglinear parametrization of the conditional distribution of capture configurations given the latent class}
  \item{tauv}{estimate of the weights of each stratum}
  \item{phiv}{estimate of the probability of being never captured for each stratum}
  \item{seN}{standard error for the estimate of N}
  \item{sebe}{standard error for the estimate of beta}
  \item{sela}{standard error for the estimate of lambda}
  \item{lk1}{component of the log-likelihood based on the binomial factor in N}
  \item{lk2}{component of the log-likelihood involving N and the overall probability of being never captured}
  \item{lk3}{component of the log-likelihood involving the capture probabilities}
  \item{lk4}{component of the log-likelihood involving the stratum weights}
}

\references{
Forcina, A. and Bartolucci, F. (2021). Estimating the size of a closed population by modeling latent and observed heterogeneity, \emph{arXiv:2106.03811}.

Liu, Y., Li, P., and Qin, J. (2017). Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data. \emph{Biometrika}, \bold{104}, 527-543.
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{design_matrix_logit}}, \code{\link{design_matrix_loglin}}, \code{\link{simLCCR}}
}

\examples{
\donttest{
# estimate latent class model with 2 classes having the same weight on 5 lists
data(data_sim1)
est = estLCCR(Y=data_sim1,H=2)
est

# estimate latent class model with 2 classes, one covariate affecting the weights and bivariate 
# loglinear interactions between consecutive lists
data(data_sim2)
est = estLCCR(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
est

# estimate latent class model with 3 classes, one covariate affecting the logits of each response,
# and lag dependence
data(data_sim3)
est = estLCCR(Y=data_sim3$Y,H=3,model="logit",X=data_sim3$X,flag="atleast")
est

# estimate latent class model with 2 classes and covariates affecting both the class weights and
# conditional capture probabilities given the latent class
data(data_sim4)
est = estLCCR(Y=data_sim4$Y,H=2,X=data_sim4$X,W=data_sim4$W)
est
}
}
