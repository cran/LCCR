\name{estLCCRcon}

\alias{estLCCRcon}

\title{CML estimation of latent class models for capture-recapture data with individual covariates}

\description{For a latent class model for stratified capture-recapture data with individual covariates, it estimates the model on the basis of observed data by the conditional likelihood method, exploiting weights associated to the different strata. Estimation of the model parameters, included the population size, is based on an EM algorithm.
}

\usage{
estLCCRcon(Y, H, model = c("loglin", "logit"), W = NULL, X = NULL, N = NULL, biv = NULL,
        flag = c("no", "prev", "sum", "atleast"),
        main = c("LC", "same", "Rasch"),
        free_cov = c("no", "class", "resp", "both"),
        free_biv = c("no", "class", "int", "both"),
        free_flag = c("no", "class", "resp", "both"),
        beta0 = NULL, lambda0 = NULL, control = list(),
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
  \item{beta0}{initial value of the parameters affecting the class weights}
  \item{lambda0}{initial value of the parameters affecting the conditional distribution of capture configurations given the latent class}
  \item{control}{a list containing control arguments for the maximum number of iterations of the EM algorithm (\code{maxit}, default value 5000) and relative tollerance (\code{reltol}, default value \eqn{10^{-10}})}
  \item{verb}{to have partial output during the model fitting}
  \item{init_rand}{to use a random initialization of the parameters}
  \item{se_out}{to require computation of the standard errors}
}

\value{
  \item{beta}{estimate of the parameters affecting the class weights}
  \item{lambda}{estimate of the parameters affecting the conditional distribution of capture configurations given the latent class}
  \item{lk}{final log-likelihood value}
  \item{dev}{final value of the deviance that is used to build confidence intervals}
  \item{N}{estimate of the population size}
  \item{np}{number of free parameters}
  \item{AIC}{value of AIC for model selection}
  \item{BIC}{value of BIC for model selection}
  \item{M}{design matrices used for the recursive logit or loglinear parametrization of the conditional distribution of capture configurations given the latent class}
  \item{phiv}{estimate of the probability of being never captured for each stratum}
  \item{Piv}{matrix of the probabilities of the latent classes for each stratum}
  \item{Q}{array of the conditional probabilities of the capture configurations given each latent class and stratum}
  \item{seN}{standard error for the estimate of N}
  \item{sebeta}{standard error for the estimate of beta}
  \item{selambda}{standard error for the estimate of lambda}
}

\references{
Bartolucci, F. and Forcina, A. (2006). Estimating the size of a closed population by modeling latent and observed heterogeneity. \emph{Journal of the American Statistical Association}, \bold{101}, 786-794.
}

\author{Francesco Bartolucci, Antonio Forcina}

\seealso{
\code{\link{design_matrix_logit}}, \code{\link{design_matrix_loglin}}, \code{\link{simLCCR}}
}

\examples{
\donttest{
# estimate latent class model with 2 classes having the same weight on 5 lists
data(data_sim1)
est = estLCCRcon(Y=data_sim1,H=2)
est

# estimate latent class model with 2 classes, one covariate affecting the weights and bivariate 
# loglinear interactions between consecutive lists
data(data_sim2)
est = estLCCRcon(Y=data_sim2$Y,H=2,W=data_sim2$W,biv=matrix(c(1,2,3,4,2,3,4,5),4),main="same")
est

# estimate latent class model with 3 classes, one covariate affecting the logits of each response,
# and lag dependence
data(data_sim3)
est = estLCCRcon(Y=data_sim3$Y,H=3,model="logit",X=data_sim3$X,flag="atleast")
est

# estimate latent class model with 2 classes and covariates affecting both the class weights and
# conditional capture probabilities given the latent class
data(data_sim4)
est = estLCCRcon(Y=data_sim4$Y,H=2,X=data_sim4$X,W=data_sim4$W)
est
}
}
