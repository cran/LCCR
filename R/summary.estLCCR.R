summary.estLCCR <-
function(object, ...){
  
  cat("\nEstimation of latent class models for capture recapture data\n")
  cat("Call:\n")
  print(object$call)
  cat("\nAvailable objects:\n")
  print(names(object))
  cat("\n")
  print(c(LogLik=object$lk,np=object$np,AIC=object$AIC,BIC=object$BIC))
  cat("\nEstimate of population size:\n")
  print(c(Nh=object$N))
  cat("\nEstimate of the parameters affecting the class weights:\n")
  print(c(object$be))
  cat("\nEstimate of the parameters affecting the conditional capture probabilities given the latent class:\n")
  print(c(object$la))
  
}
