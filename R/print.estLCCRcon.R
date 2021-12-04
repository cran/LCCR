print.estLCCRcon <-
function(x, ...){
  
  cat("\nEstimation of latent class models for capture-recapture data\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\nAvailable objects:\n")
  print(names(x))
  cat("\n")
  print(c(LogLik=x$lk,np=x$np,AIC=x$AIC,BIC=x$BIC))
  cat("\nEstimate of population size:\n")
  print(c(Nh=x$N))
  
}
