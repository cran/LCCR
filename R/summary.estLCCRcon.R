summary.estLCCRcon <-
function(object, ...){
  
  cat("\nEstimation of latent class models for capture-recapture data\n")
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\nAvailable objects:\n")
  print(names(object))
  cat("\n")
  print(c(LogLik=object$lk,np=object$np,AIC=object$AIC,BIC=object$BIC))
  cat("\nPopulation size:\n")
  if(object$se_out){
    Tab = cbind("est."=object$N,"s.e."=object$seN)
    rownames(Tab) = "N"
    print(Tab)
  }else{
    print(c(Nh = object$N))
  }
  if(object$H>1){
    if(object$se_out){
      cat("\nParameters affecting the class weights:\n")
      tt = object$be/object$sebe
      pv = 2*(1-pnorm(abs(tt)))
      Tab = cbind("est."=object$be,"s.e."=object$sebe,"t-test"=tt,"p-value"=pv)
      print(Tab)
    }else{
      print(object$be)
    }
  }
  cat("\nParameters affecting the conditional capture probabilities given the latent class:\n")
  if(object$se_out){
    tt = object$la/object$sela
    pv = 2*(1-pnorm(abs(tt)))
    Tab = cbind("est."=object$la,"s.e."=object$sela,"t-test"=tt,"p-value"=pv)
    print(Tab)
  }else{
    print(object$la)
  }

}