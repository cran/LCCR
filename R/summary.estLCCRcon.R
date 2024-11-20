summary.estLCCRcon = function(object, ...){

  #---- summarize output of estLCCRcon ---
  cat("\nEstimation of latent class models for capture-recapture data with CML method\n")
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
      tt = object$beta/object$sebeta
      pv = 2*(1-pnorm(abs(tt)))
      Tab = cbind("est."=object$beta,"s.e."=object$sebeta,"t-test"=tt,"p-value"=pv)
      print(Tab)
    }else{
      print(object$beta)
    }
  }
  cat("\nParameters affecting the conditional capture probabilities given the latent class:\n")
  if(object$se_out){
    tt = object$lambda/object$selambda
    pv = 2*(1-pnorm(abs(tt)))
    Tab = cbind("est."=object$lambda,"s.e."=object$selambda,"t-test"=tt,"p-value"=pv)
    print(Tab)
  }else{
    print(object$lambda)
  }

}