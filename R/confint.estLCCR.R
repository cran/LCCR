confint.estLCCR = function(object,parm=list(),level=0.95,...){

#---- preliminaries ----
  q = qchisq(level,1)

#---- build confidence interval ----
  lkh = object$lk
  Nh = object$N
  beta = object$beta
  lambda = object$lambda
  if(is.null(parm$step)){
    if(Nh<=50) parm$step = 0.1
    if(Nh>50 & Nh<=100) parm$step = 0.2
    if(Nh>100 & Nh<=200) parm$step = 0.5
    if(Nh>200) parm$step = 1
  }
  if(is.null(parm$mult)) parm$mult = 1.5
  if(is.null(parm$max)) parm$max = 5
  Nv = Nh
  lkv = lkh
  Nv1 = Nh
  lkv1 = lkh
  N = round(Nh/parm$step)*parm$step
  lk = lkh
  n = sum(object$Y)
  cat(paste("step =", parm$step))
  cat("\n")
  cat("------------|-------------|-------------|\n")
  cat("      N     |      lk     |  lk-max(lk) |\n")
  cat("------------|-------------|-------------|\n")
  it = 0
# iterate for values of N smaller than the point estimate
  while(2*(lkh-lk)<q*parm$mult & N>n){
    it = it+1
    N = N-parm$step
    out = estLCCR(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,beta0=beta,lambda0=lambda,
                  control=list(maxit=500),verb=FALSE)
    lk = out$lk
    be = out$be
    la = out$la
    Nv = c(N,Nv)
    lkv = c(lk,lkv)
    Nv1 = c(N,Nv1)
    lkv1 = c(lk,lkv1)
    if(it%%10==0) cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  }
  if(it%%10>0) cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  cat("------------|-------------|-------------|\n")
  N = round(Nh/parm$step)*parm$step
  lk = lkh
  Nv2 = Nh
  lkv2 = lkh
  be = object$be
  la = object$la
  it = 0
# iterate for values of N greater than the point estimate
  while(2*(lkh-lk)<q*parm$mult & N<parm$max*Nh){
    it = it+1
    N = N+parm$step
    out = estLCCR(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,beta0=beta,lambda0=lambda,
                  control=list(maxit=500),verb=FALSE)
    lk = out$lk
    be = out$be
    la = out$la
    Nv = c(Nv,N)
    lkv = c(lkv,lk)
    Nv2 = c(Nv2,N)
    lkv2 = c(lkv2,lk)
    if(it%%10==0) cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  }
  if(it%%10>0) cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  cat("------------|-------------|-------------|\n")

#---- find limits ----
  ta = lkh-q/2
  err = FALSE
  ind1 = NULL
  if(length(Nv1)>1) for(j in 1:(length(Nv1)-1)) if(lkv1[j]<=ta & ta<lkv1[j+1]) ind1 = j
  if(is.null(ind1)){
    ind1 = which.min(abs(lkh-lkv1-q/2))
    N1 = Nv1[ind1]
    lk1 = lkv1[ind1]
    err = TRUE
  }else{
    N1 = Nv1[ind1]+(ta-lkv1[ind1])/(lkv1[ind1+1]-lkv1[ind1])*(Nv1[ind1+1]-Nv1[ind1])
    lk1 = ta
  }
  ind2= NULL
  if(length(Nv2)>1) for(j in 1:(length(Nv2)-1)) if(lkv2[j]>ta & ta>=lkv2[j+1]) ind2 = j
  if(is.null(ind2)){
    ind2 = which.min(abs(lkh-lkv2-q/2))
    N2 = Nv2[ind2]
    lk2 = lkv2[ind2]
    err = TRUE
  }else{
    N2 = Nv2[ind2+1]+(ta-lkv2[ind2+1])/(lkv2[ind2]-lkv2[ind2+1])*(Nv2[ind2]-Nv2[ind2+1])
    lk2 = ta
  }
  conf = c(N1,N2)

#---- output ----
  out = list(conf=conf,Nv=Nv,lkv=lkv,level=level,Nh=Nh,lkh=lkh,lk1=lk1,lk2=lk2,
             step=parm$step,err=err,call = match.call())
  class(out) = "confLCCR"
  return(out)

}
