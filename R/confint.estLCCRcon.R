confint.estLCCRcon = function(object,parm=list(),level=0.95,...){

#---- preliminaries ----
  q = qchisq(level,1)

#---- build confidence interval ----
  devh = object$dev
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
  devv = devh
  Nv1 = Nh
  devv1 = devh
  N = round(Nh/parm$step)*parm$step
  dev = devh
  n = sum(object$Y)
  cat(paste("step =", parm$step))
  cat("\n")
  cat("------------|-------------|-------------|\n")
  cat("      N     |     dev     | dev-min(dev)|\n")
  cat("------------|-------------|-------------|\n")
  it = 0
# iterate for values of N smaller than the point estimate
  while(dev-devh<q*parm$mult & N>n){
    it = it+1
    N = N-parm$step
    out = estLCCRcon(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,beta0=beta,lambda0=lambda,
                  control=list(maxit=500),verb=FALSE)
    dev = out$dev
    be = out$be
    la = out$la
    Nv = c(N,Nv)
    devv = c(dev,devv)
    Nv1 = c(N,Nv1)
    devv1 = c(dev,devv1)
    if(it%%10==0) cat(sprintf("%11g", c(N,dev,dev-devh)), "\n", sep = " | ")
  }
  if(it%%10>0) cat(sprintf("%11g", c(N,dev,dev-devh)), "\n", sep = " | ")
  cat("------------|-------------|-------------|\n")
  N = round(Nh/parm$step)*parm$step
  dev = devh
  Nv2 = Nh
  devv2 = devh
  be = object$be
  la = object$la
  it = 0
# iterate for values of N greater than the point estimate
  while(dev-devh<q*parm$mult & N<parm$max*Nh){
    it = it+1
    N = N+parm$step
    out = estLCCRcon(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,beta0=beta,lambda0=lambda,control=list(maxit=500),
                  verb=FALSE)
    dev = out$dev
    be = out$be
    la = out$la
    Nv = c(Nv,N)
    devv = c(devv,dev)
    Nv2 = c(Nv2,N)
    devv2 = c(devv2,dev)
    if(it%%10==0) cat(sprintf("%11g", c(N,dev,dev-devh)), "\n", sep = " | ")
  }
  if(it%%10>0) cat(sprintf("%11g", c(N,dev,dev-devh)), "\n", sep = " | ")
  cat("------------|-------------|-------------|\n")

#---- find limits ----
  ta = devh+q
  err = FALSE
  ind1 = NULL
  if(length(Nv1)>1) for(j in 1:(length(Nv1)-1)) if(devv1[j]>=ta & ta>devv1[j+1]) ind1 = j
  if(is.null(ind1)){
    ind1 = which.min(abs(devv1-devh-q))
    N1 = Nv1[ind1]
    dev1 = devv1[ind1]
    err = TRUE
  }else{
    N1 = Nv1[ind1]+(ta-devv1[ind1])/(devv1[ind1+1]-devv1[ind1])*(Nv1[ind1+1]-Nv1[ind1])
    dev1 = ta
  }
  ind2 = NULL
  if(length(Nv2)>1) for(j in 1:(length(Nv2)-1)) if(devv2[j]<ta & ta<=devv2[j+1]) ind2 = j
  if(is.null(ind2)){
    ind2 = which.min(abs(devv2-devh-q))
    N2 = Nv2[ind2]
    dev2 = devv2[ind2]
    err = TRUE
  }else{
    N2 = Nv2[ind2+1]+(ta-devv2[ind2+1])/(devv2[ind2]-devv2[ind2+1])*(Nv2[ind2]-Nv2[ind2+1])
    dev2 = ta
  }
  conf = c(N1,N2)

#---- output ----
  out = list(conf=conf,Nv=Nv,devv=devv,level=level,Nh=Nh,devh=devh,dev1=dev1,dev2=dev2,
             step=parm$step,err=err,call = match.call())
  class(out) = "confLCCRcon"
  return(out)

}
