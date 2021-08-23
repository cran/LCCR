confint.estLCCR <-
function(object,parm=0.1,level=0.95,...){

# ---- preliminaries ----
  q = qchisq(level,1)
  step = parm

# ---- build confidence interval ----
  lkh = object$lk
  Nh = object$N
  be = object$be; la = object$la
  Nv = Nh; lkv = lkh
  Nv1 = lkv1 = NULL
  N = Nh; lk = lkh
  n = sum(object$Y)
  cat("------------|-------------|-------------|\n")
  cat("      N     |      lk     |  lk-max(lk) |\n")
  cat("------------|-------------|-------------|\n")
  while(2*(lkh-lk)<q*1.5 & N-step>n){
    N = N-step
    out = estLCCR(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,be,la,500,FALSE)
    lk = out$lk
    be = out$be; la = out$la
    Nv = c(N,Nv); lkv = c(lk,lkv)
    Nv1 = c(N,Nv1); lkv1 = c(lk,lkv1)
    cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  }
  cat("------------|-------------|-------------|\n")
  N = Nh; lk = lkh
  Nv2 = lkv2 = NULL
  be = object$be; la = object$la
  while(2*(lkh-lk)<q*1.5){
    N = N+step
    out = estLCCR(object$Y,object$H,object$model,object$W,object$X,
                  N,object$biv,object$flag,object$main,object$free_cov,
                  object$free_biv,object$free_flag,be,la,500,FALSE)
    lk = out$lk
    be = out$be; la = out$la
    Nv = c(Nv,N); lkv = c(lkv,lk)
    Nv2 = c(Nv2,N); lkv2 = c(lkv2,lk)
    cat(sprintf("%11g", c(N,lk,lk-lkh)), "\n", sep = " | ")
  }
  cat("------------|-------------|-------------|\n")

# ---- find limits ----
  ind1 = which.min(abs(lkh-lkv1-q/2))
  N1 = Nv1[ind1]; lk1 = lkv1[ind1]
  ind2 = which.min(abs(lkh-lkv2-q/2))
  N2 = Nv2[ind2]; lk2 = lkv2[ind2]
  conf = c(N1,N2)
  
# ---- output ----
  out = list(conf=conf,Nv=Nv,lkv=lkv,level=level,Nh=Nh,lkh=lkh,lk1=lk1,lk2=lk2,
             call = match.call())
  class(out) = "confLCCR"
  return(out)

}
