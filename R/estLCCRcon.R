estLCCRcon = function(Y,H,model=c("loglin","logit"),W=NULL,X=NULL,N=NULL,biv=NULL,
                      flag=c("no","prev","sum","atleast"),
                      main=c("LC","same","Rasch"),
                      free_cov=c("no","class","resp","both"),
                      free_biv=c("no","class","int","both"),
                      free_flag=c("no","class","resp","both"),
                      beta0=NULL,lambda0=NULL,control=list(),
                      verb=TRUE,init_rand=FALSE,se_out=FALSE){

#---- check and fix input arguments ----
  if(H<1) stop("H cannot be smaller than 1")
  if(!is.numeric(Y)) stop("Y must be numeric")
  if(!is.null(W)) if(!is.numeric(W)) stop("W must be numeric")
  if(!is.null(X)) if(!is.numeric(X)) stop("X must be numeric")
  if(!is.null(N)) if(!is.numeric(N)) stop("N must be numeric")
  if(!is.null(biv)) if(!is.numeric(biv)) stop("biv must be numeric")
  if(is.null(control$maxit)) control$maxit = 5000
  if(is.null(control$reltol)) control$reltol = 10^-10

#---- preliminaries ----
  model = match.arg(model)
  flag = match.arg(flag)
  main = match.arg(main)
  free_cov = match.arg(free_cov)
  free_biv = match.arg(free_biv)
  free_flag = match.arg(free_flag)
  estN = is.null(N)
  S = nrow(Y)
  J = round(log(ncol(Y))/log(2))
  Y1 = matrix(Y[,-1],S,2^J-1)
  nv = rowSums(Y)
# if there are no covariates W
  if(is.null(W)){
    ncov2 = 0
# if there are covariates W
  }else{
    if(is.vector(W)) W = as.matrix(W)
    ncov2 = ncol(W)
  }
  np2 = (H-1)*(1+ncov2)
  if(H>1){
    if(is.null(W)){
      Wnames = NULL
    }else{
      Wnames = colnames(W)
      if(is.null(Wnames)) Wnames = paste("W",1:ncov2,sep="")
    }
    beta_names = NULL
    for(h in 2:H) if(ncov2>0){
      beta_names = c(beta_names,paste("class",h,".int",sep=""),paste("class",h,".",Wnames,sep=""))
    }else{
      beta_names = c(beta_names,paste("class",h,".int",sep=""))
    }
  }
  if(model=="loglin") out = design_matrix_loglin(J,H,main,X,free_cov,biv,free_biv)
  if(model=="logit") out = design_matrix_logit(J,H,main,X,free_cov,flag,free_flag)
  L = array(0,c(S,H,(H-1)*(1+ncov2)))
  Tmp = diag(H)[,-1,drop=FALSE]
  for(s in 1:S) L[s,,] = Tmp%x%t(c(1,W[s,]))
  if(model=="logit"){
    A = out$A
    B = out$B
  }
  M = out$M
  np1 = length(out$par_list)
  np = np1+np2

#---- starting values ----
  n = sum(nv)
# starting value of lambda
  if(is.null(lambda0)){
    if(init_rand){
      lambda = rnorm(np1)
    }else{
      if(H==1){
        lambda = rep(0,np1)
      }else{
        est = estLCCR(Y,1,model=model,main=main,verb=verb,control=control)
        if(main=="LC" || main=="same"){
          lambda = NULL
          for(h in 1:H) lambda = c(lambda,est$lambda+h-H/2)
        }
        if(main=="Rasch") lambda = c(1:(H-1),est$lambda-(H-1)/2)
        lambda = c(lambda,rep(0,np1-length(lambda)))
      }
    }
  }else{
    lambda = lambda0
  }
# starting value of beta
  if(is.null(beta0)){
    if(H==1) beta = NULL else if(init_rand) beta = rnorm(np2) else beta = rep(0,np2)
  }else{
    beta = beta0
  }
  if(H==1){
    Piv = matrix(1,S,1)
  }else{
    Piv = matrix(0,S,H)
    for(s in 1:S){
      Ls = as.matrix(L[s,,])
      tmp = exp(Ls%*%beta)
      Piv[s,] = tmp/sum(tmp)
    }
  }
# compute probabilities
  Pj = Q = array(0,c(S,2^J,H))
  Pm = matrix(0,S,2^J)
  for(s in 1:S){
    for(h in 1:H){
      if(is.null(X)) Msh = as.matrix(M[,,h,1]) else Msh = as.matrix(M[,,h,s])
      if(model=="loglin"){
        tmp = exp(c(Msh%*%lambda))
        Q[s,,h] = tmp/sum(tmp)
      }
      if(model=="logit") Q[s,,h] = exp(A%*%Msh%*%lambda-B%*%log(1+exp(Msh%*%lambda)))
    }
    if(H>1){
      Pj[s,,] = Q[s,,]*(rep(1,2^J)%o%Piv[s,])
      Pm[s,] = Q[s,,]%*%Piv[s,]
    }
  }
  if(H==1){
    Pj = Q
    Pm = matrix(Q[,,1],S,2^J)
  }
  phiv = Pm[,1]
  Pm1 = Pm[,-1]
  Qm1 = (1/(1-phiv))*Pm1
# starting values of Nv
  if(estN){
    Nv = nv/(1-phiv)
    N = sum(Nv)
  }else{
    Nv = nv*N/n
  }
  if(!estN){
    psi = 0
    ta = N-sum(nv/(1-exp(psi)*phiv))
  }
# compute log-likelhood and deviance
  lk = sum(Y1*log(Qm1))
  Yh1 = pmax(cbind(Nv-nv,Y1),10^-300)
  dev = 2*sum(Yh1*log(Yh1/(Nv*Pm)))
  if(verb){
    cat("------------|-------------|-------------|-------------|\n")
    if(estN){
      cat("  iteration |      lk     |    lk-lko   |      N      |\n")
      cat("------------|-------------|-------------|-------------|\n")
      cat(sprintf("%11g", c(0,lk,NA,N)), "\n", sep = " | ")
    }else{
      cat("  iteration |      dev    |   dev-devo  |      N      |\n")
      cat("------------|-------------|-------------|-------------|\n")
      cat(sprintf("%11g", c(0,dev,NA,N)), "\n", sep = " | ")
    }
  }

#---- cycle ----
  dM1 = dim(M)[1]
  dM2 = dim(M)[2]
  it = 0
  devo = dev
  while((abs(dev-devo)/devo>control$reltol | it==0) & it<control$maxit){
    it = it+1
    lko = lk
    devo = dev
# E-step to update beta and lambda
    if(H==1){
      Pp = array(1,c(S,2^J,H))
    }else{
      Pp = array(0,c(S,2^J,H))
      for(s in 1:S) Pp[s,,] = (1/Pm[s,])*Pj[s,,]
    }
    Tmp = Y
    Tmp[,1] = Nv-nv
    Yh = array(0,c(S,2^J,H))
    for(h in 1:H) Yh[,,h] = Tmp*Pp[,,h]
# update beta
    if(H>1){
      Cl = apply(Yh,c(1,3),sum)
      tot = rowSums(Cl)
      sc = rep(0,np2)
      Fi = matrix(0,np2,np2)
      for(s in 1:S){
        Ls = as.matrix(L[s,,])
        sc = sc+t(Ls)%*%(Cl[s,]-tot[s]*Piv[s,])
        Om = diag(Piv[s,])-Piv[s,]%o%Piv[s,]
        Fi = Fi+tot[s]*t(Ls)%*%Om%*%Ls
      }
      if(rcond(Fi)>10^-15) dbeta = c(solve(Fi,sc)) else dbeta = c(ginv(Fi)%*%sc)
      mdbeta = max(abs(dbeta))
      if(mdbeta>0.1) dbeta = dbeta/mdbeta*0.25
      beta = beta+dbeta
    }
# update lambda
    sc = rep(0,np1)
    Fi = rep(0,np1,np1)
    for(h in 1:H) for(s in 1:S){
      if(is.null(X)) Mhs = M[,,h,1] else Mhs = M[,,h,s]
      Mhs = matrix(Mhs,dM1,dM2)
      if(model=="loglin"){
        sc = sc+t(Mhs)%*%(Yh[s,,h]-sum(Yh[s,,h])*Q[s,,h])
        Tmp = sum(Yh[s,,h])*(diag(Q[s,,h])-Q[s,,h]%o%Q[s,,h])
        Fi = Fi+t(Mhs)%*%Tmp%*%Mhs
      }
      if(model=="logit"){
        rho = exp(c(Mhs%*%lambda))
        rho = rho/(1+rho)
        sc = sc+t(Mhs)%*%(t(A)-rho*t(B))%*%Yh[s,,h]
        Tmp = sum(Yh[s,,h])*(diag(Q[s,,h])-Q[s,,h]%o%Q[s,,h])
        tmp = rho*(1-rho)*c(t(B)%*%Yh[s,,h])
        Fi = Fi+t(Mhs)%*%(tmp*Mhs)
      }
    }
    if(rcond(Fi)>10^-15) dlambda = c(solve(Fi,sc)) else dlambda = c(ginv(Fi)%*%sc)
    mdlambda = max(abs(dlambda))
    if(mdlambda>0.1) dlambda = dlambda/mdlambda*0.1
    lambda = lambda+dlambda
# compute probabilities
    if(H>1){
      Piv = matrix(0,S,H)
      for(s in 1:S){
        Ls = as.matrix(L[s,,])
        tmp = exp(Ls%*%beta)
        Piv[s,] = tmp/sum(tmp)
      }
    }
    Pj = Q = array(0,c(S,2^J,H))
    Pm = matrix(0,S,2^J)
    for(s in 1:S){
      for(h in 1:H){
        if(is.null(X)) Msh = as.matrix(M[,,h,1]) else Msh = as.matrix(M[,,h,s])
        if(model=="loglin"){
          tmp = exp(c(Msh%*%lambda))
          Q[s,,h] = tmp/sum(tmp)
        }
        if(model=="logit") Q[s,,h] = exp(A%*%Msh%*%lambda-B%*%log(1+exp(Msh%*%lambda)))
      }
      if(H>1){
        Pj[s,,] = Q[s,,]*(rep(1,2^J)%o%Piv[s,])
        Pm[s,] = Q[s,,]%*%Piv[s,]
      }
    }
    if(H==1){
      Pj = Q
      Pm = matrix(Q[,,1],S,2^J)
    }
    phiv = Pm[,1]
    Pm1 = Pm[,-1]
    Qm1 = (1/(1-phiv))*Pm1
# estimate N
    if(estN){
      Nv = nv/(1-phiv)
      N = sum(Nv)
    }else{
      it1 = 0
      tao = ta
      while((abs(ta-tao)/N>10^-10 & it1<10) | it1==0){
        it1 = it1+1
        tao = ta
        ta = N-sum(nv/(1-exp(psi)*phiv))
        d1 = sum(nv*exp(psi)*phiv/(1-exp(psi)*phiv)^2)
        tmp = ta/d1
        dtmp = abs(tmp)
        if(dtmp>0.1) tmp = tmp/dtmp*0.1
        psi = psi + tmp
      }
      Nv = nv/(1-exp(psi)*phiv)
    }
# compute log-likelihood and deviance
    lk = sum(Y1*log(Qm1))
    Yh = pmax(cbind(Nv-nv,Y1),10^-300)
    dev = 2*sum(Yh*log(Yh/(Nv*Pm)))
# display partial results
    if(verb & it%%10==0) if(estN){
      cat(sprintf("%11g", c(it,lk,lk-lko,N)), "\n", sep = " | ")
    }else{
      cat(sprintf("%11g", c(it,dev,dev-devo,N)), "\n", sep = " | ")
    }
  }
  if(verb){
    if(estN){
      if(it%%10>0) cat(sprintf("%11g", c(it,lk,lk-lko,N)), "\n", sep = " | ")
    }else{
      if(it%%10>0) cat(sprintf("%11g", c(it,dev,dev-devo,N)), "\n", sep = " | ")
    }
    cat("------------|-------------|-------------|-------------|\n")
  }
# names to the parameters
  if(H>1) names(beta) = beta_names
  names(lambda) = out$par_list

#---- standard errors ----
  if(se_out){
    th = c(beta,lambda)
    out = lkLCCRcon(th,np1,np2,model,H,J,S,L,M,X,nv,n,Y1,A,B,sc=TRUE)
    lke = out$lk
    st = out$st
    dphiv = out$dphiv
# compute numerical derivative
    tol = 10^-6
    DD = matrix(0,np1+np2,np1+np2)
    for(j in 1:(np1+np2)){
      th1 = th
      th1[j] = th1[j]+tol
      DD[,j] = (lkLCCRcon(th1,np1,np2,model,H,J,S,L,M,X,nv,n,Y1,A,B,sc=TRUE)$st-st)/tol
    }
    DD = (DD+t(DD))/2
    if(rcond(DD)>10^-15) VV = solve(-DD) else VV = ginv(-DD)
    tmp = diag(VV)
    if(any(tmp<0)) warning("negative diagonal elements in the observed information matrix")
    se = sqrt(tmp)
    if(H==1){
      sebeta = NULL
    }else{
      sebeta = se[(1:np2)]
      names(sebeta) = names(beta)
    }
    selambda = se[np2+(1:np1)]
    names(selambda) = names(lambda)
    seN = sqrt(t(nv)%*%((1/(1-phiv)^2)*dphiv)%*%VV%*%t((1/(1-phiv)^2)*dphiv)%*%nv)
    seN = c(seN)
  }

#---- output ----
  AIC = -2*lk+2*(np)
  BIC = -2*lk+log(n)*(np)
  out = list(beta=beta,lambda=lambda,lk=lk,dev=dev,N=N,np=np,AIC=AIC,BIC=BIC,M=M,phiv=phiv,Piv=Piv,Q=Q,
             call = match.call(),Y=Y,H=H,model=model,W=W,X=X,biv=biv,flag=flag,main=main,
             free_cov=free_cov,free_biv=free_biv,free_flag=free_flag,se_out=se_out)
  if(se_out){
    out$seN = seN
    out$sebeta = sebeta
    out$selambda = selambda
  }
  class(out) = "estLCCRcon"
  return(out)

}
