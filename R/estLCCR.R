estLCCR <-
function(Y,H,model=c("loglin","logit"),W=NULL,X=NULL,N=NULL,biv=NULL,
                    flag=c("no","prev","sum","atleast"),
                    main=c("LC","same","Rasch"),
                    free_cov=c("no","class","resp","both"),
                    free_biv=c("no","class","int","both"),
                    free_flag=c("no","class","resp","both"),N0=NULL,be0=NULL,la0=NULL,
                    maxit=5000,verb=TRUE,init_rand=FALSE,se_out=FALSE){

# ---- preliminaries ----
  model = match.arg(model)
  flag = match.arg(flag)
  main = match.arg(main)
  free_cov = match.arg(free_cov)
  free_biv = match.arg(free_biv)
  free_flag = match.arg(free_flag)
  estN = is.null(N)
  S = nrow(Y)
  J = round(log(ncol(Y))/log(2))
  Y1 = Y[,-1]
  nv = rowSums(Y)
  if(is.null(W)){
    ncov2 = 0
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
    be_names = NULL
    for(h in 2:H) if(ncov2>0){
      be_names = c(be_names,paste("class",h,".int",sep=""),paste("class",h,".",Wnames,sep=""))
    }else{
      be_names = c(be_names,paste("class",h,".int",sep=""))
    }
  }
  if(model=="loglin") out = design_matrix_loglin(J,H,main,X,free_cov,biv,free_biv)
  if(model=="logit") out = design_matrix_logit(J,H,main,X,free_cov,flag,free_flag)
  L = array(0,c(S,H,(H-1)*(1+ncov2)))
  Tmp = diag(H)[,-1,drop=FALSE]
  for(s in 1:S) L[s,,] = Tmp%x%t(c(1,W[s,]))
  if(model=="logit"){
    A = out$A; B = out$B
  }
  M = out$M
  np1 = length(out$par_list)
  np = np1+np2

#---- starting values ----
  n = sum(nv)
  if(estN){
    if(is.null(N0)){
      if(init_rand) N = n*(1+runif(1)) else N = n*1.25 
    }else{
      N = N0
    }
  }
  if(is.null(la0)){
    if(init_rand){
      la = rnorm(np1)
    }else{
      if(H==1){
        la = rep(0,np1)
      }else{
        est = estLCCR(Y,1,model=model,main=main,verb=verb)
        if(main=="LC" || main=="same"){
          la = NULL
          for(h in 1:H) la = c(la,est$la+h-H/2)
        }
        if(main=="Rasch") la = c(1:(H-1),est$la-(H-1)/2)
        la = c(la,rep(0,np1-length(la)))
      }
    }
  }else{
    la = la0
  }
  if(is.null(be0)){
    if(H==1) be = NULL else if(init_rand) be = rnorm(np2) else be = rep(0,np2)
  }else{
    be = be0
  }
  if(H==1){
    Piv = matrix(1,S,1)
  }else{
    Piv = matrix(0,S,H)
    for(s in 1:S){
      Ls = as.matrix(L[s,,])
      tmp = exp(Ls%*%be)
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
        tmp = exp(c(Msh%*%la))
        Q[s,,h] = tmp/sum(tmp)
      }
      if(model=="logit") Q[s,,h] = exp(A%*%Msh%*%la-B%*%log(1+exp(Msh%*%la)))
    }
    if(H>1){
      Pj[s,,] = Q[s,,]*(rep(1,2^J)%o%Piv[s,])
      Pm[s,] = Q[s,,]%*%Piv[s,]
    }
  }
  if(H==1){
    Pj = Q
    Pm = Q[,,1]
  }
  phiv = Pm[,1]
  tauv = nv/n
  phi = sum(tauv*phiv)
  Pm1 = Pm[,-1]
# compute log-likelhood
  lk = lgamma(N+1)-lgamma(N-n+1)+(N-n)*log(phi)+sum(Y1*log(Pm1))+sum(nv*log(tauv))
  if(verb){
    cat("------------|-------------|-------------|-------------|\n")
    cat("  iteration |      lk     |    lk-lko   |      N      |\n")
    cat("------------|-------------|-------------|-------------|\n")
    cat(sprintf("%11g", c(0,lk,NA,N)), "\n", sep = " | ")
  }

# ---- cycle ----
  it = 0; lko = lk
  while((abs(lk-lko)/abs(lko)>10^-10 | it==0) & it<maxit){
    it = it+1; lko = lk
# E-step to update be and La
    if(H==1){
      Pp = array(1,c(S,2^J,H))
    }else{
      Pp = array(0,c(S,2^J,H))
      for(s in 1:S) Pp[s,,] = (1/Pm[s,])*Pj[s,,]
    }
    Tmp = Y; Tmp[,1] = (N-n)*phiv*tauv/phi
    Yh = array(0,c(S,2^J,H))
    for(h in 1:H) Yh[,,h] = Tmp*Pp[,,h]
# update be
    if(H>1){
      Cl = apply(Yh,c(1,3),sum)
      tot = rowSums(Cl)
      sc = rep(0,np2); Fi = matrix(0,np2,np2)
      for(s in 1:S){
        Ls = as.matrix(L[s,,])
        sc = sc+t(Ls)%*%(Cl[s,]-tot[s]*Piv[s,])
        Om = diag(Piv[s,])-Piv[s,]%o%Piv[s,]
        Fi = Fi+tot[s]*t(Ls)%*%Om%*%Ls
      }
      dbe = c(solve(Fi,sc))
      mdbe = max(abs(dbe))
      if(mdbe>0.1) dbe = dbe/mdbe*0.25
      be = be+dbe
    }
# update La
    sc = rep(0,np1); Fi = rep(0,np1,np1)
    for(h in 1:H) for(s in 1:S){
      if(is.null(X)) Mhs = M[,,h,1] else Mhs = M[,,h,s]
      if(model=="loglin"){
        sc = sc+t(Mhs)%*%(Yh[s,,h]-sum(Yh[s,,h])*Q[s,,h])
        Tmp = sum(Yh[s,,h])*(diag(Q[s,,h])-Q[s,,h]%o%Q[s,,h])
        Fi = Fi+t(Mhs)%*%Tmp%*%Mhs
      }
      if(model=="logit"){
        rho = exp(c(Mhs%*%la)); rho = rho/(1+rho)
        sc = sc+t(Mhs)%*%(t(A)-rho*t(B))%*%Yh[s,,h]
        Tmp = sum(Yh[s,,h])*(diag(Q[s,,h])-Q[s,,h]%o%Q[s,,h])
        tmp = rho*(1-rho)*c(t(B)%*%Yh[s,,h])
        Fi = Fi+t(Mhs)%*%(tmp*Mhs)
      }
    }
    if(rcond(Fi)>10^-15) dla = c(solve(Fi,sc)) else dla = c(ginv(Fi)%*%sc)
    mdla = max(abs(dla))
    if(mdla>0.1) dla = dla/mdla*0.1
    la = la+dla
# compute probabilities
    if(H>1){
      Piv = matrix(0,S,H)
      for(s in 1:S){
        Ls = as.matrix(L[s,,])
        tmp = exp(Ls%*%be)
        Piv[s,] = tmp/sum(tmp)
      }
    }
    Pj = Q = array(0,c(S,2^J,H))
    Pm = matrix(0,S,2^J)
    for(s in 1:S){
      for(h in 1:H){
        if(is.null(X)) Msh = as.matrix(M[,,h,1]) else Msh = as.matrix(M[,,h,s])
        if(model=="loglin"){
          tmp = exp(c(Msh%*%la))
          Q[s,,h] = tmp/sum(tmp)
        }
        if(model=="logit") Q[s,,h] = exp(A%*%Msh%*%la-B%*%log(1+exp(Msh%*%la)))
      }
      if(H>1){
        Pj[s,,] = Q[s,,]*(rep(1,2^J)%o%Piv[s,])
        Pm[s,] = Q[s,,]%*%Piv[s,]
      }
    }
    if(H==1){
      Pj = Q
      Pm = Q[,,1]
    }
    phiv = Pm[,1]
    tauv = nv/n
    phi = sum(tauv*phiv)
    Pm1 = Pm[,-1]
# update tau
    dtauv = Inf
    while(max(abs(dtauv))>10^-10){
      dtauv = (nv+(N-n)*phiv*tauv/phi)/N-tauv
      tauv = tauv+dtauv
      phi = c(tauv%*%phiv)
    }
# update N
    if(estN){
      dN = Inf
      while(abs(dN)/N>10^-10){
        d1 = digamma(N+1)-digamma(N-n+1)+log(phi)
        d2 = trigamma(N+1)-trigamma(N-n+1)
        dN = -d1/d2/max(1,it-10)
        while(N+dN<n) dN = dN/2
        N = N+dN
      }
    }
# compute log-likelihood
    lk1 = lgamma(N+1)-lgamma(N-n+1)
    lk2 = (N-n)*log(phi)
    lk3 = sum(Y1*log(Pm1))
    lk4 = sum(nv*log(tauv))
    lk = lk1+lk2+lk3+lk4
# display partial results
    if(verb & it%%10==0) cat(sprintf("%11g", c(it,lk,lk-lko,N)), "\n", sep = " | ")
  }
  if(verb){
    if(it%%10>0) cat(sprintf("%11g", c(it,lk,lk-lko,N)), "\n", sep = " | ")
    cat("------------|-------------|-------------|-------------|\n")
  }
  if(H>1) names(be) = be_names
  names(la) = out$par_list


#---- standard errors ----
  if(se_out){
    th = c(N,be,la)
    out = lkLCCR(th,np1,np2,model,H,J,S,L,M,tauv,X,nv,n,Y1,A,B,sc=TRUE)
    lke = out$lk; st = out$st
# compute numerical derivative
    tol = 10^-6
    DD = matrix(0,1+np1+np2,1+np1+np2)
    for(j in 1:(1+np1+np2)){
      th1 = th; th1[j] = th1[j]+tol
      DD[,j] = (lkLCCR(th1,np1,np2,model,H,J,S,L,M,tauv,X,nv,n,Y1,A,B,sc=TRUE)$st-st)/tol
    }
    DD = (DD+t(DD))/2
    if(rcond(DD)>10^-15) tmp = diag(solve(-DD)) else tmp = diag(ginv(-DD))
    if(any(tmp<0)) warning("negative diagonal elements in the observed information matrix")
    se = sqrt(tmp)
    seN = se[1]
    if(H==1){
      sebe = NULL
    }else{
      sebe = se[1+(1:np2)]
      names(sebe) = names(be)
    }
    sela = se[1+np2+(1:np1)]
    names(sela) = names(la)
  }

#---- output ----
  AIC = -2*lk+2*(np)
  BIC = -2*lk+log(n)*(np)
  out = list(be=be,la=la,lk=lk,N=N,np=np,AIC=AIC,BIC=BIC,M=M,tauv=tauv,phiv=phiv,
             lk1=lk1,lk2=lk2,lk3=lk3,lk4=lk4,call = match.call(),
             Y=Y,H=H,model=model,W=W,X=X,biv=biv,flag=flag,main=main,
             free_cov=free_cov,free_biv=free_biv,free_flag=free_flag,se_out=se_out)
  if(se_out){
    out$seN = seN; out$sebe = sebe; out$sela = sela
  }
  class(out) = "estLCCR"
  return(out)

}
