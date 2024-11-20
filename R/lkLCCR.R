lkLCCR = function(th,np1,np2,model,H,J,S,L,M,tauv,X,nv,n,Y1,A,B,sc=FALSE){

# Internal function that compute the log-likelihood and the corresponding score used by
# function estLCCR

#---- separate parameters ----
  N = th[1]
  beta = th[1+(1:np2)]
  lambda = th[1+np2+(1:np1)]

#---- compute probabilities ----
  if(H>1){
    Piv = matrix(0,S,H)
    for(s in 1:S){
      Ls = as.matrix(L[s,,])
      tmp = exp(Ls%*%beta)
      Piv[s,] = tmp/sum(tmp)
    }
  }
  Q = array(0,c(S,2^J,H))
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
    if(H>1) Pm[s,] = Q[s,,]%*%Piv[s,]
  }
  if(H==1) Pm = matrix(Q[,,1],S,2^J)
  phiv = Pm[,1]
  phi = sum(tauv*phiv)
  Pm1 = matrix(Pm[,-1],S,2^J-1)

#---- update tau ----
  dtauv = Inf
  while(max(abs(dtauv))>10^-10){
    dtauv = (nv+(N-n)*phiv*tauv/phi)/N-tauv
    tauv = tauv+dtauv
    phi = c(tauv%*%phiv)
  }

#---- compute log-likelihood ----
  lk1 = lgamma(N+1)-lgamma(N-n+1)
  lk2 = (N-n)*log(phi)
  lk3 = sum(Y1*log(Pm1))
  lk4 = sum(nv*log(tauv))
  lk = lk1+lk2+lk3+lk4
  out = lk

#----- compute score ----
  if(sc){
    # compute score with respect to N
    sN = digamma(N+1)-digamma(N-n+1)+log(phi)
# compute derivatives of the probabilities
    if(H>1){
      dPiv = array(0,c(S,H,np1+np2))
      for(s in 1:S){
        Ls = as.matrix(L[s,,])
        dPiv[s,,1:np2] = (diag(Piv[s,])-Piv[s,]%o%Piv[s,])%*%Ls
      }
    }
    dQ = array(0,c(S,2^J,H,np1+np2))
    dPm = array(0,c(S,2^J,np1+np2))
    for(s in 1:S){
      for(h in 1:H){
        if(is.null(X)) Msh = as.matrix(M[,,h,1]) else Msh = as.matrix(M[,,h,s])
        if(model=="loglin") dQ[s,,h,np2+(1:np1)] = (diag(Q[s,,h])-Q[s,,h]%o%Q[s,,h])%*%Msh
        if(model=="logit"){
          tmp = c(exp(Msh%*%lambda)/(1+exp(Msh%*%lambda)))
          dQ[s,,h,np2+(1:np1)] = Q[s,,h]*(A%*%Msh-B%*%(tmp*Msh))
        }
      }
      if(H>1) for(j in 1:(np1+np2)) dPm[s,,j] = dQ[s,,,j]%*%Piv[s,]+Q[s,,]%*%dPiv[s,,j]
    }
    if(H==1) dPm = array(dQ[,,1,],c(S,2^J,np1+np2))
    dphiv = matrix(dPm[,1,],S,np1+np2)
    dPm1 = array(dPm[,-1,],c(S,2^J-1,np1+np2))
# compute score with respect to beta and lambda
    Dyb = 0
    for(s in 1:S) Dyb = Dyb+t(dPm1[s,,])%*%(Y1[s,]/Pm1[s,])
    sp = (N-n)/phi*t(dphiv)%*%tauv+Dyb
# overall score
    st = c(sN,sp)

#---- output ----
    out = list(lk=lk,st=st)
  }

# output
  return(out)
  
}