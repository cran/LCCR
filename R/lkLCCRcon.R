lkLCCRcon <- function(th,np1,np2,model,H,J,S,L,M,X,nv,n,Y1,A,B,sc=FALSE){

  # conditional likelihood version

# separate parameters
  be = th[1:np2]
  la = th[np2+(1:np1)]

# compute probabilities
  if(H>1){
    Piv = matrix(0,S,H)
    for(s in 1:S){
      Ls = as.matrix(L[s,,])
      tmp = exp(Ls%*%be)
      Piv[s,] = tmp/sum(tmp)
    }
  }
  Q = array(0,c(S,2^J,H))
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
    if(H>1) Pm[s,] = Q[s,,]%*%Piv[s,]
  }
  if(H==1) Pm = Q[,,1]
  phiv = Pm[,1]
  Pm1 = Pm[,-1]
  Qm1 = (1/(1-phiv))*Pm1

# compute log-likelihood
  lk = sum(Y1*log(Qm1))
  out = lk

# compute score
  if(sc){
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
          tmp = c(exp(Msh%*%la)/(1+exp(Msh%*%la)))
          dQ[s,,h,np2+(1:np1)] = Q[s,,h]*(A%*%Msh-B%*%(tmp*Msh))
        }
      }
      if(H>1) for(j in 1:(np1+np2)) dPm[s,,j] = dQ[s,,,j]%*%Piv[s,]+Q[s,,]%*%dPiv[s,,j]
    }
    if(H==1) dPm = array(dQ[,,1,],c(S,2^J,np1+np2))
    dphiv = dPm[,1,]
    dPm1 = array(dPm[,-1,],c(S,2^J-1,np1+np2))
    dQm1 = array(0,c(S,2^J-1,np1+np2))
    for(s in 1:S) for(j in 1:(2^J-1)){
      dQm1[s,j,] = (1/(1-phiv[s])^2)*dphiv[s,]*Pm1[s,j]+(1/(1-phiv[s]))*dPm1[s,j,]
    }

# compute score with respect to be and la
    Dyb = 0
    for(s in 1:S) Dyb = Dyb+t(dQm1[s,,])%*%(Y1[s,]/Qm1[s,])

#output 
    out = list(lk=lk,st=Dyb,dphiv=dphiv)
  }

# output
  return(out)
  
}