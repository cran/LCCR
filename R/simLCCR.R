simLCCR <-
function(H,J,be,la,N,model=c("loglin","logit"),Wc=NULL,Xc=NULL,
                    biv=NULL,flag=c("no","prev","sum","atleast"),
                    main=c("LC","same","Rasch"),
                    free_cov=c("no","class","resp","both"),
                    free_biv=c("no","class","int","both"),
                    free_flag=c("no","class","resp","both")){

# ---- preliminaries ----
  model = match.arg(model)
  flag = match.arg(flag)
  main = match.arg(main)
  free_cov = match.arg(free_cov)
  free_biv = match.arg(free_biv)
  free_flag = match.arg(free_flag)
  if(length(N)>1){
    S = length(N)
    aggr_data = TRUE
  }else{
    S = N
    aggr_data = FALSE
  }
  if(is.null(Wc)){
    ncov2 = 0
  }else{
    if(is.vector(Wc)) Wc = as.matrix(Wc)
    ncov2 = ncol(Wc)
  }
  np2 = (H-1)*(1+ncov2)
  if(model=="loglin") out = design_matrix_loglin(J,H,main,Xc,free_cov,biv,free_biv)
  if(model=="logit") out = design_matrix_logit(J,H,main,Xc,free_cov,flag,free_flag)
  L = array(0,c(S,H,np2))
  Tmp = diag(H)[,-1,drop=FALSE]
  if(ncov2==0){
    for(s in 1:S) L[s,,] = Tmp
  }else{
    if(ncov2>0) for(s in 1:S) L[s,,] = Tmp%x%t(c(1,Wc[s,]))
  }
  if(model=="logit"){
    A = out$A; B = out$B
  }
  M = out$M
  np1 = length(out$par_list)
  np = np1+np2

#---- probabilities ----
# class weights
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
# response probabilities
  Q = array(0,c(S,2^J,H))
  if(aggr_data) Pm = matrix(0,S,2^J)
  for(s in 1:S){
    for(h in 1:H){
      if(is.null(Xc)) Msh = as.matrix(M[,,h,1]) else Msh = as.matrix(M[,,h,s])
      if(model=="loglin"){
        tmp = exp(c(Msh%*%la))
        Q[s,,h] = tmp/sum(tmp)
      }
      if(model=="logit") Q[s,,h] = exp(A%*%Msh%*%la-B%*%log(1+exp(Msh%*%la)))
    }
    if(aggr_data) Pm[s,] = Q[s,,]%*%Piv[s,]
  }

# ---- simulate ----
# individual data
  if(!aggr_data){
    SQ = sq(J)
    Rc = matrix(0,S,J)
    Uc = rep(0,N)
  }
  Yc = matrix(0,S,2^J)
  if(aggr_data){
    for(i in 1:S) Yc[i,] = rmultinom(1,N[i],Pm[i,])
    Y = Yc; Y[,1] = 0
    ind = which(rowSums(Y[,-1])>0)
    Y = Y[ind,]
  }else{
    for(i in 1:S){
      Uc[i] = which(rmultinom(1,1,Piv[i,])==1)
      Yc[i,] = rmultinom(1,1,Q[i,,Uc[i]])
      Rc[i,] = SQ[Yc[i,]==1]
    }
    ind = which(Yc[,1]==0)
    U = Uc[ind]
    Y = Yc[ind,]
    R = Rc[ind,]
  }
  if(!is.null(Wc)){
    if(is.vector(Wc)) W = Wc[ind]
    if(is.matrix(Wc)) W = Wc[ind,]
  }
  if(!is.null(Xc)){
    if(is.vector(Xc)) X = Xc[ind]
    if(is.matrix(Xc)) X = Xc[ind,]
  }

  #---- output ----
  out = list(Y=Y,Yc=Yc,Piv=Piv,Q=Q)
  if(aggr_data){
    Pm = Pm
  }else{
    out$R=R; out$U=U; out$Rc=Rc
  }
  if(!is.null(Wc)) out$W = W
  if(!is.null(Xc)) out$X = X
  return(out)

}
