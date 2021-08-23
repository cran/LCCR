design_matrix_loglin <-
function(J,H=1,main=c("LC","same","Rasch"),
                                 X=NULL,free_cov=c("no","class","resp","both"),
                                 biv=NULL,free_biv=c("no","class","int","both")){

# ---- preliminaries ----
  if(J<2) stop("J must be at least 2")
  main = match.arg(main)
  free_cov = match.arg(free_cov)
  free_biv = match.arg(free_biv)
  if(is.null(X)){
    S = 1; ncov = 0
  }else{
    if(is.vector(X)){
      S = length(X); ncov = 1
      X = array(X,c(S,1,J))
    }
    if(is.matrix(X)){
      S = nrow(X); ncov = ncol(X)
      X = array(X,c(S,ncov,J))
    }
    if(is.array(X)){
      S = dim(X)[1]; ncov = dim(X)[2]
      if(dim(X)[3]!=J) stop("wrong dimension of covariate array")
    }
  }
  if(is.null(biv)){
    nbiv = 0
  }else{
    if(is.vector(biv)) biv = t(biv)
    nbiv = nrow(biv)
  }

# ---- main and bivariate effects ----
  Main = sq(J)
  if(nbiv==0){
    Int = NULL
  }else{
    Int = matrix(0,2^J,nbiv)
    for(j in 1:nbiv) Int[,j] = Main[,biv[j,1]]*Main[,biv[j,2]]
  }

# ----- number and list of parameters ----
# main effects
  if(main=="LC"){
    np = J*H
    par_list = NULL
    if(H==1){
      par_list = paste("main",1:J,sep="")
    }else{
      for(h in 1:H) par_list = c(par_list,paste("class",h,".main",1:J,sep=""))
    }
  }
  if(main=="same"){
    np = H
    par_list = NULL
    if(H==1){
      par_list = "main"
    }else{
      for(h in 1:H) par_list = c(par_list,paste("class",h,".main",sep=""))
    }
  }
  if(main=="Rasch"){
    np = J+(H-1)
    if(H==1){
      par_list = NULL
    }else{
      par_list = paste("class",2:H,sep="")
    }
    par_list = c(par_list,paste("main",1:J,sep=""))
  }
# covariate effects
  if(ncov>0){
    if(free_cov=="no"){
      np = np+ncov
      par_list = c(par_list,paste("cov",1:ncov,sep=""))
    }
    if(free_cov=="class"){
      np = np+H*ncov
      if(H==1){
        par_list = c(par_list,paste("cov",1:ncov,sep=""))
      }else{
        for(h in 1:H) par_list = c(par_list,paste("class",h,".cov",1:ncov,sep=""))
      }
    }
    if(free_cov=="resp"){
      np = np+J*ncov
      for(j in 1:J) par_list = c(par_list,paste("resp",j,".cov",1:ncov,sep=""))
    }
    if(free_cov=="both"){
      np = np+H*J*ncov
      if(H==1){
        for(j in 1:J) par_list = c(par_list,paste("resp",j,".cov",1:ncov,sep=""))
      }else{
        for(h in 1:H) for(j in 1:J) par_list = c(par_list,paste("class",h,".resp",j,".cov",1:ncov,sep=""))
      }
    }
  }
# interaction effect
  if(nbiv>0){
    if(free_biv=="no"){
      np = np+1
      par_list = c(par_list,"biv")
    }
    if(free_biv=="class"){
      np = np+H
      if(H==1){
        par_list = c(par_list,"biv")
      }else{
        for(h in 1:H) par_list = c(par_list,paste("class",h,".biv",sep=""))
      }
    }
    if(free_biv=="int"){
      np = np+nbiv
      for(j in 1:nbiv) par_list = c(par_list,paste("biv",biv[j,1],"-",biv[j,2],sep=""))
    }
    if(free_biv=="both"){
      np = np+H*nbiv
      if(H==1){
        for(j in 1:nbiv) par_list = c(par_list,paste("biv",biv[j,1],"-",biv[j,2],sep=""))
      }else{
        for(h in 1:H) for(j in 1:nbiv) par_list = c(par_list,paste("class",h,".biv",biv[j,1],"-",biv[j,2],sep=""))
      }
    }
  }

# ----- design matrices ----
  M = array(0,c(2^J,np,H,S))
  for(s in 1:S){
    for(h in 1:H){
# main effects
      tmp = matrix(0,1,H); tmp[h] = 1
      if(main=="LC") Tmp1 = tmp%x%Main
      if(main=="same") Tmp1 = tmp%x%as.matrix(rowSums(Main))
      if(main=="Rasch") Tmp1 = cbind(rep(1,2^J)%o%c(tmp[-1]),Main)
# covariate effects
      Tmp2 = NULL
      if(ncov>0){
        Xs = matrix(X[s,,],ncov,J)
        if(free_cov=="no") Tmp2 = Main%*%t(Xs)
        if(free_cov=="class") Tmp2 = tmp%x%(Main%*%t(Xs))
        if(free_cov=="resp"){
          Tmp2 = NULL
          for(j in 1:J) Tmp2 = cbind(Tmp2,Main[,j]%o%Xs[,j])
        }
        if(free_cov=="both"){
          Tmp2 = NULL
          for(j in 1:J) Tmp2 = cbind(Tmp2,Main[,j]%o%Xs[,j])
          Tmp2 = tmp%x%Tmp2
        }
      }
# interaction effects
      Tmp4 = NULL
      if(nbiv>0){
        if(free_biv=="no") Tmp4 = rowSums(Int)
        if(free_biv=="class") Tmp4 = tmp%x%rowSums(Int)
        if(free_biv=="int") Tmp4 = Int
        if(free_biv=="both") Tmp4 = tmp%x%Int
      }
      M[,,h,s] = cbind(Tmp1,Tmp2,Tmp4)
    }
  }
  colnames(M) = par_list

#---- output ----
  out = list(M=M,par_list=par_list,Main=Main)

}
