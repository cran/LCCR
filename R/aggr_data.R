aggr_data <-
function(Y,W=NULL,X=NULL){
  
# ---- aggregate ----
  if(is.null(W) & is.null(X)){
    Ya = t(colSums(Y))
  }else{
    n = nrow(Y); J = ncol(Y)
    if(is.null(W)){
      nc1 = 0
    }else{
      if(is.vector(W)) nc1 = 1
      if(is.matrix(W)) nc1 = ncol(W)
    }
    if(is.null(X)){
      nc2 = 0
      XX = NULL
    }else{
      if(is.vector(X)){
        nc2 = 1
        XX = X
      }
      if(is.matrix(X)){
        nc2 = ncol(X)
        XX = X
      }
      if(is.array(X) & !is.matrix(X)){
        di = dim(X)
        nc2 = di[2]*di[3]
        XX = matrix(X,di[1],di[2]*di[3])
      }
    }
    WX = cbind(W,XX)
    WXa = unique(WX)
    I = nrow(WXa)
    Ya = matrix(0,I,J)
    for(i in 1:n){
      ind = apply(WXa,1,function(x) identical(x,WX[i,]))
      Ya[ind,] = Ya[ind,]+Y[i,]
    }
    if(nc1>0){
      Wa = WXa[,1:nc1]
      if(is.matrix(W)) Wa = as.matrix(Wa)
    }
    if(nc2>0){
      Xa = WXa[,nc1+(1:nc2)]
      if(is.matrix(X)) Xa = as.matrix(Xa)
      if(is.array(X) & !is.matrix(X)) Xa = array(Xa,c(I,di[2],di[3]))
    }
    
  }

# ---- output ----
  out = list(Ya=Ya)
  if(!is.null(W)) out$Wa=Wa
  if(!is.null(X)) out$Xa=Xa
  return(out)

}
