freq_data = function(R,count=rep(1,nrow(R))){

#---- convert capture histories in frequencies ----
  n = nrow(R)
  J = ncol(R)
  ss = 2^((J-1):0)
  ind = R%*%ss+1
  Y = matrix(0,n,2^J)
  Y[cbind(1:n,ind)] = count

#---- output ----
  return(Y)

}
