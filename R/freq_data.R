freq_data <-
function(R){

  n = nrow(R)
  J = ncol(R)
  ss = 2^((J-1):0)
  ind = R%*%ss+1
  Y = matrix(0,n,2^J)
  Y[cbind(1:n,ind)] = 1
  return(Y)
  
}
