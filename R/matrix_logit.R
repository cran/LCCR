matrix_logit = function(J){

#---- preliminaries ----
  S = sq(J)

#---- build desisgn matrices ----
  A = as.matrix(S[,1])
  B = matrix(1,2^J,1)
  if(J>1) for(j in 2:J){
    if(j==2) ind = S[,1:(j-1)]+1 else ind = S[,1:(j-1)]%*%(2^((j-2):0))+1
    Tmp = matrix(0,2^J,2^(j-1))
    Tmp[cbind(1:(2^J),ind)] = 1
    A = cbind(A,S[,j]*Tmp)
    B = cbind(B,Tmp)
  }

#---- output ----
  out = list(A=A,B=B)

}
