sq <-
function(J,t=NULL){

# generate matrix
  M = NULL
  if(J>=1){
    if(is.null(t)){   # se non viene dato t
      if(J==1){
        M = rbind(0,1)
      }else{
        TT = sq(J-1);  nt = nrow(TT)
        M = rbind(cbind(0,TT),cbind(1,TT));
      }
    }else{
      if(t==J){          # se sono tutti elementi 1
        M = matrix(1,1,J)
      }else if(t>1){       # se sono sia 1 che 0
        for(i in 1:(J-t+1)){
          S = sq(J-i,t-1)
          r = nrow(S)
          TT = cbind(matrix(0,r,i-1),1,S)
          M = rbind(TT,M)
        }
      }else if(t==1){
        M = matrix(0,J,J)
        for(j in 1:J) M[j,J-j+1] = 1
      }else{
        M = matrix(0,1,J)   # se non ci sono elementi 1
      }
    }
  }

# output
  return(M)

}
