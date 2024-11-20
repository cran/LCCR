plot.confLCCR = function(x,...){

#---- plot of confidence interval starting from output of confLCCR ---
  mlk = min(x$lkv)
  plot(x$Nv,x$lkv,type="l",xlab=expression(italic(N)),ylab="LogLik",col=4)
  lines(c(x$conf[1],x$conf[1]),c(mlk,x$lk1),col=2)
  lines(c(x$conf[2],x$conf[2]),c(mlk,x$lk2),col=2)
  lines(c(x$conf),c(x$lk1,x$lk2),col=2)
  lines(c(x$Nh,x$Nh),c(mlk,x$lkh),col=2)
  text(x$conf[1],x$lk1,round(x$conf[1],0),pos=3)
  text(x$conf[2],x$lk2,round(x$conf[2],0),pos=3)
  text(x$Nh,mean(x$lk1,x$lk2),round(x$Nh,0),pos=3)
  
}
