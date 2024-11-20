plot.confLCCRcon = function(x,...){

#---- plot of confidence interval starting from output of confLCCRcon ---
  mdev = min(x$devv)
  plot(x$Nv,x$devv,type="l",xlab=expression(italic(N)),ylab="Deviance",col=4)
  lines(c(x$conf[1],x$conf[1]),c(mdev,x$dev1),col=2)
  lines(c(x$conf[2],x$conf[2]),c(mdev,x$dev2),col=2)
  lines(c(x$conf),c(x$dev1,x$dev2),col=2)
  lines(c(x$Nh,x$Nh),c(mean(x$dev1,x$dev2),x$devh),col=2)
  text(x$conf[1],x$dev1,round(x$conf[1],0),pos=3)
  text(x$conf[2],x$dev2,round(x$conf[2],0),pos=3)
  text(x$Nh,mean(x$dev1,x$dev2),round(x$Nh,0),pos=3)

}
