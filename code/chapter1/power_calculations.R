#####power to detect main snp effects associated with survival#####################
get.delta=function(PA, propEvents, N, alpha, beta)
  exp((qnorm(1-alpha)+qnorm(beta))/(sqrt(N*propEvents*PA*(1-PA))))
propEvents=seq(0.1, 0.7, by=0.01)
delta=mapply(get.delta, 
             propEvents,
             MoreArgs=list(PA=0.40, N=2580, alpha=0.05/1000000, beta=0.20))
delta.lb=mapply(get.delta,
                propEvents,
                MoreArgs=list(PA=0.05,
                              N=2580,
                              alpha=0.05/1000000,
                              beta=0.20))
plot(propEvents,
     delta,type="l",
     axes=F,
     ylim=c(0,4.0)
     ,lwd=2,
     main=c("Power Calculations"),
     xlab="Proportion of Events",
     ylab=c("Minimum Detectable", "Hazard Ratio"))
abline(h=seq(0,4.0,by=0.25),col="gray72")
abline(v=seq(0.1,0.7,by=0.025),col="gray72")
axis(1);axis(2)
lines(propEvents,delta,lwd=2)
lines(propEvents,delta.lb,lwd=2,lty=2)
legend(0.276,
       4.0,
       lty=1:2,
       legend=c("MAF=0.40","MAF=0.05"),
       lwd=2,
       cex=1.2,
       bty="o",
       bg="white")
