iic.lcp.groups <- function(n.groups,M,lambda.i.c,lambda.i.p,lambda.i){
  palette("R3")
  groups.legend <- NULL
  pdf(paste("iic_lcp.pdf",sep=""))
  plot.new()
  par(new=T,pty="s",bty="n")
  axis(side=1,labels=T,tick=T,pos=0,xaxp=c(0,1,5),mgp=c(2,0.5,0),cex.axis=1.5)
  axis(side=2,labels=T,tick=T,pos=0,yaxp=c(0,1,5),at=seq(0,1,0.2),cex.axis=1.5)
  axis(side=3,labels=F,pos=1,tck=0,at=seq(0,1,0.2))
  axis(side=4,labels=F,pos=1,tck=0,at=seq(0,1,0.2))
  title(mgp=c(2,0.5,0),xlab=expression(G[c](lambda)),ylab=expression(G[p](lambda)),cex.lab=2.2,cex.main=2.8)
  for(j in 1:n.groups){
    print(j)
    Gc <- Gp <- H <- NULL
    for(i in 1:M){	
      Gc[i] <- sum(lambda.i.c[,j]<=lambda.i.c[i,j])/M
      Gp[i] <- sum(lambda.i.p[,j]<=lambda.i.p[i,j])/M
      H[i] <- sum(lambda.i[,j]<=lambda.i[i,j])/M
    }
    plot.Gc <- Gc[order(Gc)]
    plot.Gp <- vector("numeric",M)
    indices <- findInterval(lambda.i.c[order(Gc),j],lambda.i.p[order(Gp),j])
    plot.Gp[indices==0] <- 0
    plot.Gp[indices!=0] <- Gp[order(Gp)][indices]
    points(plot.Gc,plot.Gp,type="l",col=(j+1))
    points(plot.Gc[findInterval(lambda.i[order(H),j][findInterval(quantile(lambda.i[,j],probs=seq(0.1,0.9,0.1)),lambda.i[order(H),j])],lambda.i.c[order(Gc),j])],plot.Gp[findInterval(lambda.i[order(H),j][findInterval(quantile(lambda.i[,j],probs=seq(0.1,0.9,0.1)),lambda.i[order(H),j])],lambda.i.c[order(Gc),j])],pch=8,col=(j+1))
    groups.legend <- c(groups.legend,paste("Group",j))
    }
  legend(0.05,0.9,groups.legend,col=seq(2,n.groups+1),lty=1)
  dev.off()
}