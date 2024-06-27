mcmc <- function(data,M,lambda.0,mu.0,m.mu,s2.mu,nu.sigma,s.2.sigma,nu.tau,s.2.tau){
  for(i in 1:M){
    cat("Iteration number:",i,"\n")
    k <- dim(data)[1]
    n <- dim(data)[2]
    lambda <- lambda.0
    mu <- mu.0
    tau2 <- 1/rgamma(1,shape=(nu.tau+k)/2,rate=(nu.tau*s.2.tau+sum((lambda-mu)^2))/2)
    sigma2 <- 1/rgamma(1,shape=(nu.sigma+k*n)/2,rate=(nu.sigma*s.2.sigma+sum((data-lambda)^2))/2)
    mu <- rnorm(1,mean=((m.mu/s2.mu)+(k*mean(lambda)/tau2))/((1/s2.mu)+(k/tau2)),sd=sqrt(1/((1/s2.mu)+(k/tau2))))
    lambda <- rnorm(k,mean=((mu/tau2)+(n*apply(data,1,mean)/sigma2))/((1/tau2)+(n/sigma2)),sd=sqrt(rep(1/((1/tau2)+(n/sigma2)),k)))
    write(lambda,"trace.lambda.txt",append=T,ncol=k)
  }
}
mcmc.c <- function(alldata,M,lambda.0,nu.sigma,s.2.sigma){
  m <- dim(alldata)[1]
  for(j in 1:m){
    cat("Lambda number:",j,"\n")
    data <- alldata[j,,drop=F]
    lambda <- lambda.0
    for(i in 1:M){
      k <- dim(data)[1]
      n <- dim(data)[2]
      sigma2 <- 1/rgamma(1,shape=(nu.sigma+k*n)/2,rate=(nu.sigma*s.2.sigma+sum((data-lambda)^2))/2)
      lambda <- rnorm(1,mean=apply(data,1,mean),sd=sqrt(sigma2/n))
      write(lambda,paste("trace_trace.lambda.c_",j,sep=""),append=T)
    }
  }
}

mcmc.p <- function(alldata,M,lambda.0,mu.0,m.mu,s2.mu,nu.sigma,s.2.sigma,nu.tau,s.2.tau){
  m <- dim(alldata)[1]
  for(j in 1:m){
    cat("Lambda number:",j,"\n")
    data <- alldata[-j,]
    lambda <- lambda.0
    mu <- mu.0
    for(i in 1:M){
      k <- dim(data)[1]
      n <- dim(data)[2]
      tau2 <- 1/rgamma(1,shape=(nu.tau+k)/2,rate=(nu.tau*s.2.tau+sum((lambda-mu)^2))/2)
      sigma2 <- 1/rgamma(1,shape=(nu.sigma+k*n)/2,rate=(nu.sigma*s.2.sigma+sum((data-lambda)^2))/2)
      mu <- rnorm(1,mean=((m.mu/s2.mu)+(k*mean(lambda)/tau2))/((1/s2.mu)+(k/tau2)),sd=sqrt(1/((1/s2.mu)+(k/tau2))))
      lambda <- rnorm(k,mean=((mu/tau2)+(n*apply(data,1,mean)/sigma2))/((1/tau2)+(n/sigma2)),sd=sqrt(rep(1/((1/tau2)+(n/sigma2)),k)))
      thislambda <- rnorm(1,mean=mu,sd=sqrt(tau2))
      write(thislambda,paste("trace.lambda.p_",j,sep=""),append=T)
    }
  }
}

ohagan.data <- as.matrix(read.table("data.txt",sep=","))
M <- 10000
mcmc(ohagan.data,M,lambda.0=rep(1,5),mu.0=1,m.mu=2,s2.mu=10,nu.sigma=20,s.2.sigma=22/20,nu.tau=20,s.2.tau=6/20)
mcmc.c(ohagan.data,M,lambda.0=1,nu.sigma=20,s.2.sigma=22/20)
mcmc.p(ohagan.data,M,lambda.0=rep(1,4),mu.0=1,m.mu=2,s2.mu=10,nu.sigma=20,s.2.sigma=22/20,nu.tau=20,s.2.tau=6/20)

n.groups <- 5
lambda.i.c <- lambda.i.p <- lambda.i <- matrix(NA,M,n.groups)
for(j in 1:n.groups){
  lambda.i.c[,j] <- scan(paste("trace.lambda.c_",j,sep=""))
  lambda.i.p[,j] <- scan(paste("trace.lambda.p_",j,sep=""))
}
lambda.i <- as.matrix(read.table("trace.lambda.txt",sep=""))

source("iic-lcp-groups.R")
iic.lcp.groups(n.groups,M,lambda.i.c,lambda.i.p,lambda.i)

