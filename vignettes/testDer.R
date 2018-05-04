require("numDeriv")
source("testKF.R")

Short <- 1

dP <- function(theta,y){
    lambda <- exp(theta)
    n <- length(y)
    res <- numeric(n)
    for (i in 1:n) res[i] <- dpois(y[i],lambda[i],T)
    sum(res)
}

jac.num <- function(theta,y) jacobian(dP,theta,y=y)

hess.num <- function(theta,y){
    hess <- hessian(dP,theta,y=y)
    zapsmall(hess)
}

jac_<-  function(theta,y,Pars=list()) (y-exp(theta))
hess_<- function(theta,y,Pars=list()) diag(-exp(theta))

ifelse(Short,load("GoalsData10Teams.RData"),load("10YearsGoals.RData"))

yt <- Goals
testNLSignalSmoothing <- function(yt=Goals,MaxErr=1e-6){
    d <- nrow(yt)
    n <- ncol(yt)
    m <- d
    #gp= matrix(0,d,n)
    
    phi <- 0.9975
    sig2 <- 0.000205
    P0 <- diag(sig2,m,m)
    a0 <- numeric(m)
    Zt <- diag(1,d,m)
    dt <- numeric(m)
    ct <- numeric(d)
    Tt <- diag(phi,m,m)
    Qt <- diag(sig2,m,m)

    PosteriorSignalMode(a0,P0,dt,ct,Tt,Zt,Qt,yt,
                        jac_,hess_,list(),
                        0,MaxErr,10,TRUE)
}

testSampling <- function(yt=Goals,seed=1200){
    d <- nrow(yt)
    n <- ncol(yt)
    m <- d
        
    phi <- 0.9975
    sig2 <- 0.205
    P0 <- diag(sig2,m,m)
    a0 <- numeric(m)
    Zt <- diag(1,d,m)
    dt <- numeric(m)
    ct <- numeric(d)
    Tt <- diag(phi,m,m)
    Qt <- diag(sig2,m,m)

    res <- PosteriorSignalSampling(a0,P0,dt,ct,Tt,Zt,Qt,yt,
                                   list(jac=jac_,hess=hess_),
                                   1,seed,"NLnonGaussian")
    list(theta=res$theta[,,1],
         lambda=exp(res$theta[,,1]),
         logLik=res$logLik)
}

testSampling.plot <- function(yt=Goals,Team="Chelsea",seed=1200){
    res <- testSampling(yt,seed)
    lambda <- res$lambda
    rownames(lambda) <- rownames(yt)
    n <- ncol(yt)
    plot(1:n,yt[Team,])
    points(1:n,lambda[Team,],col="red")
    leg <- c(paste("data_", Team),paste("lambda_", Team))
    legend("topleft",leg,col=c("black","red"),lwd=c(2,2),cex=0.7,bty="n")
}
