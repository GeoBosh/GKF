require("RcppArmadillo")
require("Rcpp")
require("FKF")
require("dlm")

library("GKF")

Plot <- 0
Arima21 <- 1
LocalLevel <- 0
LinearGrowth <- 0
DynamicCAPM <- 0
    

if (Arima21){
    ## <--------------------------------------------------------------------------->
    ## Example 1: ARMA(2, 1) model estimation.
    ## <--------------------------------------------------------------------------->
    ## This example shows how to fit an ARMA(2, 1) model using this Kalman
    ## filter implementation (see also stats' makeARIMA and KalmanRun).
    n <- 1000
    ## Set the AR parameters
    ar1 <- 0.6
    ar2 <- 0.2
    ma1 <- -0.2
    sigma <- sqrt(0.2)
    
    ## Sample from an ARMA(2, 1) process
    a <- arima.sim(model = list(ar = c(ar1, ar2), ma = ma1), n = n,innov = rnorm(n) * sigma)
    
    ## Create a state space representation out of the four ARMA parameters
    arma21ss <- function(ar1, ar2, ma1, sigma) {
        Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
        Zt <- matrix(c(1, 0), ncol = 2)
        ct <- matrix(0)
        dt <- matrix(0, nrow = 2)
        GGt <- matrix(0)
        H <- matrix(c(1, ma1), nrow = 2) * sigma
        HHt <- H %*% t(H)
        a0 <- c(0, 0)
        P0 <- matrix(1e6, nrow = 2, ncol = 2)
        return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,HHt=HHt))
    }
    
    sp <- arma21ss(ar1,ar2,ma1,sigma)
    yt <- rbind(a)
    
    
fkfPack <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)

fkfCpp <- FKF(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
              Zt = sp$Zt, Qt = sp$HHt, Ht = sp$GGt, yt = yt)

}

if (LocalLevel){
    ## <--------------------------------------------------------------------------->
    ## Example 2: Local level model for the Nile's annual flow.
    ## <--------------------------------------------------------------------------->
    ## Transition equation:
    ## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)
    ## Measurement equation:
    ## y[t] = alpha[t] + eps[t], eps[t] ~ N(0, GGt)
    y <- c(Nile,Nile + runif(length(Nile),-1,1))
    y[c(3, 10)] <- NA # NA values can be handled
    ## Set constant parameters:
    dt <- ct <- matrix(0)
    Zt <- Tt <- matrix(1)

    a0 <- y[1] # Estimation of the first year flow
    P0 <- matrix(100) # Variance of 'a0'

    HHt <- matrix(var(y, na.rm = TRUE) * .5)
    GGt <- matrix(var(y, na.rm = TRUE) * .5)
    
    ## Filter Nile data with estimated parameters:
    fkfPack <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = HHt,GGt = GGt, yt = rbind(y))
    
    fkfCpp <- FKF_Rcpp(a0=a0, P0=P0, dt=dt, ct=ct, Tt=Tt, Zt=Zt,Qt=HHt,Ht=GGt,yt= rbind(y))
}

if (LinearGrowth){
    ## <--------------------------------------------------------------------------->
    ## Example 3: Multivariate Linear growth model (See Petris p128)
    ## <--------------------------------------------------------------------------->
    ## alpha_t=c(mu1_t,mu2_t,beta1_t,beta2_t)'
    ## Tt=[1,0,1,0|0,1,0,1|0,0,1,0|0,0,0,1]
    ## Rt=0 , Qt=diag(Q_mu,Q_beta), Q_mu=0,Q_beta=[49,155|155,437266]
    ## Zt=[1,0,0,0|0,1,0,0], Ht=[72,1018|1018,14353]
    ## Transition equation:
    ## alpha_t+1 = Tt*alpha_t + eta[t], eta[t] ~ N(0, Qt)
    ## Measurement equation:
    ## y_t = Zt*alpha_t + eps_t, eps_t ~ N(0, Ht)
    ## data used are annual groth for spain and Danemark collected in invest2.dat

    ## read data
    invest <- read.table("invest2.dat")
    ## Prepare model
    mod <- dlmModPoly(2)
    mod$FF <- mod$FF %x% diag(2)
    mod$GG <- mod$GG %x% diag(2)
    W1 <- matrix(c(0.5,0,0,0.5), 2, 2)
    W2 <- diag(c(49, 437266))
    W2[1, 2] <- W2[2, 1] <- 155
    mod$W <- bdiag(W1, W2)
    V <- diag(c(72, 14353))
    V[1, 2] <- V[2, 1] <- 1018
    mod$V <- V
    mod$m0 <- rep(0, 4)
    mod$C0 <- diag(4) * 1e4
    ## dlm computation
    filtered <- dlmFilter(invest, mod)
    logLikdlm <- dlmLL(as.matrix(invest), mod) ## up to a cte
    ## Smoothed values
    smoothed <- dlmSmooth(filtered) 
    alpha.Smoothed <- dropFirst(smoothed$s)
    theta.Smoothed <- t(mod$FF %*% t(alpha.Smoothed))
    ## Sampled values
    alpha.Sampled <- dlmBSample(filtered)
    theta.Sampled <- t(mod$FF %*% t(dropFirst(alpha.Sampled)))
    ## Rcpp
    l <- .checkKFInputs(a0=mod$m0,P0=mod$C0,dt=rep(0, 4),ct=rep(0, 2),
                        Tt=mod$GG,Zt=mod$FF,Qt=mod$W,Ht=mod$V,
                        yt= t(as.matrix(invest)))
    a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
    Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    yt=t(as.matrix(invest))
        
    ##fkfCpp <- FKF_Rcpp(a0= a0,P0=P0,dt=dt,ct=ct,Tt=Tt,Zt=Zt,Qt=Qt,
      ##                 Ht=Ht,yt=yt,checkInputs=FALSE)

    SmoothCpp <- GaussianSignalSmoothing(a0_=a0,P0_=P0,dt_=dt,ct_=ct,Tt_=Tt,
                                         Zt_=Zt,Ht_=Ht,Qt_=Qt,yt_=yt)
                                          
    SampleCpp <- GaussianthetaSampling(a0_=a0,P0_=P0,dt_=dt,ct_=ct,Tt_=Tt,
                                       Zt_=Zt,Ht_=Ht,Qt_=Qt,yt_=yt,
                                       M=1,seedVal=200)
        
    ##------------------- plot the result
    ## Extract Relevent part of the Data
    RcppSmooth <- t(SmoothCpp$theta_hat)
    RcppSample <- t(res$theta_tilda[,,1])
    n <- nrow(RcppSmooth)
    i=2
    
    if (Plot) pdf("dlmSamplingCompare_2.pdf")
    ylim=c(-50,max(c(theta.Smoothed[,i],theta.Sampled[,i],RcppSmooth[,i],RcppSample[,i]))+50)
    plot(1:n,theta.Smoothed[,i],type="l",col="black",ylim=ylim,lwd=2)
    points(1:n,theta.Sampled[,i],type="l",col="blue",ylim=ylim,lwd=2)
    points(1:n,RcppSmooth[,i],type="l",col="red",ylim=ylim,lwd=2)
    points(1:n,RcppSample[,i],type="l",col="green",ylim=ylim,lwd=2)

    leg <- c("dlm smoothed","dlm sampled","Rcpp smoothed", "Rcpp sampled")
    legend("topleft",leg,col=c("black","blue","red","green"),lwd=c(2,2,2,2),cex=0.7,bty="n")
    if (Plot) dev.off()
}

if (DynamicCAPM){
    ## <--------------------------------------------------------------------------->
    ## Example 4: Dynamic CAPM model (See Petris p132)
    ## <--------------------------------------------------------------------------->
    ## alpha_t=c(mu1_t,mu2_t,beta1_t,beta2_t)'
    ## Tt=[1,0,1,0|0,1,0,1|0,0,1,0|0,0,0,1]
    ## Rt=0 , Qt=diag(Q_mu,Q_beta), Q_mu=0,Q_beta=[49,155|155,437266]
    ## Zt=[1,0,0,0|0,1,0,0], Ht=[72,1018|1018,14353]
    ## Transition equation:
    ## alpha_t+1 = Tt*alpha_t + eta[t], eta[t] ~ N(0, Qt)
    ## Measurement equation:
    ## y_t = Zt*alpha_t + eps_t, eps_t ~ N(0, Ht)
    ## data used are annual groth for spain and Danemark collected in invest2.dat


    tmp <- ts(read.table("P.dat",header = TRUE),
              start = c(1978, 1), frequency = 12) * 100
    y <- tmp[, 1 : 4] - tmp[, "RKFREE"]
    colnames(y) <- colnames(tmp)[1 : 4]
    market <- tmp[, "MARKET"] - tmp[, "RKFREE"]
    rm("tmp")
    m <- NCOL(y)
    ## Set up the model
    CAPM <- dlmModReg(market)
    CAPM$FF <- CAPM$FF %x% diag(m)
    CAPM$GG <- CAPM$GG %x% diag(m)
    CAPM$JFF <- CAPM$JFF %x% diag(m)
    CAPM$W <- CAPM$W %x% matrix(0, m, m)
    CAPM$W[-(1 : m), -(1 : m)] <- c(8.153e-07, -3.172e-05, -4.267e-05,-6.649e-05,
                                    -3.172e-05, 0.001377, 0.001852, 0.002884,
                                    -4.267e-05, 0.001852, 0.002498, 0.003884,
                                    -6.649e-05, 0.002884, 0.003884, 0.006057)
    
    CAPM$V <- CAPM$V %x% matrix(0, m, m)
    CAPM$V[] <- c(41.06, 0.01571, -0.9504, -2.328,
                  0.01571, 24.23, 5.783, 3.376,
                  -0.9504, 5.783, 39.2, 8.145,
                  -2.328, 3.376, 8.145,39.29)
    
    CAPM$m0 <- rep(0, 2 * m)
    CAPM$C0 <- diag(1e7, nr = 2 * m)

    fkfdlm <- dlmFilter(y, CAPM)
    logLikdlm <- dlmLL(y,CAPM)
    ## Rcpp
    ## get the matrix Zt
    n <- nrow(y)
    Zt <- array(0,dim=c(m,2*m,n))
    Vmarket <- as.numeric(market)
    yt <- t(as.matrix(y))
    for (i in 1:n) Zt[,,i] <- cbind(1,Vmarket[i]) %x% diag(m)
    l <- .checkKFInputs(a0=CAPM$m0,P0=CAPM$C0,dt=rep(0,8),ct=rep(0, 4),
                        Tt=CAPM$GG,Zt=Zt,Qt=CAPM$W,Ht=CAPM$V ,
                        yt= yt)
    a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
    Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    
    fkfCpp <- FKF_Rcpp(a0= a0,P0= P0,dt=dt,ct=ct,Tt=Tt,Zt=Zt,Qt=Qt,
                       Ht=Ht,yt= yt,checkInputs=FALSE)
    UpdateCpp <- NRUpdatingStep(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt)
    
    fkfPack <- fkf(a0,P0,dt,ct,Tt,Zt,HHt=Qt,GGt=Ht,yt= yt)
}
