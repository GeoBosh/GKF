## ==============================================================================================================
## -------------------------Posterior Signal mode computation ---------------------------------------------------
## ==============================================================================================================
##  Iterates until convergence NRUpdatingStep starting with a given initial guess
##  The function requires the following Kalman Filter output: a0,P0,dt,ct,Tt,Zt,Qt and yt
##  where here yt is the actual observation matrix at each time step.
##  In also requires routines to compute first (jac) and second (hess) derivatives of log p(yt|thetat)
##  at each time step t (this is due to the block-diagonal form of the different matrices)
##  We choose the following signatures for those function:
##  ***** jac=function(thetat,yt,Pars) where Pars is a list with the extra parameters (could be empty) 
##        and returns a (dx1) vector (numeric), d is length(yt)=length(thetat)
##  ***** hess=function(thetat,yt,Pars) similar to jac but returns a matrix (dxd)
##  We also require a tol parameter (default=1e-6) to declare convergence and max-number of iterations (default=10)

PosteriorSignalMode <- function(a0,P0,dt,ct,Tt,Zt,Qt,yt,jac,hess,Pars,
                                g=0,tol=1e-6,maxIter=10,checkInputs=TRUE){

    d <- nrow(yt)
    n <- ncol(yt)
    m <- length(a0)

    ## convert initial guess to appriopriate size
    if (g==0) gp <- matrix(0,d,n)
    else gp <- matrix(gp,d,n)

    iter <- 0
    err <- Inf

    while(iter < maxIter & err > tol){
        gpPrev <- gp
        gp <- .ComputeGpIter(d,n,a0,P0,dt,ct,Tt,Zt,Qt,
                             gpPrev,yt,jac,hess,Pars,
                             checkInputs)
        err <- .CompteGpErr(gp,gpPrev)
        iter <- iter+1
    }

    if (!is.null(rownames(yt))) rownames(gp) <- rownames(yt)
    list(theta_hat=gp,Iter=iter,err=err)
}

## ==============================================================================================================
## -------------------------Posterior Signal Sampling ---------------------------------------------------
## ==============================================================================================================
##  Wrapper for thetaSampling and GaussianthetaSampling allowing sampling from the posterior signal dis
##  The function requires the following Kalman Filter output: a0,P0,dt,ct,Tt,Zt,Qt and yt
##  where here yt is the actual observation matrix at each time step.
##  In also requires type specific parameters stored in the SpecificPars list
##  ********** NonLinearNonGaussian case:
##             routines to compute first (jac) and second (hess) derivatives of log p(yt|thetat)
##             at each time step t (this is due to the block-diagonal form of the different matrices)
##             as well as  tol parameter (default=1e-6) to declare convergence, initial guess g (default 0)
##             and max-number of iterations (default=10), see PosteriorSignalMode
## **********  LinearGaussian case:
##             Only Ht array (matrix) is required
## User can also choose the number of sample required 
## The function returns a list :
## ------- the sample matrices (an array where each slice corresponds to one sample)
## ------- the vector of (signal) cond log-likelihood l(theta|Yn).
## ------- the  (obsevation) Uncond log-likelihood l(Yn)
## ------- the vector of (signal) Uncond log-likelihood l(theta^i),i=1,..,M

PosteriorSignalSampling <- function(a0,P0,dt,ct,Tt,Zt,Qt,yt,
                                    SpecificPars,M=1,
                                    type=c("LinGaussian","NLnonGaussian")
                                    ){
    type <- match.arg(type)
    switch(type,
           LinGaussian={
               ## check we have all needed parameters
               stopifnot(!is.null(Ht <- SpecificPars$Ht),
                         is.numeric(Ht))
               
               resSignal <- GaussianthetaSampling(a0,P0,dt,ct,Tt,Zt,Ht,
                                                  Qt,yt,M,TRUE)
               thetaSample <- resSignal$theta
               thetaLogLik <- resSignal$logLik
               resKalmanFilter <- FKF(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,
                                      ComputeThetaLik=TRUE,
                                      ThetaVal=list(theta=thetaSample),
                                      checkInputs=TRUE)
               return(list(theta=thetaSample,
                           CondthetaLogLik=thetaLogLik,
                           UncondobsLogLik=resKalmanFilter$logLik,
                           UncondthetaLogLik=resKalmanFilter$thetalogLik)
                      )
           },
           NLnonGaussian={
               ## check we have all needed parameters
               stopifnot(!is.null(jac <- SpecificPars$jac),
                         !is.null(hess <- SpecificPars$hess)
                         )

               ## set not provided parameters to default
               if (is.null(SpecificPars$Pars)) Pars <- list()
               else Pars <- SpecificPars$Pars
               if (is.null(SpecificPars$g)) g <- 0
               else g <- SpecificPars$g
               if (is.null(SpecificPars$tol)) tol <- 1e-6
               else tol <- SpecificPars$tol
               if (is.null(SpecificPars$maxIter)) maxIter <- 10
               else maxIter <- SpecificPars$maxIter

               theta_hat <- PosteriorSignalMode(a0,P0,dt,ct,Tt,Zt,Qt,
                                                yt,jac,hess,Pars,
                                                g,tol,maxIter,TRUE)$theta_hat
               Ax <- .computeAx(nrow(yt),ncol(yt),theta_hat,yt,jac,hess,Pars)
               resSignal <- thetaSampling(a0,P0,dt,ct,Tt,Zt,Ax$A,Qt,Ax$x,
                                          theta_hat,M)
               thetaSample <- resSignal$theta
               thetaLogLik <- resSignal$logLik
               resKalmanFilter <- FKF(a0,P0,dt,ct,Tt,Zt,Ax$A,Qt,Ax$x,
                                      ComputeThetaLik=TRUE,
                                      ThetaVal=list(theta=thetaSample),
                                      checkInputs=TRUE)
               return(list(theta=thetaSample,
                           CondthetaLogLik=thetaLogLik,
                           UncondobsLogLik=resKalmanFilter$logLik,
                           UncondthetaLogLik=resKalmanFilter$thetalogLik)
                      )     
           })
}

## Mean squared error computation
.CompteGpErr <- function(gp,gpPrev) mean(colMeans((gp-gpPrev)**2))

## Compute one iteration of the Newton-Raphson algorithm
.ComputeGpIter <- function(d,n,a0,P0,dt,ct,Tt,Zt,Qt,
                           gpPrev,yt,jac,hess,Pars,
                           checkInputs=TRUE){
    
    Ax <- .computeAx(d,n,gpPrev,yt,jac,hess,Pars)
    NRUpdatingStep(a0,P0,dt,ct,Tt,Zt,Ax$A,
                   Qt,Ax$x,checkInputs)
}

## Compute A and x as defined in Jungbacker-Koopman(2007, Eq14)
.computeAx <- function(d,n,gpPrev,yt,jac,hess,Pars){
    A <- array(0,dim=c(d,d,n))
    x <- matrix(0,d,n)
    
    for (t in 1:n){
        th <- gpPrev[,t]
        y <- yt[,t]
        At <- -InvMat(hess(th,y,Pars))
        xt <- th + At%*%jac(th,y,Pars)
        A[,,t] <- At
        x[,t] <- xt
    }
    
    list(A=A,x=x)
}
