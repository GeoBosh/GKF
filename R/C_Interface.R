InvMat <- function(input) {
    .Call('_GKF_InvMat_Rcpp', PACKAGE = 'GKF', input)
}

nearPD <- function(x, w, corr = FALSE, keepDiag = FALSE, EnforcePosDef = TRUE,
                   doSym = FALSE, ensureSymmetry = FALSE,
                   eig_tol = 1e-6, conv_tol = 1e-7, posd_tol = 1e-3, maxit = 100,
                   chol = TRUE){

    if (missing(w)) w <- 1 + numeric(nrow(x))

    .Call('_GKF_nearPD_Rcpp', PACKAGE = 'GKF', x,w, corr, keepDiag,
          EnforcePosDef, doSym, ensureSymmetry, eig_tol,
          conv_tol, posd_tol, maxit, chol
          # 2018-03-28 removing this: , trace   - since there is no `trace' arg.
          #    apparently this was not picked up before the compulsory registration
         )
}

Robust_chol <- function(B, eig_tol = 0.001) {
    .Call('_GKF_Robust_chol_Rcpp', PACKAGE = 'GKF', B, eig_tol)
}

FKF <- function(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,ComputeThetaLik=FALSE,
                ThetaVal=list(),checkInputs=TRUE){
    if (checkInputs) {
        l <- .checkKFInputs(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt)
        if (ComputeThetaLik){
            stopifnot(!is.null(ThetaVal$theta),
                      all(dim(yt)==dim(ThetaVal$theta)[1:2]))
        }
        a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
        Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    }

    .Call('_GKF_FKF', PACKAGE = 'GKF',a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,
          ComputeThetaLik,ThetaVal)
}

NRUpdatingStep <- function(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,
                           checkInputs=TRUE){
    if (checkInputs) {
        l <- .checkKFInputs(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt)
        a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
        Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    }

    .Call('_GKF_NRUpdatingStep', PACKAGE = 'GKF',a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt)
}

thetaSampling <- function(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,thetaHat,
                          M,checkInputs=TRUE){
    if (checkInputs) {
        l <- .checkKFInputs(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt)
        a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
        Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    }

    .Call('_GKF_thetaSampling', PACKAGE = 'GKF',a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,
          thetaHat,M)
}

GaussianSignalSmoothing <- function(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,checkInputs=TRUE){
    if (checkInputs) {
        l <- .checkKFInputs(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt)
        a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
        Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    }

    .Call('_GKF_GaussianSignalSmoothing', PACKAGE = 'GKF',a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt)
}

GaussianthetaSampling <- function(a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,
                                  M,checkInputs=TRUE){
    if (checkInputs) {
        l <- .checkKFInputs(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt)
        a0 = l$a0 ; P0=l$P0 ; dt =l$dt ; ct=l$ct
        Tt= l$Tt  ; Zt=l$Zt ; Qt =l$Qt ; Ht=l$Ht
    }

    .Call('_GKF_GaussianthetaSampling', PACKAGE = 'GKF',a0,P0,dt,ct,Tt,Zt,Ht,Qt,yt,M)
}

## 2018-03-28 Commenting out since this is done in RcppExports.R
##
## # Register entry points for exported C++ functions
## methods::setLoadAction(function(ns) {
##     .Call('GKF_RcppExport_registerCCallable', PACKAGE = 'GKF')
## })
