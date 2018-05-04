.checkKFInputs <- function(a0,P0,dt,ct,Tt,Zt,Qt,Ht,yt){
    n <- ncol(yt)
    d <- nrow(yt)
    m <- length(a0)

    ## -------------- check sizes ---------------------------------
    ## --------- dt
    if (length(dim(dt)) < 2){
        stopifnot(!is.null(dt))
        dt <- matrix(dt,nrow=m,ncol=n)
    }
    else if (length(dim(dt))==2){
        if (dim(dt)[2]==1) dt <- matrix(dt,nrow=m,ncol=n)
        else if (dim(dt)[2] > 1) stopifnot(dim(dt)[1]==m,dim(dt)[2]==n)
    }

    ## --------- ct
    if (length(dim(ct)) < 2){
        stopifnot(!is.null(ct))
        ct <- matrix(ct,nrow=d,ncol=n)
    }
    else if (length(dim(ct))==2){
        if (dim(ct)[2]==1) ct <- matrix(dt,nrow=d,ncol=n)
        else if (dim(ct)[2] > 1) stopifnot(dim(ct)[1]==d,dim(ct)[2]==n)
    }
    ## --------- Tt
    stopifnot(length(dim(Tt)) >= 2)
    if (length(dim(Tt))== 2){
        stopifnot(dim(Tt)[1]==m,dim(Tt)[2]==m)
        Tt <- array(Tt,dim=c(m,m,n))
    }
    else if (length(dim(Tt)) ==3){
        stopifnot(dim(Tt)[1]==m,dim(Tt)[2]==m,dim(Tt)[3]==n)
    }
    ## --------- Zt
    stopifnot(length(dim(Zt)) >= 2)
    if (length(dim(Zt))== 2){
        stopifnot(dim(Zt)[1]==d,dim(Zt)[2]==m)
        Zt <- array(Zt,dim=c(d,m,n))
    }
    else if (length(dim(Zt)) ==3){
        stopifnot(dim(Zt)[1]==d,dim(Zt)[2]==m,dim(Zt)[3]==n)
    }
    ## --------- Qt
    stopifnot(length(dim(Qt)) >= 2)
    if (length(dim(Qt))== 2){
        stopifnot(dim(Qt)[1]==m,dim(Qt)[2]==m)
        Qt <- array(Qt,dim=c(m,m,n))
    }
    else if (length(dim(Qt)) ==3){
        stopifnot(dim(Qt)[1]==m,dim(Qt)[2]==m,dim(Qt)[3]==n)
    }
    ## --------- Ht
    stopifnot(length(dim(Ht)) >= 2)
    if (length(dim(Ht))== 2){
        stopifnot(dim(Ht)[1]==d,dim(Ht)[2]==d)
        Ht <- array(Ht,dim=c(d,d,n))
    }
    else if (length(dim(Ht)) ==3){
        stopifnot(dim(Ht)[1]==d,dim(Ht)[2]==d,dim(Ht)[3]==n)
    }
    list(a0=a0,P0=P0,ct=ct,dt=dt,Zt=Zt,Tt=Tt,Qt=Qt,Ht=Ht)
}

