.extractscores <-
  function(datalist,niter=25, truncvec=rep(FALSE,2*length(datalist)), level=0.05){
    nmat<-length(datalist)
    dims<-lapply(datalist,dim)
    dimok <- TRUE
    if(length(level)==1) level <- rep(level, 2*length(datalist))
    
    for(ii in 1:(nmat-1)){
      dimok <- ifelse(dims[[ii]][2]==dims[[(ii+1)]][1],TRUE, FALSE)
    }
    if(!dimok) stop("Dimension mismatch\n")
    scorelist<-list()
    for(ii in 1:nmat){
      scorelist[[ii]]<-list(t1=matrix(0,nrow=dims[[ii]][1],ncol=1),
                            t2=matrix(0,nrow=dims[[ii]][2],ncol=1))
    }
    
    #Initiating the NIPALS algorithm
    scorelist[[1]]$t1<-as.matrix(rnorm(dims[[1]][1],0,1))

    #NIPALS
    for(k in 1:niter){
      t2<-.projectonto(datalist[[1]],scorelist[[1]]$t1)
      if(truncvec[2]) t2 <- .lenth.trunc(t2, alpha=level[2])$w
      scorelist[[1]]$t2 <- .lnorm(t2)
      
      for(j in 2:nmat){
        t2<-.projectonto(datalist[[j]],scorelist[[(j-1)]]$t2)
        if(truncvec[2*j]) t2 <- .lenth.trunc(t2, alpha=level[2*j])$w
        scorelist[[j]]$t2<-.lnorm(t2)
      }
      t1<-.projectonto(datalist[[nmat]],scorelist[[nmat]]$t2)
      if(truncvec[2*j-1]) t1 <- .lenth.trunc(t1, alpha=level[2*j-1])$w
      scorelist[[nmat]]$t1 <- .lnorm(t1)
      for(j in (nmat-1):1){
        t1<-.projectonto(datalist[[j]],scorelist[[(j+1)]]$t1)
        if(truncvec[2*j-1]) t1 <- .lenth.trunc(t1, alpha=level[2*j-1])$w
        scorelist[[j]]$t1 <- .lnorm(t1)
      }
    }
    return(scorelist)
  }
.lenth.trunc <- function(w, alpha){
  #This version will always select at least one variable
  s0<-1.5*median(abs(w))
  w0<-w[abs(w)<2.5*s0]
  sd<-1.5*median(abs(w0)) 
  upper<-qnorm((1-alpha/2),0,sd)
  lower<-qnorm(alpha/2,0,sd)
  lo<-sum(w>lower&w<upper)
  if(lo==0)lower <- upper <- abs(max(w)); lower <- -lower
  w[w>lower&w<upper]<-rep(0,lo)
  return(list(sd=sd,upper=upper,lower=lower,w=w))
}
.lnorm<-function(vec){
  vec/sqrt(drop(crossprod(vec)))
}
.projectonto <-
  function(A,b){
    A.na <- any(is.na(A))  
    if(is.null(dim(A))) A<-matrix(A,ncol=1)
    if(is.null(dim(b))) b <- matrix(b, ncol=1)
    
    if(!any(dim(A)==dim(b)[1])) stop("Non-matching dimensions")
    q <- dim(b)[2]

    #Orthogonalize b if multiple b's
    if(A.na && q>1){
      svd1 <- .svd.nipals(b, 10, q)
      b <- b%*%svd1$v
    }
  
    
    if(dim(A)[1]==dim(b)[1]) A <- t(A)
    if(!A.na){ 
      proj<-A%*%b%*%solve(crossprod(b))
    }else{
      miss <- which(is.na(A))
      A[miss]<-0    
      if(q==1){
        bb <- t(b^2)
        ones <- matrix(1,dim(A)[1],dim(A)[2])
        ones[miss] <- 0
        ones <- t(ones)
        btbinv <- 1/drop(bb%*%ones)
        proj <- (A%*%b)*btbinv      
      }else if(q>1){
        proj <- numeric()
        for(i in 1:q){
          bb <- t(b[,i]^2)
          ones <- matrix(1,dim(A)[1],dim(A)[2])
          ones[miss] <- 0
          ones <- t(ones)
          btbinv <- 1/drop(bb%*%ones)
          proj <- cbind(proj,(A%*%b[,i])*btbinv)      
        }
        proj <- proj%*%t(svd1$v)
      }
    }
    return(proj)
  }
.svd.imp <-
  function(X, max.niter=50, expl.min=0.98, interactive=FALSE, tol=1e-3, ploteval=FALSE){
    
    ncomp=min(dim(X))
    missing <- which(is.na(X))
    X.scaled <- scale(X,scale=F)
    colmeans <- attr(X.scaled,"scaled:center")
    meanmat <- matrix(1,nrow=dim(X)[1],ncol=1)%*%colmeans
    X.imp <- X
    impvals <- meanmat[missing]
    X.imp[missing] <- impvals
    initiate <- TRUE
    j<-1
    relchange <- 1
    change <- numeric()
    while(relchange > tol & j<=max.niter){        
      uvd <- .svd(X.imp)
      
      if(initiate){
        ssx <- rep(0,ncomp)
        for(i in 1:ncomp){
          D <- matrix(0,i,i)
          diag(D) <- uvd$d[1:i]
          ssx[i] <- sum((uvd$u[,1:i,drop=F]%*%D%*%t(uvd$v[,1:i]))^2)
          initiate <- FALSE
        }
        ssxtot <- sum(X.imp^2)
        explvar <- ssx/ssxtot
        if(interactive){
          plot(1:ncomp,explvar)
          ncomp <- readline("Choose number of components for imputation \n")
          ncomp <- as.numeric(ncomp)
        }else{
          ncomp <- min(which(explvar > expl.min))  
        }
      }
      D <- matrix(0,ncomp,ncomp)
      diag(D) <- uvd$d[1:ncomp]        
      Xhat <- uvd$u[,1:ncomp,drop=F]%*%D%*%t(uvd$v[,1:ncomp,drop=F])
      impvallength <- sqrt(sum((impvals)^2))
      change[j] <- relchange <- sqrt(sum((impvals-Xhat[missing])^2))
      impvals <- Xhat[missing]
      X.imp[missing] <- impvals
      if(ploteval){
        plot(0:j,c(1,change/impvallength),ylab="Change",main=paste("Relative change in imputed values using",ncomp,"components\n"),xlab="iteration")
      }
      j <- j+1
    }
    X.imp
  }
.svd.nipals <- function(X, niter, ncomp)
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    miss.x1 <- miss.x2 <- NULL
    x.na <- any(is.na(X))
    if(x.na)
    {
      miss.x1 <- which(is.na(X))    
      miss.x2 <- which(is.na(t(X)))
    }
    
    
    #Initiation
    w.x <- matrix(0,nrow=p,ncol=1)
    p.x <- matrix(0,nrow=p,ncol=1)
    
    u <- matrix(0,nrow=n,ncol=ncomp)
    d <- rep(0,ncomp)
    v <- matrix(0, nrow=p, ncol=ncomp)
    
    for(j in 1:ncomp){
      
      t.x <- .lnorm(rnorm(n,0,1))
      #NIPALS
      for(k in 1:niter)
      {
        w.x <- .lnorm(.projectonto(X,t.x))#,miss.x1))
        t.x <- .projectonto(t(X),w.x)##,miss.x2)
        x.scores <- t.x
        dj <- sqrt(drop(crossprod(t.x)))
        t.x <- .lnorm(t.x)
        p.x <- .projectonto(X,x.scores)##,miss.x1)
      }
      X <- X - x.scores%*%t(p.x)
      u[,j] <- t.x
      v[,j] <- w.x
      d[j] <- dj
    }
    list(u = u, d = d, v = v, resX = X)
  }
