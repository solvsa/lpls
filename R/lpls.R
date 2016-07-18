lpls <-
function(X1,X2,X3,npc=2,doublecenter=TRUE,scaledata=c(FALSE,FALSE,FALSE),
         type=c("exo"), impute=FALSE, niter=25, subsetX1=NULL, subsetX3=NULL,...){

  
#The three matrices are assumed positioned in this manner
#                 _________ 
#                |         |
#                |         |
#                |  t(X3)  |
#                |         |
#                |_________|
#                           
#                           
#   _______       _________ 
#  |       |     |         |
#  |       |     |         |
#  |  X1   |     |   X2    |
#  |       |     |         |
#  |_______|     |_________|
#                           

#  if(type!="exo_ort") impute=TRUE  
  if(impute){
    if(any(is.na(X1))) X1 <-.svd.imp(X1)
    if(any(is.na(X2))) X2 <-.svd.imp(X2)
    if(any(is.na(X3))) X3 <-.svd.imp(X3)
    }
  
    #Extracting subsets
    if(!is.null(subsetX1)){
        X1 <- X1[subsetX1,,drop=F]
        X2 <- X2[subsetX1,,drop=F]
        }

    if(!is.null(subsetX3)){
        X3 <- X3[subsetX3,,drop=F]
        X2 <- X2[,subsetX3,drop=F]
        }
  
    #For saving 
    X1save <- X1
    X2save <- X2
    X3save <- X3

    #Dimensions
    X1dim<-dim(X1)
    X2dim<-dim(X2)
    X3dim<-dim(X3)

    #Centering and scaling
    X1 <- scale(X1,scale=scaledata[1])
    mX1 <- attr(X1,"scaled:center")
    X3 <- scale(X3,scale=scaledata[3])
    mX3 <- attr(X3,"scaled:center")
 
    rowmX2 <- apply(X2,1,mean, na.rm=TRUE)  
    colmX2 <- apply(X2,2,mean, na.rm=TRUE)
    grandmX2 <- mean(X2, na.rm=TRUE)    
    
    #Do double centering of X2?
        if(doublecenter){
            X2 <- X2 - t(matrix(rep(1,X2dim[2]),ncol=1)%*%rowmX2) - 
            matrix(rep(1,X2dim[1]),ncol=1)%*%colmX2 + 
            matrix(grandmX2,nrow=X2dim[1],ncol=X2dim[2])  
        }else{
            X2 <- X2 - matrix(grandmX2,nrow=X2dim[1],ncol=X2dim[2])
        }
    #Column scaling of X2
    X2 <- scale(X2,scale=scaledata[2])
  
    attributes(X1save) <- attributes(X1)
    attributes(X2save) <- attributes(X2)
    attributes(X3save) <- attributes(X3)  

    dlist <- list(X1=as.matrix(t(X1)), X2=as.matrix(X2), X3=as.matrix(X3))


    T11 <- T12 <- T21 <- T22 <- T31 <- T32 <- numeric(0)
    W21 <- W22 <- numeric(0)
    P1  <- P3 <- P21 <- P22 <- numeric(0)
    D <- diag(rep(1,npc))
    X1totvar <- X2totvar <- X3totvar <- rep(0,npc+1)   

  for(a in 1:npc){
        
        X1totvar[a] <-sum(dlist$X1^2, na.rm=TRUE)
        X2totvar[a] <-sum(dlist$X2^2, na.rm=TRUE) 
        X3totvar[a] <-sum(dlist$X3^2, na.rm=TRUE)
        
        latent <- .extractscores(dlist, niter=niter)
                
        
        #Score vectors   
        T11<-cbind(T11, .lnorm(latent[[1]]$t1))
        T12<-cbind(T12, .lnorm(latent[[1]]$t2))        
        T21<-cbind(T21, .lnorm(latent[[2]]$t1))        
        T22<-cbind(T22, .lnorm(latent[[2]]$t2))        
        T31<-cbind(T31, .lnorm(latent[[3]]$t1))        
        T32<-cbind(T32, .lnorm(latent[[3]]$t2))         
        
        
        #Fitting an exo-LPLS model
        if(type=="exo_ort"){
        #Orthogonal scores
            P1 <- cbind(P1, .projectonto(dlist$X1,T21[,a]))
            P3 <- cbind(P3, .projectonto(dlist$X3,T22[,a]))
          
            P21 <- cbind(P21, .projectonto(dlist$X2,T21[,a]))
            P22 <- cbind(P22, .projectonto(dlist$X2,T22[,a]))
          
            d <- drop(solve(crossprod(T21[,a]))%*%t(T21[,a])%*%P22[,a])
            D[a,a] <- d
            
            #Deflation
            dlist$X1 <- dlist$X1 - P1[,a]%*%t(T21[,a])
            dlist$X3 <- dlist$X3 - T22[,a]%*%t(P3[,a])  
            dlist$X2 <- dlist$X2 - T21[,a]%*%t(P21[,a]) - P22[,a]%*%t(T22[,a]) + T21[,a]%*%t(T22[,a])*drop(d)
        }#end type=="exo"
        
        
        if(type=="exo"){
        #Simpler model for X2, but non-orthogonal scores
            P1 <- .projectonto(t(X1),T21)
            P3 <- .projectonto(t(X3),T22)
            
            P21 <- .projectonto(t(X2),T21)
            P22 <- .projectonto(X2,T22)
            
            D <- solve(crossprod(T21))%*%t(T21)%*%P22
            
            #Deflation
            dlist$X1 <- t(X1 - T21%*%t(P1))
            dlist$X3 <- X3 - T22%*%t(P3) 
            dlist$X2 <- X2 - T21%*%D%*%t(T22)
            
        }#end type=="exo"        
        
         
        if(type=="endo"){
            P1 <- .projectonto(t(X1),T12)
            P3 <- .projectonto(t(X3),T31)
            D <- solve(crossprod(T12))%*%(t(T12)%*%.projectonto(X2, T31))

            #Deflation
            dlist$X1 <- t(X1 - T12%*%t(P1))
            dlist$X3 <- X3 - T31%*%t(P3)                         
            dlist$X2 <- X2 - T12%*%(D%*%t(T31))
        
        }#end type=="endo"

    }
    
        X1totvar[npc+1] <-sum(dlist$X1^2, na.rm=TRUE)
        X2totvar[npc+1] <-sum(dlist$X2^2, na.rm=TRUE) 
        X3totvar[npc+1] <-sum(dlist$X3^2, na.rm=TRUE)
        
        X1varprop <- diff((X1totvar[1]-X1totvar)/(X1totvar[1]))
        X2varprop <- diff((X2totvar[1]-X2totvar)/(X2totvar[1]))
        X3varprop <- diff((X3totvar[1]-X3totvar)/(X3totvar[1]))
        
  

    #Various output
        B1 <- B3 <- Ca <- NULL
    if(type!="endo"){ 
        B1 <- T31%*%(solve(crossprod(P21,T31))%*%t(P1))
        B3 <- T12%*%(solve(crossprod(P22,T12))%*%t(P3))
        options(warn=-1)
        R1  <- cor(X1,T21, use="pairwise.complete.obs") 
        R21 <- t(cor(T21,X2, use="pairwise.complete.obs")) 
        R22 <- t(cor(T22,t(X2), use="pairwise.complete.obs"))  
        R3  <- cor(X3,T22, use="pairwise.complete.obs")  
        R2rmean <- cor(rowmX2,T21, use="pairwise.complete.obs")  
        R2cmean <- cor(colmX2,T22, use="pairwise.complete.obs")
        options(warn=0)
    }else if(type=="endo"){
        V1 <- T11%*%(solve(crossprod(P1,T11)))
        V3 <- T32%*%(solve(crossprod(P3,T32)))
        Ca <- V1%*%D%*%t(V3)
        options(warn=-1)
        R1  <- cor(X1,T12, use="pairwise.complete.obs") 
        R21 <- t(cor(T12,X2, use="pairwise.complete.obs")) 
        R22 <- t(cor(T31,t(X2), use="pairwise.complete.obs"))  
        R3  <- cor(X3,T31, use="pairwise.complete.obs")  
        R2rmean <- cor(rowmX2,T12, use="pairwise.complete.obs")  
        R2cmean <- cor(colmX2,T31, use="pairwise.complete.obs")
        options(warn=0)
    }
  
    res <- list(call=match.call())
    res$ncomp <- npc
    res$coefficients <- list(B1=B1, B3=B3, C=Ca)
    res$scores <- list(T11=T11, T12=T12, T21=T21, T22=T22, T31=T31, T32=T32)
    res$loadings <- list(P1=P1, P3=P3, P21=P21, P22=P22, D=D)
    res$corloadings <- list(R1=R1, R21=R21, R22=R22, R3=R3, R2rmean=R2rmean, R2cmean=R2cmean)
    res$means <- list(mX1=mX1, mX3=mX3, grandmX2=grandmX2, rowmX2=rowmX2, colmX2=colmX2)  
    res$data <- list(X1=X1save, X2=X2save, X3=X3save)
    res$residuals <- dlist
    res$options <- list(doublecenter=doublecenter, scaledata=scaledata, type=type)
    res$vars <- list(X1varprop=X1varprop, X2varprop=X2varprop, X3varprop=X3varprop)
    class(res)<-c("lpls")
    return(res)
}
