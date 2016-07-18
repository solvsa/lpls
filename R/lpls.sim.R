lpls.sim <-
function(q=5,n=30,p=20,m=6,comps=2){
    #Data-simulation for LPLS testing (exo-type)

    #Simulations
    X2rowm<-rnorm(n,0,1)
    X2colm<-rnorm(p,2,2)
    X1colm<-rnorm(q,5,1)
    X3colm<-rnorm(p,0,1)
    X3rowm<-rnorm(m,2,1)

    X2<-matrix(rnorm(n*p,0,1),n,p)
    ULV<-svd(X2)
    T21<-ULV$u[,1:comps]%*%diag(sqrt(ULV$d[1:comps]))
    T22<-ULV$v[,1:comps]%*%diag(sqrt(ULV$d[1:comps]))
    P1<-matrix(rnorm(comps*q,4,1),q,comps)
    P3<-matrix(rnorm(comps*m,2,2),m,comps)
    X2<-T21%*%t(T22)
    X1<-T21%*%t(P1)
    X3<-T22%*%t(P3)

 #   X3<-t(X3)

    dimnames(X1)<-list(paste("n",1:n,sep=""),paste("class",1:q))
    dimnames(X2)<-list(paste("n",1:n,sep=""),as.character(1:p))
    dimnames(X3)<-list(as.character(1:p),paste("clust",1:m))

    return(list(X1=X1,X2=X2,X3=t(X3)))
}
