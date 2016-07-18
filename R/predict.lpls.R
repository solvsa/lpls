predict.lpls <-
function(obj, X1new=NULL, X2new=NULL, X3new=NULL, exo.direction=c("X1","X3")){

#Prediction for endo-LPLS
if(obj$options$type=="endo"){

    if(any(is.na(X1new)) | any(is.na(X3new))) stop("Prediction requires complete predictor data.\n In case of cross-validation, use 'impute=TRUE' in model fit\n")
    k <- dim(X1new)[2]
    ntest <-dim(X1new)[1]
    l <- dim(X3new)[2]
    ptest <- dim(X3new)[1]
    
    #Centering and scaling of new observations
    if(obj$options$scaledata[1]){
        X1new <- scale(X1new,obj$means$mX1,attr(obj$data$X1,"scaled:scale"))
    }else{
        X1new <- scale(X1new,obj$means$mX1, scale=F)
    }

    if(obj$options$scaledata[3]){  
        X3new <- scale(X3new,obj$means$mX3, attr(obj$data$X3,"scaled:scale"))
    }else{
        X3new <- scale(X3new,obj$means$mX3, scale=F)
    }
    
    pred <- matrix(1,nrow=ntest,ncol=ptest)*obj$means$grandmX2 + X1new%*%obj$coefficients$C%*%t(X3new)
    if(!is.null(attr(obj$data$X2,"scaled:scale"))) cat("Warning: Only grand mean adjusted and column-scaled X2 has been predicted")
    
    res<-list(pred=pred)

}else if(obj$options$type=="exo"){
    n <- dim(obj$data$X2)[1]
    p <- dim(obj$data$X2)[2]    
    if(any(is.na(X2new))) stop("Prediction requires complete predictor data.\n In case of cross-validation, use 'impute=TRUE' in model fit\n")
    X2new <- as.matrix(X2new)
    ntest <- dim(X2new)[1]
    ptest <- dim(X2new)[2]
    rowm1 <- apply(X2new,1,mean)
    colm1 <- apply(X2new,2,mean)

        if(exo.direction=="X1"){   #Prediction of X1
        

                #Double centering of new obs
                if(!obj$options$doublecenter){
                    X2new <- X2new - matrix(obj$means$grandmX2,nrow=ntest,ncol=p)
                    }else{
                        X2new <- X2new-
                        t(matrix(rep(1,p),ncol=1)%*%rowm1) - 
                        matrix(rep(1,ntest),ncol=1)%*%obj$means$colmX2 + 
                        matrix(obj$means$grandmX2,nrow=ntest,ncol=p)    
                }
    
       
                if(is.null(attr(obj$data$X1,"scaled:scale"))){
                    pred<- matrix(1,ntest,1)%*%obj$means$mX1 + X2new%*%obj$coefficients$B1
                }else{
                    pred<- matrix(1,ntest,1)%*%obj$means$mX1 + X2new%*%obj$coefficients$B1*(matrix(1,ntest,1)%*%attr(obj$data$X1,"scaled:scale"))
                }
                
                res<-list(pred=pred)
                
 
        }
        
        if(exo.direction=="X3"){   #Prediction of X3
        

                #centering of new obs
                if(!obj$options$doublecenter){
                    X2new <- X2new - matrix(obj$means$grandmX2,nrow=n,ncol=ptest)
                }else{
                    X2new <- X2new-
                    t(matrix(rep(1,ptest),ncol=1)%*%obj$means$rowmX2) - 
                    matrix(rep(1,n),ncol=1)%*%colm1 + 
                    matrix(obj$means$grandmX2,nrow=n,ncol=ptest)    
                }

                if(is.null(attr(obj$data$X3,"scaled:scale"))){
                    pred <- matrix(1,ptest,1)%*%obj$means$mX3 + t(X2new)%*%obj$coefficients$B3
                }else{
                    pred <- matrix(1,ptest,1)%*%obj$means$mX3 + t(X2new)%*%obj$coefficients$B3*(matrix(1,ptest,1)%*%attr(obj$data$X3,"scaled:scale"))
                }                
                
                res<-list(pred=pred)
                
 
        }        

    }else if(obj$options$type=="exo_ort"){
        cat("No prediction method is implemented for method exo_ort\n")
        res <- NULL
    }
    res

}
