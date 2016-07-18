lplsCV <-
function(fit, segments1=NULL, segments2=NULL, trace=TRUE){
  X1 <- fit$data$X1
  X2 <- fit$data$X2
  X3 <- fit$data$X3
  n <- dim(fit$data$X2)[1]
  if(is.null(segments1) & is.null(segments2)) segments1 <- as.list(1:n)
  
  if(fit$options$type=="exo_ort") stop("No prediction method for type='exo_ort'. Use type='exo'.\n")
  if(fit$options$type=="exo"){
    if(!is.null(segments1)&!is.null(segments2)) stop("CV can only be run in one direction for type='exo'")
      if(!is.null(segments1)){
        nsegs <- length(segments1)
        
        rowpred <- dim(X1)[1]
        colpred <- dim(X1)[2]
        pred <- array(0,dim=c(rowpred, colpred,fit$ncomp))
        for(i in 1:nsegs){
          testX1 <- X1[segments1[[i]],,drop=F]
          testX2 <- X2[segments1[[i]],,drop=F]
          for(j in 1:fit$ncomp){
            trainfit <- lpls(X1, X2, X3, npc=j, 
                          doublecenter=fit$options$doublecenter,
                          scaledata=fit$options$scaledata,
                          type=fit$options$type,
                          subsetX1=-segments1[[i]])
            pred[segments1[[i]],,j] <- predict(trainfit,X2new=testX2, exo.direction="X1")$pred 
          }
          if(trace)cat(paste("Segment",i,"of",nsegs,"completed\n"))
        }
        dimnames(pred) <- list(dimnames(X1)[[1]],dimnames(X1)[[2]],paste("Comp",1:fit$ncomp,sep=""))
        rmsep<- apply(pred,3,function(x){sqrt(mean((x-X1)^2))})
      }else{
        nsegs <- length(segments2)
        
        rowpred <- dim(X3)[1]
        colpred <- dim(X3)[2]
        pred <- array(0,dim=c(rowpred, colpred,fit$ncomp))
                
        for(i in 1:nsegs){
          testX3 <- X3[segments2[[i]],,drop=F]
          testX2 <- X2[,segments2[[i]],drop=F]
                    
          for(j in 1:fit$ncomp){
            trainfit <- lpls(X1, X2, X3, npc=j, 
                             doublecenter=fit$options$doublecenter,
                             scaledata=fit$options$scaledata,
                             type=fit$options$type,
                             subsetX3=-segments2[[i]])
            pred[segments2[[i]],,j] <- predict(trainfit,X2new=testX2, exo.direction="X3")$pred
          }
          if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
        }

        dimnames(pred) <- list(dimnames(X3)[[1]],dimnames(X3)[[2]],paste("Comp",1:fit$ncomp,sep=""))
        rmsep<- apply(pred,3,function(x){sqrt(mean((x-X3)^2))})
        
      }
    
    }else if(fit$options$type=="endo"){
      if(!is.null(segments1)&!is.null(segments2)) stop("CV can only be run in one direction for type='exo'")
      if(!is.null(segments1)){
        nsegs <- length(segments1)
        
        rowpred <- dim(X2)[1]
        colpred <- dim(X2)[2]
        pred <- array(0,dim=c(rowpred, colpred,fit$ncomp))
        for(i in 1:nsegs){
          testX1 <- X1[segments1[[i]],,drop=F]
          testX2 <- X2[segments1[[i]],,drop=F]
          for(j in 1:fit$ncomp){
            trainfit <- lpls(X1, X2, X3, npc=j, 
                             doublecenter=fit$options$doublecenter,
                             scaledata=fit$options$scaledata,
                             type=fit$options$type,
                             subsetX1=-segments1[[i]])
            pred[segments1[[i]],,j] <- predict(trainfit,X1new=testX1,X3new=X3)$pred 
          }
          if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
        }
        dimnames(pred) <- list(dimnames(X2)[[1]],dimnames(X2)[[2]],paste("Comp",1:fit$ncomp,sep=""))
        rmsep<- apply(pred,3,function(x){sqrt(mean((x-X2)^2))})
      }else{
        nsegs <- length(segments2)
        
        rowpred <- dim(X2)[1]
        colpred <- dim(X2)[2]
        pred <- array(0,dim=c(rowpred, colpred,fit$ncomp))
        for(i in 1:nsegs){
          testX3 <- X3[segments2[[i]],,drop=F]
          testX2 <- X2[,segments2[[i]],drop=F]
          for(j in 1:fit$ncomp){
            trainfit <- lpls(X1, X2, X3, npc=j, 
                             doublecenter=fit$options$doublecenter,
                             scaledata=fit$options$scaledata,
                             type=fit$options$type,
                             subsetX3=-segments2[[i]])
            pred[,segments2[[i]],j] <- predict(trainfit,X1new=X1,X3new=testX3)$pred
          }
          if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
        }
        dimnames(pred) <- list(dimnames(X2)[[1]],dimnames(X2)[[2]],paste("Comp",1:fit$ncomp,sep=""))
        rmsep<- apply(pred,3,function(x){sqrt(mean((x-X2)^2))})
      }
  }
  names(rmsep)<-paste("Comp",1:fit$ncomp,sep="")
  return(list(rmsep=rmsep,pred=pred))
}
