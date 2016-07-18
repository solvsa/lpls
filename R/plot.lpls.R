plot.lpls <-
function(fit,comps=c(1,2),doplot=c(TRUE,TRUE,TRUE),level=c(2,2,2),
    arrow=c(1,0,1),xlim=c(-1,1),ylim=c(-1,1), samplecol=4, pathcol=2, varcol="grey70",
    X1varsize=1, X2varsize=1, X3varsize=1, sampleindex=1:dim(fit$corloadings$R22)[1], pathindex=1:dim(fit$corloadings$R3)[1], varindex=1:dim(fit$corloadings$R21)[1]){
    
#  require(car)
    plottype<-c("p","n")
    
    plot(xlim,ylim,"n",
    xlab=paste("Comp",comps[1]),
    ylab=paste("Comp",comps[2]),
    main=paste("Correlation loading plot from ",fit$options$type,"-LPLS analysis",sep=""))
    
    ellipse(c(0,0),matrix(c(1,0,0,1),2,2),radius=1,lty=1,lwd=1,col=1,center.pch=F)
    abline(h=0,lty=3)
    abline(v=0,lty=3)   
        
    if(doplot[2]){
        points(fit$corloadings$R21[varindex,comps[1]],fit$corloadings$R21[varindex,comps[2]],type=plottype[level[2]],pch=20,col=varcol,cex=X2varsize)
        points(fit$corloadings$R22[sampleindex,comps[1]],fit$corloadings$R22[sampleindex,comps[2]],type=plottype[level[2]],pch=20,col=samplecol,cex=2)
        if(level[2]==2){
           text(fit$corloadings$R21[varindex,comps[1]],fit$corloadings$R21[varindex,comps[2]],labels=colnames(fit$data$X2),cex=X2varsize,col=varcol)
           text(fit$corloadings$R22[sampleindex,comps[1]],fit$corloadings$R22[sampleindex,comps[2]],labels=rownames(fit$data$X2),cex=0.7,col=samplecol)
        }
        if(arrow[2]==1){
            arrows(0,0,fit$corloadings$R21[varindex,comps[1]],fit$corloadings$R21[varindex,comps[2]],lwd=2,col="grey70",length = 0.1)
            arrows(0,0,fit$corloadings$R22[sampleindex,comps[1]],fit$corloadings$R22[sampleindex,comps[2]],lwd=2,col=4,length = 0.1)    
        }
    }
    
    if(doplot[1]){
        points(fit$corloadings$R1[,comps[1]],fit$corloadings$R1[,comps[2]],type=plottype[level[1]],pch=20,col=3,cex=X1varsize)
        if(level[1]==2){
            text(fit$corloadings$R1[,comps[1]],fit$corloadings$R1[,comps[2]],labels=colnames(fit$data$X1),cex=X1varsize,col=3)    
        }
        if(arrow[1]==1){arrows(0,0,fit$corloadings$R1[,comps[1]],fit$corloadings$R1[,comps[2]],col=3,length = 0.1)}
    }
    
    if(doplot[3]){
        points(fit$corloadings$R3[pathindex,comps[1]],fit$corloadings$R3[pathindex,comps[2]],type=plottype[level[3]],pch=20,col=2,cex=X3varsize)
        if(level[3]==2){
            text(fit$corloadings$R3[pathindex,comps[1]],fit$corloadings$R3[pathindex,comps[2]],labels=colnames(fit$data$X3),cex=X3varsize,col=pathcol)    
        }
        if(arrow[3]==1){arrows(0,0,fit$corloadings$R3[pathindex,comps[1]],fit$corloadings$R3[pathindex,comps[2]],col=2,length = 0.1)}
    }
}
