## These are the R scripts and numerical results accompanying 
## Krzysztof Bartoszek
## Cubic stochastic operators for Markov chains and mappings of distributions

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This script generates Figure 1. To recreate run in R  as
## set.seed(1)
## source("CSO_3D_Fig1.R")
## fplotlimitstraj()

f_VCSO_3Dex<-function(x,numprec=sqrt(.Machine$double.eps)){
    xm1<-x[1]
    x0<-x[2]
    x2<-x[3]


    Vm1<-xm1^3+3*xm1^2*(x0+x2)+2*xm1*x0*x2
    V0<-x0^3+3*x0^2*(xm1+x2)+4*xm1*x0*x2
    V2<-x2^3+3*x2^2*(xm1+x0)
    res<-c(Vm1,V0,V2)
    res[x==0]<-0
    res[x==1]<-1
    res[x==0.5]<-0.5
    if (!is.null(numprec)){
	res[res<numprec]<-0
	res[res>1]<-1	
	res[abs(res-0.5)<numprec]<-0.5
    }
    res<-res/sum(res)    
    res
}

r3d<-function(x2=NULL){
    res<-rep(NA,3)
    if (is.null(x2)){
        u<-sort(runif(2))
        res<-c(u,1)-c(0,u)
    }else{
        u<-runif(1,min=0,max=1-x2)
        res<-c(u,1-x2-u,x2)
    }
    res
}


f_rotateSimplexTriangle<-function(x){
    msqrt3d2<- (-sqrt(3)/2)
    rbind(c(-0.5,0.5),c(msqrt3d2,msqrt3d2))%*%x+c(0.5,-msqrt3d2)
}

fplotlimitstraj<-function(numplot=100000,numiter=20,filename="CSO_3D_basins",M=NULL,b_eqtria=TRUE){
    if(is.null(M)){M<-sapply(1:numplot,function(i){x<-r3d();x0<-x;for (i in 1:numiter){x<-f_VCSO_3Dex(x)};c(x0,x)},simplify=TRUE)}
    if (b_eqtria){M[1:2,]<-apply(M[1:2,],2,f_rotateSimplexTriangle)}
    vcolours<-c("#e87ba4","#1baf7a","#eb6834")
    png(paste0(filename,".png"));
    vcolM<-apply(M,2,function(x0x){switch(which.max(x0x[4:6]),vcolours[1],vcolours[2],vcolours[3])})
    plot(M[1,],M[2,],pch=19,col=vcolM,cex=0.5,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5,main="",xaxt="n",yaxt="n");
    legend("topright",pch=19,bty="n",col=vcolours,legend=c("(1,0,0)","(0,1,0)","(0,0,1)"),cex=1.5);
    dev.off()    
    invisible(NULL)
}


