
source('Wadsworth 2016 source files/JointMLEFunctions.r')


library(ismev)


NHPP.diag<-function(x, u= NULL, k, q1, q2=1, par=NULL, M=NULL, nbs=1000, alpha=0.05, plot.out=TRUE, plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE,pmar=c(5.5,7,3,3),...)
{
  unull<-is.null(u)
  if(!unull){k<-length(u)}
  
  if(is.null(M))
  {
    if(unull){M<-length(x[x>quantile(x,q1)])/3}
    else{M<-length(x[x>min(u)])/3}
  }
  if(is.null(par))
  {
    require(ismev)
    {
      if(unull){ppf<-pp.fit(xdat=x,thresh=quantile(x,q1),npy=length(x)/M,show=F)}
      else{ppf<-pp.fit(xdat=x,thresh=min(u),npy=length(x)/M,show=F)}
      par<-ppf$mle}
  }
  
  #par(mfrow=c(length(plots), 1))
  par(mar=pmar)
  
  par(las=1)
  
  J1<-Joint.MLE.NHPP(x=x, u=u, k=k, q1=q1, q2=q2, par=par,M=M)
  warn<-any(eigen(J1$Cov.xi)$val<=.Machine$double.eps)
  if(!unull&&warn){
    ustar<-NA #added
    theta.hat <- rep(NA,3)
    return(list(thresh=ustar, mle.u=theta.hat)) #added
  }#stop("Estimated covariance matrix for xi not positive definite: try different thresholds")} #added
  else{ #added
    while(any(eigen(J1$Cov.xi)$val<=.Machine$double.eps))
    {
      k<-k-1
      J1<-Joint.MLE.NHPP(x=x,k=k, q1=q1, q2=q2, par=par,M=M)
    }
    if(warn){warning(paste("Estimated covariance matrix for xi not positive definite for initial k. Final k:", k))}
    
    if(unull){u<-quantile(x,seq(q1,q2,len=k+1))} 
    
    wn<-C1(k)%*%J1$mle[,3]/sqrt(diag(C1(k)%*%J1$Cov.xi%*%t(C1(k))))
    nl<-norm.LRT(x=wn,u=u[-c(1,k+1)])
    
    nlt<-NULL
    for(j in 1:nbs)
    {
      nlt[j]<-max(norm.LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
    }
    
    pval<-length(nlt[nlt>max(nl[,2])])/nbs
    
    if(pval<alpha){ustar<-unname(nl[nl[,2]==max(nl[,2]),1])} #added unname()
    else{ustar<-min(u)}
    ind<-u[-(k+1)]==ustar
    theta.hat<-J1$mle[ind,]
    
    if(unull){qs<- seq(q1, q2, len=k+1)[-(k+1)]}
    
    if(plot.out==TRUE){
      if(is.element("LRT",plots)){
        if(UseQuantiles && unull)
        {plot(qs, c(rep(NA,2),nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)}
        else{plot(u[-c(k+1)], c(rep(NA,2),nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)}
      }
      
      if(is.element("WN",plots)){
        if(UseQuantiles && unull){plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
          abline(h=0,col=2)
          abline(v=mean(x<=ustar),col=4)
        }
        else{plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
          abline(h=0,col=2)
          abline(v=ustar,col=4)}
      }
      
      if(is.element("PS",plots)){  
        TradCI<-cbind(J1$mle[,3]-qnorm(0.975)*sqrt(diag(J1$Cov.xi)),J1$mle[,3]+qnorm(0.975)*sqrt(diag(J1$Cov.xi)))
        if(UseQuantiles && unull)
        {plot(qs,J1$mle[,3],ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(xi)),...)
          lines(qs,TradCI[,1], lty=2)
          lines(qs,TradCI[,2], lty=2)
          abline(v=mean(x<=ustar),col=4)
        }
        else{plot(u[-(k+1)],J1$mle[,3],ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(xi)),...)
          lines(u[-(k+1)],TradCI[,1], lty=2)
          lines(u[-(k+1)],TradCI[,2], lty=2)
          abline(v=ustar,col=4)}
      }
    }
  }
  return(list(MLEall=J1$mle, Cov.xi=J1$Cov.xi, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}

# Expl.diag<-function(x, u= NULL, k, q1, q2=1, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE, param="InvRate",pmar=c(5.5,7,3,3),...)
# {
#   unull<-is.null(u)
#   if(!unull){k<-length(u)}
#   
#   par(mfrow=c(length(plots), 1))
#   par(mar=pmar)
#   
#   par(las=1)
#   
#   J1<-Joint.MLE.Expl(x=x,u=u,k=k, q1=q1, q2=q2, param=param)
#   warn<-any(eigen(J1$Cov)$val<=.Machine$double.eps)
#   if(!unull&&warn){stop("Estimated covariance matrix for eta not positive definite: try different thresholds")}
#   
#   while(any(eigen(J1$Cov)$val<=.Machine$double.eps))
#   {
#     k<-k-1
#     J1<-Joint.MLE.Expl(x=x,k=k, q1=q1, q2=q2, param=param)
#   }
#   if(warn){warning(paste("Estimated covariance matrix for 1/eta not positive definite for initial k. Final k:", k))}
#   
#   if(unull){u<-quantile(x,seq(q1,1,len=k+1))} 
#   
#   wn<-C1(k)%*%J1$mle/sqrt(diag(C1(k)%*%J1$Cov%*%t(C1(k))))
#   nl<-norm.LRT(x=wn,u=u[-c(1,k+1)])
#   
#   nlt<-NULL
#   for(j in 1:nbs)
#   {
#     nlt[j]<-max(norm.LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
#   }
#   
#   pval<-length(nlt[nlt>max(nl[,2])])/nbs
#   
#   if(pval<alpha){ustar<-nl[nl[,2]==max(nl[,2]),1]}
#   else{ustar<-min(u)}
#   ind<-u[-(k+1)]==ustar
#   theta.hat<-J1$mle[ind]
#   
#   if(unull){qs<- seq(q1, q2, len=k+1)[-(k+1)]}
#   
#   if(is.element("LRT",plots)){
#     if(unull&&UseQuantiles)
#     {plot(qs, c(NA,NA,nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)}
#     else{plot(u[-c(k+1)], c(NA,NA,nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)}
#   }
#   
#   if(is.element("WN",plots)){
#     if(unull&&UseQuantiles){plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
#       abline(h=0,col=2)
#       abline(v=mean(x<=ustar),col=4)
#     }
#     else{plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
#       abline(h=0,col=2)
#       abline(v=ustar,col=4)}
#   }
#   
#   if(is.element("PS",plots)){  
#     TradCI<-cbind(J1$mle-qnorm(0.975)*sqrt(diag(J1$Cov)),J1$mle+qnorm(0.975)*sqrt(diag(J1$Cov)))
#     if(UseQuantiles)
#     {
#       if(param=="InvRate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(eta)),...)}
#       else if(param=="Rate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(theta)),...)}
#       lines(qs,TradCI[,1], lty=2)
#       lines(qs,TradCI[,2], lty=2)
#       abline(v=mean(x<=ustar),col=4)
#     }
#     else{
#       if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(eta)),...)}
#       else if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(theta)),...)}
#       lines(u[-(k+1)],TradCI[,1], lty=2)
#       lines(u[-(k+1)],TradCI[,2], lty=2)
#       abline(v=ustar,col=4)
#     }
#   }
#   
#   return(list(MLE=J1$mle, Cov=J1$Cov, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
# }

#############################################################################################################
# NHPP examples
###############

# x <- rgpd(1000, shape=0.3, scale=0.4, mu=1.0)
# y <- runif(100,0.5,1)
# dat <- c(y,x)
# hist(dat)
# thresholds<- seq(0.8,2.3,by=0.05)
# thresholds <- quantile(probs=c(1:22)/23, dat)
# 
# NHPP.diag(dat,u=thresholds,UseQuantiles = FALSE, plot.out = FALSE)
# 
# 
# seeds <- seq(1,100,by=1)
# thrwads <- NULL
# thrgpd <- NULL
# threxp <- NULL
# for(j in 1:length(seeds)){
#   set.seed(seeds[j])
#   x <- rgpd(1000, shape=0.3, scale=0.4, mu=1.0)
#   y <- runif(100,0.5,1)
#   dat <- c(y,x)
#   thresholds<- seq(0.8,2.3,by=0.05)
#   wads <- NHPP.diag(dat, u=thresholds, UseQuantiles = FALSE, plot.out = FALSE)
#   reswads <- wads$thresh
#   resgpd <- distmetricfixedm(dat, thresholds = thresholds)[1]
#   resexp <- distmetricexpfixedm(dat, thresholds = thresholds)[1]
#   thrwads <- c(thrwads,reswads)
#   thrgpd <- c(thrgpd,resgpd)
#   threxp <- c(threxp,resexp)
# }
# 
# plot(thrwads[which(!is.na(thrwads))])
# 
# par(mfrow=c(1,1))
# plot(thrgpd, pch=19, col="blue")
# points(threxp,pch=0, col="red")
# points(which(!is.na(thrwads)),thrwads[which(!is.na(thrwads))], pch=20, col="black", ylim=c(0,2.3))
# points(which(is.na(thrwads)),rep(0.95, length(which(is.na(thrwads)))), pch=4, col="orange")
# library(Metrics)
# rmse(1,thrwads[which(!is.na(thrwads))])
# rmse(1,thrgpd)
# rmse(1,threxp)
# 
