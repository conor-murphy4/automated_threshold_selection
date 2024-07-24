#######################################################################################################
# Code taken from the supplementary materials of Wadsworth (2016) at https://www.tandfonline.com/doi/abs/10.1080/00401706.2014.998345 
# and adapted to count the number of samples where the method fails to estimate a threshold. 
#######################################################################################################

# NHPP.diag

# Details:

# Function to produce diagnostic plots and test statistics for the NHPP model

# Arguments:

# x - data vector
# u - optional; vector of candidate thresholds
# k - number of thresholds to consider (if u unspecified)
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost 
#       threshold will be at the (q2-1/k) quantile)
# M - number of superpositions or "blocks" / "years" the process corresponds to (can affect the optimization)
# nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
# alpha - significance level of the LRT
# plots - which plots to produce; "LRT"= likelihood ratio test, "WN" = white noise, "PS" = parameter stability
# UseQuantiles - logical; use quantiles as the thresholds in the plot?
# pmar - vector of length 4 giving the arguments for the plot margins in par(mar=c(*,*,*,*))

# Value: (results:)

# Plots the requested diagnostics. Also returns a list with components:

# MLEall - MLEs from all thresholds
# Cov.xi - Joint asymptotic covariance matrix for xi
# WN - values of the white noise process
# LRT - values of the LT test statistic vs threshold
# pval - p-value of the LR test
# k - final number of thresholds used
# thresh - threshold selected by LR procedure
# mle.u - MLE from selected threshold


NHPP.diag<-function(x, u= NULL, k, q1, q2=1, par=NULL, M=NULL, nbs=1000, alpha=0.05,plot.out=TRUE, plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE,pmar=c(5.5,7,3,3),...)
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

  par(mfrow=c(length(plots), 1))
  par(mar=pmar)

  par(las=1)
  J1<-Joint.MLE.NHPP(x=x, u=u, k=k, q1=q1, q2=q2, par=par,M=M)
  warn<-any(eigen(J1$Cov.xi)$val<=.Machine$double.eps)
  if(!unull&&warn){
    ustar<-NA #added by CM to count samples where method fails 
  return(list(thresh=ustar)) #added by CM
  }
  #stop("Estimated covariance matrix for xi not positive definite: try different thresholds")} #commented by CM
  
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
  
  if(pval<alpha){ustar<-nl[nl[,2]==max(nl[,2]),1]}
  else{ustar<-min(u)}
  ind<-u[-(k+1)]==ustar
  theta.hat<-J1$mle[ind,]

  if(unull){qs<- seq(q1, q2, len=k+1)[-(k+1)]}
  
  if(plot.out==TRUE){ #added by CM so plots are not outputted during simulation study
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
  
  return(list(MLEall=J1$mle, Cov.xi=J1$Cov.xi, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}


#############################################################################################################

# Expl.diag

# Details:

# Function to produce diagnostic plots and test statistics for the rate / inverse rate parameter of the Exponential model

# Arguments:

# x - data vector
# u - optional; vector of candidate thresholds
# k - number of thresholds to consider (if u unspecified)
# q1 -  lowest quantile for the threshold sequence
# q2 -  upper quantile limit for the threshold sequence (q2 itself is not used as a threshold, but rather the uppermost 
#       threshold will be at the (q2-1/k) quantile)
# nbs - number of simulations used to assess the null distribution of the LRT, and produce the p-value
# alpha - significance level of the LRT
# plots - which plots to produce; "LRT"= likelihood ratio test, "WN" = white noise, "PS" = parameter stability
# UseQuantiles - logical; use quantiles as the thresholds in the plot?
# param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively



# Value: (results:)

# Plots the requested diagnostics. Also returns a list with components:

# MLE - MLEs from all thresholds
# Cov - Joint asymptotic covariance matrix for xi
# WN - values of the white noise process
# LRT - values of the LT test statistic vs threshold
# pval - p-value of the LR test
# k - final number of thresholds used
# thresh - threshold selected by the LR procedure
# mle.u - MLE from selected threshold

Expl.diag<-function(x, u= NULL, k, q1, q2=1, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE, param="InvRate",pmar=c(5.5,7,3,3),...)
{
  unull<-is.null(u)
  if(!unull){k<-length(u)}
  
  par(mfrow=c(length(plots), 1))
  par(mar=pmar)
  
  par(las=1)
  
  J1<-Joint.MLE.Expl(x=x,u=u,k=k, q1=q1, q2=q2, param=param)
  warn<-any(eigen(J1$Cov)$val<=.Machine$double.eps)
  if(!unull&&warn){stop("Estimated covariance matrix for eta not positive definite: try different thresholds")}
  
  while(any(eigen(J1$Cov)$val<=.Machine$double.eps))
  {
    k<-k-1
    J1<-Joint.MLE.Expl(x=x,k=k, q1=q1, q2=q2, param=param)
  }
  if(warn){warning(paste("Estimated covariance matrix for 1/eta not positive definite for initial k. Final k:", k))}
  
  if(unull){u<-quantile(x,seq(q1,1,len=k+1))} 
  
  wn<-C1(k)%*%J1$mle/sqrt(diag(C1(k)%*%J1$Cov%*%t(C1(k))))
  nl<-norm.LRT(x=wn,u=u[-c(1,k+1)])
  
  nlt<-NULL
  for(j in 1:nbs)
  {
    nlt[j]<-max(norm.LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
  }
  
  pval<-length(nlt[nlt>max(nl[,2])])/nbs
  
  if(pval<alpha){ustar<-nl[nl[,2]==max(nl[,2]),1]}
  else{ustar<-min(u)}
  ind<-u[-(k+1)]==ustar
  theta.hat<-J1$mle[ind]
  
  if(unull){qs<- seq(q1, q2, len=k+1)[-(k+1)]}
  
  if(is.element("LRT",plots)){
    if(unull&&UseQuantiles)
    {plot(qs, c(NA,NA,nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)}
    else{plot(u[-c(k+1)], c(NA,NA,nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)}
  }
  
  if(is.element("WN",plots)){
    if(unull&&UseQuantiles){plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
                     abline(h=0,col=2)
                     abline(v=mean(x<=ustar),col=4)
    }
    else{plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
         abline(h=0,col=2)
         abline(v=ustar,col=4)}
  }
  
  if(is.element("PS",plots)){  
    TradCI<-cbind(J1$mle-qnorm(0.975)*sqrt(diag(J1$Cov)),J1$mle+qnorm(0.975)*sqrt(diag(J1$Cov)))
    if(UseQuantiles)
    {
     if(param=="InvRate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(eta)),...)}
     else if(param=="Rate"){plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(theta)),...)}
     lines(qs,TradCI[,1], lty=2)
     lines(qs,TradCI[,2], lty=2)
     abline(v=mean(x<=ustar),col=4)
    }
    else{
      if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(eta)),...)}
      else if(param=="InvRate"){plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold", ylab=expression(hat(theta)),...)}
      lines(u[-(k+1)],TradCI[,1], lty=2)
      lines(u[-(k+1)],TradCI[,2], lty=2)
      abline(v=ustar,col=4)
    }
  }
  
  return(list(MLE=J1$mle, Cov=J1$Cov, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}







#######################################################################################################

# Joint.MLE.Expl
##################

# Details: 

# Calculates the MLEs of the rate parameter, and joint asymptotic covariance matrix of these MLEs 
# over a range of thresholds as supplied by the user.

# Arguments:

# x - vector of data
# u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
# k - number of thresholds to consider if u not supplied
# q1, q2 - lower and upper quantiles to consider for threshold
# param - character specifying "InvRate" or "Rate" for either inverse rate parameter / rate parameter, respectively

# Value: (returns:)

# mle - vector of MLEs above the supplied thresholds
# cov - joint asymptotic covariance matrix of these MLEs


Joint.MLE.Expl<-function(x, u=NULL, k, q1, q2=1, param)
{ 
  if(!is.element(param, c("InvRate","Rate"))){stop("param should be one of InvRate or Rate")}
  if(!is.null(u))
  {
    k<-length(u)
    x<-x[x>u[1]]
    # add threshold above all data, to "close" final region
    u<-c(u,max(x)+1)
  }
  else{u<-quantile(x,seq(q1,q2,len=k+1))}
  
  I<-n<-m<-thetahat<-NULL
  
  for(i in 1:k)
  {
    if(param=="InvRate")
    {
    thetahat[i]<-mean(x[x>=u[i]]- u[i])
    }
    else if(param=="Rate")
    {    
    thetahat[i]<-1/mean(x[x>=u[i]]- u[i])
    }
    n[i]<-length(x[x>=u[i]&x<=u[i+1]])
  }
  for(i in 1:k)
  {
    m[i]<-sum(n[i:k])
    I[i]<- 1/thetahat[i]^2
  }
  
  Tcov<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in 1:k)
    {
      Tcov[i,j]<-1/(I[min(i,j)]*m[min(i,j)])
    }
  }
  CovT<-Tcov
  return(list(mle=thetahat, Cov=CovT))
}

#####################################################################################

# nhpp.nll 
########################

# Details:

# Negative log likelihood for the NHPP model, to minimize for MLEs

# Arguments:

# theta - parameter vector (mu, sigma, xi)
# x - data vector
# u - threshold
# M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
# but these can be changed post-hoc to correspond to any number)


# Value: (results:)

# negative log-likelihood value

nhpp.nll<-function(theta,x,u,M)
{
  x<-x[x>u]
  mu<-theta[1]; sig<-theta[2]; xi<-theta[3]
  Nu<-length(x)
  
  if(sig<=0|| any(1+xi*(x-mu)/sig < 0)){return(10e10)}
  
  else{
    if(abs(xi)>1e-10)
    {  
    nll<- (1/xi+1)*sum(log(1+xi*(x-mu)/sig)) + Nu*log(sig) + M*(1+xi*(u-mu)/sig)^(-1/xi)
    }
    else{
      nll<-Nu*log(sig)+sum((x-mu)/sig) + M*exp(-(u-mu)/sig)
    }
    return(nll)
  }
}



#####################################################################################


# Joint.MLE.NHPP
##################

# Details: 

# Calculates the MLEs of the parameters (mu,sigma,xi), and joint asymptotic covariance matrix of these MLEs 
# over a range of thresholds as supplied by the user.

# Arguments:

# x - vector of data
# u - vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
# k - number of thresholds to consider if u not supplied
# q1, q2 - lower and upper quantiles to consider for threshold
# par - starting values for the optimization
# M - number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of mu, sigma,
# but these can be changed post-hoc to correspond to any number)

# Value: (returns:)

# mle - matrix of MLEs above the supplied thresholds; columns are (mu, sigma, xi)
# Cov.all - joint asymptotic covariance matrix of all MLEs
# Cov.mu - joint asymptotic covariance matrix of MLEs for mu
# Cov.sig - joint asymptotic covariance matrix of MLEs for sig
# Cov.xi - joint asymptotic covariance matrix of MLEs for xi


Joint.MLE.NHPP<-function(x, u=NULL, k, q1, q2=1, par, M)
{
  if(!is.null(u))
  {
    k<-length(u)
    x<-x[x>u[1]]    
    # add threshold above all data, to "close" final region
    u<-c(u,max(x)+1)
  }
  else{u<-quantile(x,seq(q1,q2,len=k+1))}
  
  I<-Iinv<-list()
  thetahat<-matrix(NA,ncol=3,nrow=k)
  
  for(i in 1:k)
  {
    opt<-optim(nhpp.nll, par=par, x=x, u=u[i],M=M, hessian=F)
    thetahat[i,]<-opt$par
    
    ###
    # Deal with xi<-0.5
    if(thetahat[i,3]>-0.5)
    ###
    {
    I[[i]]<- E.Info.Mat(theta=opt$par, u=u[i], M=M)$EIM
    Iinv[[i]]<-solve(I[[i]])
    }
    else{I[[i]]<-Iinv[[i]]<-matrix(0,3,3)}
    }
    
  Wcov<-list()
  Wcov1<-NULL
  for(i in 1:k)
  {
    Wcov[[i]]<-matrix(0,3,3)
    for(j in 1:k)
    {
      Wcov[[i]]<-cbind(Wcov[[i]],Iinv[[min(i,j)]])
    }
    Wcov1<-rbind(Wcov1,Wcov[[i]])
  }
  Wcov1<-Wcov1[,-c(1:3)]
  
  CovT<-Wcov1
  
  Cov.mu<-CovT[seq(1,3*k,by=3),seq(1,3*k,by=3)]
  Cov.sig<-CovT[seq(2,3*k,by=3),seq(2,3*k,by=3)]
  Cov.xi<-CovT[seq(3,3*k,by=3),seq(3,3*k,by=3)]
  
  return(list(mle=thetahat, Cov.all=CovT, Cov.mu=Cov.mu, Cov.sig = Cov.sig, Cov.xi=Cov.xi))
}


###################################################################################

# norm.LRT

# Details:

# Evaluates the likelihood ratio statistics for testing white noise

# Arguments:

# x - vector of white noise process (WNP, usually normalized estimates of xi or the exponential rate paramter 1/eta)
# u - vector of thresholds that are associated to the WNP


norm.LRT<-function(x,u)
{
  l<-length(u)
  v<-u[-c(1)] # means two or more obs available for std dev calculation
  
  lr<-NULL
  for(i in 1:length(v))
  {
    n1<-length(x[u<=v[i]])
    num<-nll.norm(theta=c(mean(x[u<=v[i]]),sd(x[u<=v[i]])*sqrt((n1-1)/n1)),x=x[u<=v[i]])
    den<-nll.norm(theta=c(0,1),x=x[u<=v[i]])
    lr[i]<--2*(num-den)
  }
  return(cbind(v,lr))  
}



###################################################################################

# nll.norm - negative log liklihood for the normal distribution

nll.norm<-function(theta,x){
  if(theta[2]<0){return(10e10)}
  else{
    return(-sum(dnorm(x,mean=theta[1],sd=theta[2],log=T)))
  }
}



###################################################################################

# C1 

# Details: 

# Produces "Contrast matrix" with (1,-1) elements running down the two diagonals

# Arguments:

# k - number of columns (k-1 = number of rows)

# Value: (returns:)

# (k-1) x k contrast matrix

C1<-function(k)
{
  C<-matrix(0,k-1,k)
  for(i in 1:(k-1))
  {
    C[i,i]<-1
    C[i,i+1]<--1
  }
  return(C)
}


###################################################################################

# Functions for the expected information matrix used in Joint.MLE.NHPP

# Details:

# Functions with names of form "d2ldmu2" are derivatives (with expectation over the random NUMBER of points
# already incorporated). Functions with names of form "i.d2ldmu2", are in a form ready for integration; this integration
# yields the expectation over the random LOCATIONS of the points.

# Arguments:

# x - dummy variable over which to integrate
# mu, sig, xi - location, scale, shape parameters
# u - threshold above which NHPP model assumed
# M - number of superpositions or "blocks" / "years" the process corresponds to


d2ldmu2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<-M*c1*(xi*(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-2)
  p2<- -M*((1+xi)/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2)
}

i.d2ldmu2<-function(x, mu, sig, xi, u, M)
{
  d2ldmu2(x=x,mu=mu,sig=sig,xi=xi,u=u,M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



d2ldmudsig<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<-M*c1*(xi*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p2<-M*c1*(-(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p3<- M*(1/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4<- -M*((1+xi)/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2+p3+p4)
}

i.d2ldmudsig<-function(x, mu, sig, xi, u, M)
{
  d2ldmudsig(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

#########

d2ldmudxi<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<-M*c1* (1/sig)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(-(1+xi)/sig^2)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1<-(1/sig)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) + (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3<-M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  
  return(p1+p2+p3)
}

i.d2ldmudxi<-function(x, mu, sig, xi, u, M)
{
  d2ldmudxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



##########


d2ldsig2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<-M*c1*(-2*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(xi*(1+xi)/sig^4)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3<- M*(2/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4<- -M*((1+xi)/sig^4)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  p5<-M*c1*1/sig^2
  return(p1+p2+p3+p4+p5)
}

i.d2ldsig2<-function(x, mu, sig, xi, u, M)
{
  d2ldsig2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



#########

d2ldsigdxi<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<- M*c1*((x-mu)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p2<-M*c1*(-(1+xi)/sig^3)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1<-((u-mu)/sig^2)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) + (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3<-M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  
  return(p1+p2+p3)
}


i.d2ldsigdxi<-function(x, mu, sig, xi, u, M)
{
  d2ldsigdxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}


#########

d2ldxi2<-function(x, mu, sig, xi, u, M)
{
  c1<-(1+(xi/sig)*(u-mu))^(-1/xi)
  
  p1<- M*c1*(-2/xi^3)*log(1+(xi/sig)*(x-mu))#
  p2<- M*c1*2*((x-mu)/(sig*xi^2))*(1+(xi/sig)*(x-mu))^(-1)#
  p3<-M*c1*((1/xi+1)/sig^2)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)#
  
  p4.1<-((-1/xi^2)*log(1+(xi/sig)* (u-mu)) + (1/sig)*(1/xi)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))#
  p4<- -(p4.1^2)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#
  
  p5.1<-((2/xi^3)*log(1+(xi/sig)*(u-mu)) - (2/sig)*(1/xi^2)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1)- (1/sig^2)*(1/xi)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-2))#
  p5<- (p5.1)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#
  
  return(p1+p2+p3+p4+p5)
}

i.d2ldxi2<-function(x, mu, sig, xi, u, M)
{
  d2ldxi2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

## version of d2ldxi2 and i.d2ldxi2 with limit as xi \to 0. (This derivative has most trouble with small absolute
## values of xi.) Also possible to take such limits in other derivatives, but not implemented as they are generally 
## less problematic.


d2ldxi2.xi0<-function(x, mu, sig, xi, u, M)
{
  c1<-exp(-(u-mu)/sig)
  q1<-((x-mu)/sig)
  q2<-((u-mu)/sig)
  
  p1<--((2/3)*q1^3-q1^2)*M*c1
  p2<--(-(2/3)*q2^3+(1/4)*q2^4)*M*c1
  
  return(p1+p2)
}

i.d2ldxi2.xi0<-function(x, mu, sig, xi, u, M)
{
  d2ldxi2.xi0(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1/sig)*exp(-(x-u)/sig)
}



#############################################################################################

# E.Info.Mat

# Details:

# Calculates the numerically-integrated expected information matrix for an NHPP with specified parameters

# Arguments:

# theta - vector of parameters (mu,sigma, xi)
# u - threshold for NHPP
# M - number of superpositions or "blocks" / "years" the process corresponds to

# Value: (returns:)

# EIM - expected information matrix
# Errors - vector of errors from the numerical integration of the 6 unique components


E.Info.Mat<-function(theta, u, M)
{
  if(theta[3]<0)
  {
    up<-theta[1] - theta[2]/theta[3]
  }
  else{up<-Inf}
    
  Ed2ldmu2<-integrate(i.d2ldmu2,low=u,upper=up,mu=theta[1], sig=theta[2],
                      xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudsig<-integrate(i.d2ldmudsig,low=u,upper=up,mu=theta[1], sig=theta[2],
                         xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudxi<-integrate(i.d2ldmudxi,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsig2<-integrate(i.d2ldsig2,low=u,upper=up,mu=theta[1], sig=theta[2],
                       xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsigdxi<-integrate(i.d2ldsigdxi,low=u,upper=up,mu=theta[1], sig=theta[2],
                         xi=theta[3], u=u,M=M, abs.tol=0)
  if(abs(theta[3])>0.0001)
  {
    Ed2ldxi2<-integrate(i.d2ldxi2,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  }
  else{
    Ed2ldxi2<-integrate(i.d2ldxi2.xi0,low=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)  
  }
  
  Errors<-c(Ed2ldmu2$abs.err,Ed2ldmudsig$abs.err,Ed2ldmudxi$abs.err,Ed2ldmudsig$abs.err,Ed2ldsig2$abs.err,
            Ed2ldsigdxi$abs.err,Ed2ldmudxi$abs.err,Ed2ldsigdxi$abs.err,Ed2ldxi2$abs.err)
  
  EIM<--matrix(c(Ed2ldmu2$val,Ed2ldmudsig$val,Ed2ldmudxi$val,Ed2ldmudsig$val,Ed2ldsig2$val,Ed2ldsigdxi$val,
                 Ed2ldmudxi$val,Ed2ldsigdxi$val,Ed2ldxi2$val), byrow=T, nrow=3)
  
  return(list(EIM=EIM,Errors=Errors))
}