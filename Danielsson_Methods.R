source("_gpd.R")
#Danielsson et al 2019 
dan2019 <- function(data, k){
  n <- length(data)
  s_dat <- sort(data)
  max_dist <- numeric(length(k))
  T_max <- ceiling(0.2*n)
  j <- seq(1,T_max, by=1)
  for(i in 1:length(k)){
    gamma_k <- (1/k[i])*sum(log(s_dat[n:n-k[i]+1])-log(s_dat[n-k[i]]))  #Hill estimator
    q_est <- s_dat[n-k[i]+1]*((k[i]/j)^(gamma_k))
    dists <- abs(s_dat[n-j]-q_est)
    max_dist[i] <- max(dists) 
  }
  q_t <- min(max_dist)
  k_star <- k[which.min(max_dist)]
  u <- sort(data, decreasing = TRUE)[k_star]
  return(u)
}



#Danielsson 2001
# 
# install.packages("tea")
# library(tea)
# danielsson(data)

danielsson <-
  function(data,B=500,epsilon=0.9){
    n=length(data)
    n1=floor(n^epsilon)
    n2=(n1^2)/n

    Qn1=function (k) {
      xstat = sort(x1, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
      xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat-(2*xihat2^2)
    }
    Qn2=function (k) {
      xstat = sort(x2, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
      xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat-(2*xihat2^2)
    }

    qn1=matrix(nrow=B,ncol=n1-1)
    qn2=matrix(nrow=B,ncol=n2-1)
    for (l in 1:B){
      x1=sample(data,n1,replace=TRUE)
      x2=sample(data,n2,replace=TRUE)
      qn1[l,]=sapply(1:(n1-1),Qn1)
      qn2[l,]=sapply(1:(n2-1),Qn2)
    }
    qn1=qn1^2
    qn2=qn2^2
    qn1star=colMeans(qn1)
    qn2star=colMeans(qn2)
    k1star=which.min(qn1star)
    k2star=which.min(qn2star)

    Exp=(log(n1)-log(k1star))/(log(n1))
    Z=(log(k1star))^2
    N=(2*log(n1)-log(k1star))^2
    k0star=floor(k1star^2/k2star*((Z/N)^Exp))+1
    u=sort(data,decreasing=TRUE)[k0star]
    rho=log(k1star)/(-2*log(n1)+2*log(k1star))

    helphill=function (k) {
      xstat = sort(data, decreasing = TRUE)
      xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
      xihat
    }

    ti=1/helphill(k0star)
    list=list(sec.order.par=rho,k0=k0star,threshold=u,tail.index=ti)
    list
  }



#Sim study for Danielsson 2019

danthr <- danshape <- danscale <- danlen <- numeric(500)
set.seed(12345)
for(ii in 1:500){
  data <- rnorm(2000)
  k <- length(data)*seq(0, 0.95, by=0.05)
  danthr[ii] <- dan2019(data, k)
  fit.data <- data[data>danthr[ii]] - danthr[ii]
  optdan <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  danshape[ii] <- optdan$par[2]
  danscale[ii] <- optdan$par[1]
  danlen[ii] <- length(fit.data)
}


danthrG <- data.frame(thr=danthr, scale=danscale, shape=danshape, len=danlen)
saveRDS(danthrG, "SimGauss/Seeded_final/danthrG.rds")

library(evir)
data("nidd.thresh")
k <- length(nidd.thresh)*seq(0,0.95, by=0.05)
dan2019(nidd.thresh,k)


#-----------Danielsson 2001 results-------------------------------------

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')
hist(dan2001I$thr, breaks=30, main="Case 1")
hist(dan2001I1$thr, breaks=30, main="Case 2")
hist(dan2001I2$thr, breaks=30, main="Case 3")
hist(dan2001C$thr, breaks=30, main="Case 4")


#Danielsson 2019
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')
hist(danthrI$thr, breaks=30, main="Case 1")
hist(danthrI1$thr, breaks=30, main="Case 2")
hist(danthrI2$thr, breaks=30, main="Case 3")
hist(danthrC$thr, breaks=30, main="Case 4")

