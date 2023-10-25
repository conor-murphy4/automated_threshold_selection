
source("thresh_qq_metric.R")
source("simulation_study/JointMLEFunctions.r")
library(threshr)
library(tea)

#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=1200, ncol=500)
for(ii in 1:n){
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- thresh_qq_metric(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]

  #Wadsworth 2016 method
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  if(!is.na(wadsthresh[ii])){
    fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
    optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
    wadsscale[ii] <- optwads$par[1]
    wadsshape[ii] <- optwads$par[2]
    wadslen[ii] <- length(fit.data_w)
    wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  }
  else{
    wadsscale[ii] <- NA
    wadsshape[ii] <- NA
    wadslen[ii] <- NA
    wadsquantile[ii] <- NA
  }

  #Northrop 2017 method
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
}


mythr_case1 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case1 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case1 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)


#Case 2: small sample size
n1 <- 500
Time_mythr <- Time_wads <- Time_wads_x <- Time_north <- Time_north_x <- numeric(n1)
mythresh <- wadsthresh <- norththresh <- numeric(n1)
myquantile <- wadsquantile <- northquantile <- numeric(n1)
myscale <- wadsscale <- northscale <- numeric(n1)
myshape <- wadsshape <- northshape <- numeric(n1)
mylen <- wadslen <- northlen <- numeric(n1)
data_matrix <- matrix(NA, nrow=480, ncol=500)
for(ii in 1:n1){
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- thresh_qq_metric(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num

  #Wadsworth 2016 method
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  if(!is.na(wadsthresh[ii])){
    fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
    optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
    wadsscale[ii] <- optwads$par[1]
    wadsshape[ii] <- optwads$par[2]
    wadslen[ii] <- length(fit.data_w)
    wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  }
  else{
    wadsscale[ii] <- NA
    wadsshape[ii] <- NA
    wadslen[ii] <- NA
    wadsquantile[ii] <- NA
  }

  #Northrop 2017 method
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
}


mythr_case2 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case2 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case2 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)


#Case 3: negative shape parameter
n2 <- 500
Time_mythr <- Time_wads <- Time_wads_x <- Time_north <- Time_north_x <- numeric(n2)
myquantile <- wadsquantile <- northquantile <- numeric(n2)
mythresh <- wadsthresh <- norththresh <- numeric(n2)
myscale <- wadsscale <- northscale <- numeric(n2)
myshape <- wadsshape <- northshape <- numeric(n2)
mylen <- wadslen <- northlen <- numeric(n2)
data_matrix <- matrix(NA, nrow=2400, ncol=500)
for(ii in 1:n2){
  set.seed(ii)
  dat1 <- runif(400, 0.5, 1.0)
  dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- thresh_qq_metric(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num

  #Wadsworth 2016 method
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]

  if(!is.na(wadsthresh[ii])){
    fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
    optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
    wadsscale[ii] <- optwads$par[1]
    wadsshape[ii] <- optwads$par[2]
    wadslen[ii] <- length(fit.data_w)
    wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  }
  else{
    wadsscale[ii] <- NA
    wadsshape[ii] <- NA
    wadslen[ii] <- NA
    wadsquantile[ii] <- NA
  }

  #Northrop 2017 method
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
}


mythr_case3 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case3 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case3 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)


# Case 4: partially-observed GPD
n3 <- 500
u <- 1.0
mythresh <- wadsthresh <- norththresh <- numeric(n3)
myquantile <- wadsquantile <- northquantile <- numeric(n3)
myscale <- wadsscale <- northscale <- numeric(n3)
myshape <- wadsshape <- northshape <- numeric(n3)
mylen <- wadslen <- northlen <- numeric(n3)
data_matrix <- matrix(NA, nrow=1000, ncol=500)
set.seed(12345)
for(ii in 1:n3){
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-u*rbeta(length(data_all),1,0.5)
  keep <- data_all>cens_thr
  data <- sample(data_all[keep], 1000, replace=FALSE)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- thresh_qq_metric(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num

  #Wadsworth 2016 method
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  if(!is.na(wadsthresh[ii])){
    fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
    optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
    wadsscale[ii] <- optwads$par[1]
    wadsshape[ii] <- optwads$par[2]
    wadslen[ii] <- length(fit.data_w)
    wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  }
  else{
    wadsscale[ii] <- NA
    wadsshape[ii] <- NA
    wadslen[ii] <- NA
    wadsquantile[ii] <- NA
  }
  
  #Northrop 2017 method
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
}


mythr_case4 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case4 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case4 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

