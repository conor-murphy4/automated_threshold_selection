
source("src/eqd.R")
source("src/JointMLEFunctions.R")
library(threshr)

#Gaussian case
#n=2000
n <- 500
probs=seq(0.5, 0.95, by=0.05)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=2000, ncol=500)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data <- rnorm(2000)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
  
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


eqdthr_gauss <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_gauss <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_gauss <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Gaussian.csv")
write.csv(eqdthr_gauss, "output/threshold_selection/eqd_gauss.csv")
write.csv(wadsthr_gauss, "output/threshold_selection/wads_gauss.csv")
write.csv(norththr_gauss, "output/threshold_selection/north_gauss.csv")

#n=20000
n <- 500
probs=seq(0.5, 0.95, by=0.005)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=20000, ncol=500)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data <- rnorm(20000)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)

  #EQD method
  myres <- eqd(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess

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


eqdthr_gauss <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_gauss <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_gauss <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Gaussian_large.csv")
write.csv(eqdthr_gauss, "output/threshold_selection/eqd_gauss_large.csv")
write.csv(wadsthr_gauss, "output/threshold_selection/wads_gauss_large.csv")
write.csv(norththr_gauss, "output/threshold_selection/north_gauss_large.csv")