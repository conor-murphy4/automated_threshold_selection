
source("src/eqd.R")
source("src/JointMLEFunctions.R")
library(threshr)

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
  print(ii)
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
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


eqdthr_case1 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case1 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case1 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Case_1.csv")
write.csv(eqdthr_case1, "output/threshold_selection/eqd_case1.csv")
write.csv(wadsthr_case1, "output/threshold_selection/wads_case1.csv")
write.csv(norththr_case1, "output/threshold_selection/north_case1.csv")

#Case 2: small sample size
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=480, ncol=500)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
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


eqdthr_case2 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case2 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case2 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Case_2.csv")
write.csv(eqdthr_case2, "output/threshold_selection/eqd_case2.csv")
write.csv(wadsthr_case2, "output/threshold_selection/wads_case2.csv")
write.csv(norththr_case2, "output/threshold_selection/north_case2.csv")

#Case 3: negative shape parameter
myquantile <- wadsquantile <- northquantile <- numeric(n)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=2400, ncol=500)
for(ii in 1:n2){
  print(ii)
  set.seed(ii)
  dat1 <- runif(400, 0.5, 1.0)
  dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
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


eqdthr_case3 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case3 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case3 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Case_3.csv")
write.csv(eqdthr_case3, "output/threshold_selection/eqd_case3.csv")
write.csv(wadsthr_case3, "output/threshold_selection/wads_case3.csv")
write.csv(norththr_case3, "output/threshold_selection/north_case3.csv")

# Case 4: partially-observed GPD
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=1000, ncol=500)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-rbeta(length(data_all),1,2)
  data_above <- sample(data_all[data_all > 1], 279, replace=FALSE)
  data_below <- sample(data_all[data_all > cens_thr & data_all <= 1], 721, replace = FALSE)
  data <- c(data_below, data_above)
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


eqdthr_case4 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case4 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case4 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(data_matrix, "data/Case_4.csv")
write.csv(eqdthr_case4, "output/threshold_selection/eqd_case4.csv")
write.csv(wadsthr_case4, "output/threshold_selection/wads_case4.csv")
write.csv(norththr_case4, "output/threshold_selection/north_case4.csv")
