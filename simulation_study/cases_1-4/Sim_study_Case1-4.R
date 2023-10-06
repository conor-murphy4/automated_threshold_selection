#setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
source("thr_selection_final.R")
source("JointMLEFunctions.r")
library(threshr)
library(tea)

#Idealistic
n <- 500
probs=seq(0, 0.95, by=0.05)
Time_mythr <- Time_wads <- Time_wads_x <- Time_north <- Time_north_x <- numeric(n)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=1200, ncol=500)
#set.seed(12345)
for(ii in 1:n){
  #saveRDS(ii, "STORM_output/iter_2.rds")
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  #My method
  myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  # start_wads <- Sys.time()
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  #end_wads <- Sys.time()
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
  # Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  # start_north <- Sys.time()
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}


mythrI <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)#, time=Time_mythr)
wadsthrI <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)#, time=Time_wads)
norththrI <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)#, time=Time_north, timeX = Time_north_x)

saveRDS(data_matrix, "Rerun with quantiles/data_sim_study_Case_1.rds")
saveRDS(mythrI, "Rerun with quantiles/mythrI.rds")
saveRDS(wadsthrI, "Rerun with quantiles/wadsthrI.rds")
saveRDS(norththrI, "Rerun with quantiles/norththrI.rds")

#Idealistic with small sample size
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
  #My method
  myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  # start_wads <- Sys.time()
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  #end_wads <- Sys.time()
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
  # Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  # start_north <- Sys.time()
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}


mythrI1 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)#, time=Time_mythr)
wadsthrI1 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)#, time=Time_wads)
norththrI1 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)#, time=Time_north, timeX = Time_north_x)

saveRDS(data_matrix, "Rerun with quantiles/data_sim_study_Case_2.rds")
saveRDS(mythrI1, "Rerun with quantiles/mythrI1.rds")
saveRDS(wadsthrI1, "Rerun with quantiles/wadsthrI1.rds")
saveRDS(norththrI1, "Rerun with quantiles/norththrI1.rds")

#Idealistic with negative shape
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
  #My method
  myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  # start_wads <- Sys.time()
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  #end_wads <- Sys.time()
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
  # Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  # start_north <- Sys.time()
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}


mythrI2 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)#, time=Time_mythr)
wadsthrI2 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)#, time=Time_wads)
norththrI2 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)#, time=Time_north, timeX = Time_north_x)

saveRDS(data_matrix, "Rerun with quantiles/data_sim_study_Case_3.rds")
saveRDS(mythrI2, "Rerun with quantiles/mythrI2.rds")
saveRDS(wadsthrI2, "Rerun with quantiles/wadsthrI2.rds")
saveRDS(norththrI2, "Rerun with quantiles/norththrI2.rds")

# #GPD data with censoring

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
  #saveRDS(ii, "STORM_output/iter_2.rds")
  #set.seed(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-u*rbeta(length(data_all),1,0.5)
  keep <- data_all>cens_thr
  data <- sample(data_all[keep], 1000, replace=FALSE)
  #data_matrix[,ii] <- data
  thresh <- quantile(data, probs, names=F)
  #My method
  myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  # start_wads <- Sys.time()
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  #end_wads <- Sys.time()
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
  # Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  # start_north <- Sys.time()
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}


mythrC <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)#, time=Time_mythr)
wadsthrC <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)#, time=Time_wads)
norththrC <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)#, time=Time_north, timeX = Time_north_x)

saveRDS(data_matrix, "Rerun with quantiles/data_sim_study_Case_4.rds")
saveRDS(mythrC, "Rerun with quantiles/mythrC.rds")
saveRDS(wadsthrC, "Rerun with quantiles/wadsthrC.rds")
saveRDS(norththrC, "Rerun with quantiles/norththrC.rds")
