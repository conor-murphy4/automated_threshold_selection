setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
source("STORM_input/thr_selection_final.R")
source("STORM_input/Wadsworth2016.R")
library(threshr)

#Idealistic 
n <- 500
Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
for(ii in 1:n){
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, seq(0,0.95,by=0.05))
  #My method
  start_mythr <- Sys.time()
  myres <- thr_dist(data, thresh=thresh)
  end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  start_wads <- Sys.time()
  wadsres <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)
  end_wads <- Sys.time()
  wadsthresh[ii] <- wadsres$thresh
  wadsscale[ii] <- wadsres$mle.u[2]
  wadsshape[ii] <- wadsres$mle.u[3]
  wadslen[ii] <- length(data[data>wadsres$thresh])
  Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  start_north <- Sys.time()
  norththr <- summary(ithresh(data, u_vec = thresh))[3]
  end_north1 <- Sys.time()
  fit.data <- data[data>norththr]
  optnorth <- optim(GPD_LL, z=fit.data, par=c(mean(fit.data), 0.1), control=list(fnscale=-1))
  end_north <- Sys.time()
  norththresh[ii] <- norththr
  northscale[ii] <- optnorth$par[[1]]
  northshape[ii] <- optnorth$par[[2]]
  northlen[ii] <- length(fit.data)
  Time_north[ii] <- end_north-start_north
  Time_north_x[ii] <- end_north1 - start_north
}


mythrI <- data.frame(thr=mythresh, scale=myscale, shape=myshape, len=mylen, time=Time_mythr)
wadsthrI <- data.frame(thr=wadsthresh, scale=wadsscale, shape=wadsshape, len=wadslen, time=Time_wads)
norththrI <- data.frame(thr=norththresh, scale=northscale, shape=northshape, len=northlen, time=Time_north, timeX = Time_north_x)

saveRDS(mythrI, "STORM_output/mythrI.rds")
saveRDS(wadsthrI, "STORM_output/wadsthrI.rds")
saveRDS(norththrI, "STORM_output/norththrI.rds")

#Idealistic with small sample size
n1 <- 500
Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n1)
mythresh <- wadsthresh <- norththresh <- numeric(n1)
myscale <- wadsscale <- northscale <- numeric(n1)
myshape <- wadsshape <- northshape <- numeric(n1)
mylen <- wadslen <- northlen <- numeric(n1)
for(ii in 1:n1){
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, seq(0,0.95,by=0.05))
  #My method
  start_mythr <- Sys.time()
  myres <- thr_dist(data, thresh=thresh)
  end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  start_wads <- Sys.time()
  wadsres <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)
  end_wads <- Sys.time()
  wadsthresh[ii] <- wadsres$thresh
  wadsscale[ii] <- wadsres$mle.u[2]
  wadsshape[ii] <- wadsres$mle.u[3]
  wadslen[ii] <- length(data[data>wadsres$thresh])
  Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  start_north <- Sys.time()
  norththr <- summary(ithresh(data, u_vec = thresh))[3]
  end_north1 <- Sys.time()
  fit.data <- data[data>norththr]
  optnorth <- optim(GPD_LL, z=fit.data, par=c(mean(fit.data), 0.1), control=list(fnscale=-1))
  end_north <- Sys.time()
  norththresh[ii] <- norththr
  northscale[ii] <- optnorth$par[[1]]
  northshape[ii] <- optnorth$par[[2]]
  northlen[ii] <- length(fit.data)
  Time_north[ii] <- end_north-start_north
  Time_north_x[ii] <- end_north1 - start_north
}


mythrI1 <- data.frame(thr=mythresh, scale=myscale, shape=myshape, len=mylen, time=Time_mythr)
wadsthrI1 <- data.frame(thr=wadsthresh, scale=wadsscale, shape=wadsshape, len=wadslen, time=Time_wads)
norththrI1 <- data.frame(thr=norththresh, scale=northscale, shape=northshape, len=northlen, time=Time_north, timeX = Time_north_x)

saveRDS(mythrI1, "STORM_output/mythrI1.rds")
saveRDS(wadsthrI1, "STORM_output/wadsthrI1.rds")
saveRDS(norththrI1, "STORM_output/norththrI1.rds")

#Idealistic with negative shape
n2 <- 500
Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n2)
mythresh <- wadsthresh <- norththresh <- numeric(n2)
myscale <- wadsscale <- northscale <- numeric(n2)
myshape <- wadsshape <- northshape <- numeric(n2)
mylen <- wadslen <- northlen <- numeric(n2)
for(ii in 1:n2){
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, seq(0,0.95,by=0.05))
  #My method
  start_mythr <- Sys.time()
  myres <- thr_dist(data, thresh=thresh)
  end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  start_wads <- Sys.time()
  wadsres <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)
  end_wads <- Sys.time()
  wadsthresh[ii] <- wadsres$thresh
  wadsscale[ii] <- wadsres$mle.u[2]
  wadsshape[ii] <- wadsres$mle.u[3]
  wadslen[ii] <- length(data[data>wadsres$thresh])
  Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  start_north <- Sys.time()
  norththr <- summary(ithresh(data, u_vec = thresh))[3]
  end_north1 <- Sys.time()
  fit.data <- data[data>norththr]
  optnorth <- optim(GPD_LL, z=fit.data, par=c(mean(fit.data), 0.1), control=list(fnscale=-1))
  end_north <- Sys.time()
  norththresh[ii] <- norththr
  northscale[ii] <- optnorth$par[[1]]
  northshape[ii] <- optnorth$par[[2]]
  northlen[ii] <- length(fit.data)
  Time_north[ii] <- end_north-start_north
  Time_north_x[ii] <- end_north1 - start_north
}


mythrI2 <- data.frame(thr=mythresh, scale=myscale, shape=myshape, len=mylen, time=Time_mythr)
wadsthrI2 <- data.frame(thr=wadsthresh, scale=wadsscale, shape=wadsshape, len=wadslen, time=Time_wads)
norththrI2 <- data.frame(thr=norththresh, scale=northscale, shape=northshape, len=northlen, time=Time_north, timeX = Time_north_x)

saveRDS(mythrI2, "STORM_output/mythrI2.rds")
saveRDS(wadsthrI2, "STORM_output/wadsthrI2.rds")
saveRDS(norththrI2, "STORM_output/norththrI2.rds")


# #GPD data with censoring

n3 <- 500
u <- 1.0
# Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n3)
# mythresh <- wadsthresh <- norththresh <- numeric(n3)
# myscale <- wadsscale <- northscale <- numeric(n3)
# myshape <- wadsshape <- northshape <- numeric(n3)
# mylen <- wadslen <- northlen <- numeric(n3)
n_below_u <- n_above_u <- numeric(n3)
for(ii in 1:n3){
  set.seed(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-u*rbeta(length(data_all),1,0.5)
  keep <- data_all>cens_thr
  data <- sample(data_all[keep],1000, replace=FALSE)
  n_below_u[ii] <- length(data[data <= u])
  n_above_u[ii] <- length(data[data > u])
  # thresh <- quantile(data, seq(0,0.95,by=0.05))
  # #My method
  # start_mythr <- Sys.time()
  # myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  # mythresh[ii] <- myres$thresh
  # myscale[ii] <- myres$par[1]
  # myshape[ii] <- myres$par[2]
  # mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  # #Wadsworth 2016 method
  # start_wads <- Sys.time()
  # wadsres <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)
  # end_wads <- Sys.time()
  # wadsthresh[ii] <- wadsres$thresh
  # wadsscale[ii] <- wadsres$mle.u[2]
  # wadsshape[ii] <- wadsres$mle.u[3]
  # wadslen[ii] <- length(data[data>wadsres$thresh])
  # Time_wads[ii] <- end_wads-start_wads
  # #Northrop 2017 method
  # start_north <- Sys.time()
  # norththr <- summary(ithresh(data, u_vec = thresh))[3]
  # end_north1 <- Sys.time()
  # fit.data <- data[data>norththr]
  # optnorth <- optim(GPD_LL, z=fit.data, par=c(mean(fit.data), 0.1), control=list(fnscale=-1))
  # end_north <- Sys.time()
  # norththresh[ii] <- norththr
  # northscale[ii] <- optnorth$par[[1]]
  # northshape[ii] <- optnorth$par[[2]]
  # northlen[ii] <- length(fit.data)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}

mean(n_below_u)
mean(n_above_u)
hist(data)
mythrC <- data.frame(thr=mythresh, scale=myscale, shape=myshape, len=mylen, time=Time_mythr)
wadsthrC <- data.frame(thr=wadsthresh, scale=wadsscale, shape=wadsshape, len=wadslen, time=Time_wads)
norththrC <- data.frame(thr=norththresh, scale=northscale, shape=northshape, len=northlen, time=Time_north, timeX = Time_north_x)

saveRDS(mythrC, "STORM_output/mythrC.rds")
saveRDS(wadsthrC, "STORM_output/wadsthrC.rds")
saveRDS(norththrC, "STORM_output/norththrC.rds")

dev.new(width=9.17, height=6.00,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')

dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data, breaks=25, main="Case 1", xlab='')
abline(v=1.0, col="red", lwd=2)
dat1 <- runif(80, 0.5, 1.0)
dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data, breaks=25,main="Case 2", xlab='')
abline(v=1.0, col="red", lwd=2)
dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=-0.05, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data, breaks=25,main="Case 3", xlab='')
abline(v=1.0, col="red", lwd=2)
data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
cens_thr<-rbeta(length(data_all),1,0.5)
keep <- data_all>cens_thr
data <- sample(data_all[keep],1000, replace = FALSE)
hist(data, breaks=25,main="Case 4", xlab='')
abline(v=1.0, col="red", lwd=2)

