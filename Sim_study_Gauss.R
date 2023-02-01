setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
source("STORM_input/thr_selection_final.R")
source("STORM_input/Wadsworth2016.R")
library(threshr)

#Gaussian data
n4 <- 500
Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n4)
mythresh <- wadsthresh <- norththresh <- numeric(n4)
myscale <- wadsscale <- northscale <- numeric(n4)
myshape <- wadsshape <- northshape <- numeric(n4)
mylen <- wadslen <- northlen <- numeric(n4)
for(ii in 1:n4){
  set.seed(ii)
  data <- rnorm(1000)
  thresh <- quantile(data, probs=seq(0.5, 0.95, by=0.05), names=FALSE)
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

mythrG <- data.frame(thr=mythresh, scale=myscale, shape=myshape, len=mylen, time=Time_mythr)
wadsthrG <- data.frame(thr=wadsthresh, scale=wadsscale, shape=wadsshape, len=wadslen, time=Time_wads)
norththrG <- data.frame(thr=norththresh, scale=northscale, shape=northshape, len=northlen, time=Time_north, timeX = Time_north_x)

saveRDS(mythrG, "STORM_output/mythrG.rds")
saveRDS(wadsthrG, "STORM_output/wadsthrG.rds" )
saveRDS(norththrG, "STORM_output/norththrG.rds")

dev.new(width=9.17, height=4,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
hist(mythrG$thr, breaks=20, main="Our method", xlab="Threshold")
hist(wadsthrG$thr[!is.na(wadsthrG$thr)], main="Wadsworth",xlab="Threshold")
hist(norththrG$thr, breaks=20, main="Northrop", xlab="Threshold")
