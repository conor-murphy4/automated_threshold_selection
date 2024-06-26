source("src/eqd_Qpx.R")

#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
Qpx_thresh <- Qpx_quantile <- Qpx_scale <- Qpx_shape <- Qpx_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #Qpx method
  Qpx_res <- eqd_Qpx(data, thresh=thresh)
  Qpx_thresh[ii] <- Qpx_res$thresh
  Qpx_quantile[ii] <- probs[thresh==Qpx_thresh[ii]]
  Qpx_scale[ii] <- Qpx_res$par[1]
  Qpx_shape[ii] <- Qpx_res$par[2]
  Qpx_len[ii] <- Qpx_res$num_excess
}

Qpx_thr_case1 <- data.frame(thr=Qpx_thresh,quantile=Qpx_quantile, scale=Qpx_scale, shape=Qpx_shape, len=Qpx_len)

write.csv(Qpx_thr_case1, "output/threshold_selection/eqd_Qpx_case1.csv")

#Case 2: small sample size
Qpx_thresh <- Qpx_quantile <- Qpx_scale <- Qpx_shape <- Qpx_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  Qpx_res <- eqd_Qpx(data, thresh=thresh)
  Qpx_thresh[ii] <- Qpx_res$thresh
  Qpx_quantile[ii] <- probs[thresh==Qpx_thresh[ii]]
  Qpx_scale[ii] <- Qpx_res$par[1]
  Qpx_shape[ii] <- Qpx_res$par[2]
  Qpx_len[ii] <- Qpx_res$num_excess
}


Qpx_thr_case2 <- data.frame(thr=Qpx_thresh,quantile=Qpx_quantile, scale=Qpx_scale, shape=Qpx_shape, len=Qpx_len)

write.csv(Qpx_thr_case2, "output/threshold_selection/eqd_Qpx_case2.csv")

#Case 3: negative shape parameter
Qpx_quantile <- Qpx_thresh <- Qpx_scale <- Qpx_shape <- Qpx_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(400, 0.5, 1.0)
  dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  Qpx_res <- eqd_Qpx(data, thresh=thresh)
  Qpx_thresh[ii] <- Qpx_res$thresh
  Qpx_quantile[ii] <- probs[thresh==Qpx_thresh[ii]]
  Qpx_scale[ii] <- Qpx_res$par[1]
  Qpx_shape[ii] <- Qpx_res$par[2]
  Qpx_len[ii] <- Qpx_res$num_excess
}


Qpx_thr_case3 <- data.frame(thr=Qpx_thresh,quantile=Qpx_quantile, scale=Qpx_scale, shape=Qpx_shape, len=Qpx_len)

write.csv(Qpx_thr_case3, "output/threshold_selection/eqd_Qpx_case3.csv")

# Case 4: partially-observed GPD
u <- 1.0
Qpx_thresh <- Qpx_quantile <- Qpx_scale <- Qpx_shape <- Qpx_len <- numeric(n)
set.seed(12345)
for(ii in 1:n){
  print(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-rbeta(length(data_all),1,2)
  data_above <- sample(data_all[data_all > 1], 279, replace=FALSE)
  data_below <- sample(data_all[data_all > cens_thr & data_all <= 1], 721, replace = FALSE)
  data <- c(data_below, data_above)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  Qpx_res <- eqd_Qpx(data, thresh=thresh)
  Qpx_thresh[ii] <- Qpx_res$thresh
  Qpx_quantile[ii] <- probs[thresh==Qpx_thresh[ii]]
  Qpx_scale[ii] <- Qpx_res$par[1]
  Qpx_shape[ii] <- Qpx_res$par[2]
  Qpx_len[ii] <- Qpx_res$num_excess
  
}

Qpx_thr_case4 <- data.frame(thr=Qpx_thresh,quantile=Qpx_quantile, scale=Qpx_scale, shape=Qpx_shape, len=Qpx_len)

write.csv(Qpx_thr_case4, "output/threshold_selection/eqd_Qpx_case4.csv")


#Reading in results
eqd_Qpx_case1 <- read.csv("output/threshold_selection/eqd_Qpx_case1.csv", row.names = 1, header = T)
eqd_Qpx_case2 <- read.csv("output/threshold_selection/eqd_Qpx_case2.csv", row.names = 1, header = T)
eqd_Qpx_case3 <- read.csv("output/threshold_selection/eqd_Qpx_case3.csv", row.names = 1, header = T)
eqd_Qpx_case4 <- read.csv("output/threshold_selection/eqd_Qpx_case4.csv", row.names = 1, header = T)

# RMSE, bias, var results
library(Metrics)

rmse(eqd_Qpx_case4$thr, 1)
bias(eqd_Qpx_case4$thr, 1)
var(eqd_Qpx_case4$thr)


