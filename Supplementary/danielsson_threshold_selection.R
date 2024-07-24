library(tea)

source("src/helper_functions.R")

#Danielsson et al 2019 method
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
  return(list(u=u, k0=k_star))
}


#Case 1
length_data <- 1200  
dan19thr <- dan19shape <- dan19scale <- dan19len <- dan19quantile <- numeric(500)
dan01thr <- dan01shape <- dan01scale <- dan01len <- dan01quantile <- numeric(500)
probs=seq(0, 0.95, by=0.05)
k <-length_data*probs
data_matrix <- read.csv("data/Case_1.csv", row.names = 1, header=T) 
for(ii in 1:500){
  data <- data_matrix[,ii] 
  thresh <- quantile(data, probs, names=F)
  #Danielsson 2019
  dan19selection <- dan2019(data, k)
  dan19thr[ii] <- dan19selection$u
  dan19quantile[ii] <- 1 - (dan19selection$k0/length_data)
  fit.data <- data[data>dan19thr[ii]] - dan19thr[ii]
  optdan19 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan19shape[ii] <- optdan19$par[2]
  dan19scale[ii] <- optdan19$par[1]
  dan19len[ii] <- length(fit.data)
  
  #Danielsson 2001
  dan01selection <- danielsson(data)
  dan01thr[ii] <- dan01selection$threshold
  dan01quantile[ii] <- 1-(dan01selection$k0/length_data)
  fit.data <- data[data>dan01thr[ii]] - dan01thr[ii]
  optdan01 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan01shape[ii] <- optdan01$par[2]
  dan01scale[ii] <- optdan01$par[1]
  dan01len[ii] <- length(fit.data)
}

dan2019_case1 <- data.frame(thr=dan19thr,quantile=dan19quantile, scale=dan19scale, shape=dan19shape, len=dan19len)
write.csv(dan2019_case1, "output/threshold_selection/dan2019_case1.csv")

dan2001_case1 <- data.frame(thr=dan01thr,quantile=dan01quantile, scale=dan01scale, shape=dan01shape, len=dan01len)
write.csv(dan2001_case1, "output/threshold_selection/dan2019_case1.csv")

#Case 2
length_data <- 480 
dan19thr <- dan19shape <- dan19scale <- dan19len <- dan19quantile <- numeric(500)
dan01thr <- dan01shape <- dan01scale <- dan01len <- dan01quantile <- numeric(500)
probs=seq(0, 0.95, by=0.05)
k <-length_data*probs
data_matrix <- read.csv("data/Case_2.csv", row.names = 1, header=T) 
for(ii in 1:500){
  data <- data_matrix[,ii] 
  thresh <- quantile(data, probs, names=F)
  #Danielsson 2019
  dan19selection <- dan2019(data, k)
  dan19thr[ii] <- dan19selection$u
  dan19quantile[ii] <- 1 - (dan19selection$k0/length_data)
  fit.data <- data[data>dan19thr[ii]] - dan19thr[ii]
  optdan19 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan19shape[ii] <- optdan19$par[2]
  dan19scale[ii] <- optdan19$par[1]
  dan19len[ii] <- length(fit.data)
  
  #Danielsson 2001
  dan01selection <- danielsson(data)
  dan01thr[ii] <- dan01selection$threshold
  dan01quantile[ii] <- 1-(dan01selection$k0/length_data)
  fit.data <- data[data>dan01thr[ii]] - dan01thr[ii]
  optdan01 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan01shape[ii] <- optdan01$par[2]
  dan01scale[ii] <- optdan01$par[1]
  dan01len[ii] <- length(fit.data)
}

dan2019_case2 <- data.frame(thr=dan19thr,quantile=dan19quantile, scale=dan19scale, shape=dan19shape, len=dan19len)
write.csv(dan2019_case2, "output/threshold_selection/dan2019_case2.csv")

dan2001_case2 <- data.frame(thr=dan01thr,quantile=dan01quantile, scale=dan01scale, shape=dan01shape, len=dan01len)
write.csv(dan2001_case2, "output/threshold_selection/dan2019_case2.csv")


#Case 3
length_data <- 2400 
dan19thr <- dan19shape <- dan19scale <- dan19len <- dan19quantile <- numeric(500)
dan01thr <- dan01shape <- dan01scale <- dan01len <- dan01quantile <- numeric(500)
probs=seq(0, 0.95, by=0.05)
k <-length_data*probs
data_matrix <- read.csv("data/Case_3.csv", row.names = 1, header=T) 
for(ii in 1:500){
  data <- data_matrix[,ii] 
  thresh <- quantile(data, probs, names=F)
  #Danielsson 2019
  dan19selection <- dan2019(data, k)
  dan19thr[ii] <- dan19selection$u
  dan19quantile[ii] <- 1 - (dan19selection$k0/length_data)
  fit.data <- data[data>dan19thr[ii]] - dan19thr[ii]
  optdan19 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan19shape[ii] <- optdan19$par[2]
  dan19scale[ii] <- optdan19$par[1]
  dan19len[ii] <- length(fit.data)
  
  #Danielsson 2001
  dan01selection <- danielsson(data)
  dan01thr[ii] <- dan01selection$threshold
  dan01quantile[ii] <- 1-(dan01selection$k0/length_data)
  fit.data <- data[data>dan01thr[ii]] - dan01thr[ii]
  optdan01 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan01shape[ii] <- optdan01$par[2]
  dan01scale[ii] <- optdan01$par[1]
  dan01len[ii] <- length(fit.data)
}

dan2019_case3 <- data.frame(thr=dan19thr,quantile=dan19quantile, scale=dan19scale, shape=dan19shape, len=dan19len)
write.csv(dan2019_case3, "output/threshold_selection/dan2019_case3.csv")

dan2001_case3 <- data.frame(thr=dan01thr,quantile=dan01quantile, scale=dan01scale, shape=dan01shape, len=dan01len)
write.csv(dan2001_case3, "output/threshold_selection/dan2019_case3.csv")

#Case 4
length_data <- 1000 
dan19thr <- dan19shape <- dan19scale <- dan19len <- dan19quantile <- numeric(500)
dan01thr <- dan01shape <- dan01scale <- dan01len <- dan01quantile <- numeric(500)
probs=seq(0, 0.95, by=0.05)
k <-length_data*probs
data_matrix <- read.csv("data/Case_4.csv", row.names = 1, header=T) 
for(ii in 1:500){
  data <- data_matrix[,ii] 
  thresh <- quantile(data, probs, names=F)
  #Danielsson 2019
  dan19selection <- dan2019(data, k)
  dan19thr[ii] <- dan19selection$u
  dan19quantile[ii] <- 1 - (dan19selection$k0/length_data)
  fit.data <- data[data>dan19thr[ii]] - dan19thr[ii]
  optdan19 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan19shape[ii] <- optdan19$par[2]
  dan19scale[ii] <- optdan19$par[1]
  dan19len[ii] <- length(fit.data)
  
  #Danielsson 2001
  dan01selection <- danielsson(data)
  dan01thr[ii] <- dan01selection$threshold
  dan01quantile[ii] <- 1-(dan01selection$k0/length_data)
  fit.data <- data[data>dan01thr[ii]] - dan01thr[ii]
  optdan01 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan01shape[ii] <- optdan01$par[2]
  dan01scale[ii] <- optdan01$par[1]
  dan01len[ii] <- length(fit.data)
}

dan2019_case4 <- data.frame(thr=dan19thr,quantile=dan19quantile, scale=dan19scale, shape=dan19shape, len=dan19len)
write.csv(dan2019_case4, "output/threshold_selection/dan2019_case4.csv")

dan2001_case4 <- data.frame(thr=dan01thr,quantile=dan01quantile, scale=dan01scale, shape=dan01shape, len=dan01len)
write.csv(dan2001_case4, "output/threshold_selection/dan2019_case4.csv")

#Gaussian Case
length_data <- 2000 
dan19thr <- dan19shape <- dan19scale <- dan19len <- dan19quantile <- numeric(500)
dan01thr <- dan01shape <- dan01scale <- dan01len <- dan01quantile <- numeric(500)
probs=seq(0.05, 0.5, by=0.05)
k <-length_data*probs
data_matrix <- read.csv("data/Gaussian_case.csv", row.names = 1, header=T) 
for(ii in 1:500){
  data <- data_matrix[,ii]
  thresh <- quantile(data, probs, names=F)
  #Danielsson 2019
  dan19selection <- dan2019(data, k)
  dan19thr[ii] <- dan19selection$u
  dan19quantile[ii] <- 1 - (dan19selection$k0/length_data)
  fit.data <- data[data>dan19thr[ii]] - dan19thr[ii]
  optdan19 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan19shape[ii] <- optdan19$par[2]
  dan19scale[ii] <- optdan19$par[1]
  dan19len[ii] <- length(fit.data)
  
  #Danielsson 2001
  dan01selection <- danielsson(data)
  dan01thr[ii] <- dan01selection$threshold
  dan01quantile[ii] <- 1-(dan01selection$k0/length_data)
  fit.data <- data[data>dan01thr[ii]] - dan01thr[ii]
  optdan01 <- optim(GPD_LL, z=fit.data, par=c(0.1,0.1), control = list(fnscale=-1))
  dan01shape[ii] <- optdan01$par[2]
  dan01scale[ii] <- optdan01$par[1]
  dan01len[ii] <- length(fit.data)
}

dan2019_gauss <- data.frame(thr=dan19thr,quantile=dan19quantile, scale=dan19scale, shape=dan19shape, len=dan19len)
write.csv(dan2019_gauss, "output/threshold_selection/dan2019_gauss.csv")

dan2001_gauss <- data.frame(thr=dan01thr,quantile=dan01quantile, scale=dan01scale, shape=dan01shape, len=dan01len)
write.csv(dan2001_gauss, "output/threshold_selection/dan2019_gauss.csv")
