source("src/eqd_noboot.R")
library(evir)
library(Metrics)

#Case 0
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd_noboot(data, thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
}

eqdthr_case0 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case0_noboot.csv"
write.csv(eqdthr_case0, filename)


#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
    print(ii)
    set.seed(ii)
    dat1 <- runif(200, 0.5, 1.0)
    dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
    data <- c(dat1, dat2)
    thresh <- quantile(data, probs, names=F)
    
    #EQD method
    myres <- eqd_noboot(data, thresh)
    mythresh[ii] <- myres$thresh
    myquantile[ii] <- probs[thresh==mythresh[ii]]
    myscale[ii] <- myres$par[1]
    myshape[ii] <- myres$par[2]
    mylen[ii] <- myres$num_excess
}

eqdthr_case1 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
  
filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case1_noboot.csv"
write.csv(eqdthr_case1, filename)

#Case 2
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd_noboot(data, thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
}

eqdthr_case2 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case2_noboot.csv"
write.csv(eqdthr_case2, filename)

#Case 3
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(400, 0.5, 1.0)
  dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd_noboot(data, thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
}

eqdthr_case3 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case3_noboot.csv"
write.csv(eqdthr_case3, filename)

#Case 4
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-rbeta(length(data_all),1,2)
  data_above <- sample(data_all[data_all > 1], 279, replace=FALSE)
  data_below <- sample(data_all[data_all > cens_thr & data_all <= 1], 721, replace = FALSE)
  data <- c(data_below, data_above)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd_noboot(data, thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
}

eqdthr_case4 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case4_noboot.csv"
write.csv(eqdthr_case4, filename)


#Gaussian case
n <- 500
probs=seq(0.5, 0.95, by=0.05)
mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data <- rnorm(2000)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd_noboot(data, thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
}

eqdthr_gauss <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

filename <- "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_gauss_noboot.csv"
write.csv(eqdthr_gauss, filename)

library(Metrics)

rmse(eqdthr_case1$thr,1)
bias(eqdthr_case1$thr,1)
var(eqdthr_case1$thr)
hist(eqdthr_case1$thr)

rmse(eqdthr_case4$thr,1)
bias(eqdthr_case4$thr,1)
var(eqdthr_case4$thr)
hist(eqdthr_case4$thr)


#Nidd data

data("nidd.thresh")

#Candidate threshold grids
thresholds_list <- list(quantile(nidd.thresh, seq(0,0.93, by=0.01)),quantile(nidd.thresh, seq(0,0.9, by=0.01)), quantile(nidd.thresh, seq(0,0.8, by=0.01)), 
                        quantile(nidd.thresh, seq(0,0.8, by=0.2)), quantile(nidd.thresh, seq(0,0.9, by=0.3)), quantile(nidd.thresh, seq(0,0.75, by=0.25)),
                        quantile(nidd.thresh, c(0,0.1,0.4,0.7)))

EQD_thr <- numeric(7)
for(i in 1:7){
  thresholds <- thresholds_list[[i]]
  #EQD
  set.seed(11111) 
  EQD_thr[i] <- eqd_noboot(nidd.thresh, thresholds)$thresh
}

print(EQD_thr)

write.csv(EQD_thr, "output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_Nidd_noboot.csv")


#Reading in results
eqd_case0_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case0_noboot.csv", row.names = 1, header = T)
eqd_case1_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case1_noboot.csv", row.names = 1, header = T)
eqd_case2_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case2_noboot.csv", row.names = 1, header = T)
eqd_case3_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case3_noboot.csv", row.names = 1, header = T)
eqd_case4_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_case4_noboot.csv", row.names = 1, header = T)
eqd_gauss_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_gauss_noboot.csv", row.names = 1, header = T)
eqd_Nidd_noboot <- read.csv("output/threshold_selection/additional_sensitivity/no_bootstrapping/eqd_Nidd_noboot.csv", row.names = 1, header = T)


#Calculating RMSE, bias and variance
rmse(eqd_case0_noboot$thr, 1)
bias(eqd_case0_noboot$thr,1)
var(eqd_case0_noboot$thr)

rmse(eqd_case1_noboot$thr, 1)
bias(eqd_case1_noboot$thr,1)
var(eqd_case1_noboot$thr)

rmse(eqd_case2_noboot$thr, 1)
bias(eqd_case2_noboot$thr,1)
var(eqd_case2_noboot$thr)

rmse(eqd_case3_noboot$thr, 1)
bias(eqd_case3_noboot$thr,1)
var(eqd_case3_noboot$thr)

rmse(eqd_case4_noboot$thr, 1)
bias(eqd_case4_noboot$thr,1)
var(eqd_case4_noboot$thr)

#Gaussian
estimated_quantile <- function(df,p,n){
  lu <- df$len/n
  est_quants <- df$thr + (df$scale/df$shape)*((p/lu)^(-df$shape)-1)
  return(est_quants)
}

n <- 2000
p_seq <- c(1/n,1/(10*n), 1/(100*n))

eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_gauss_noboot, n=n)

true_quants <- qnorm(1-p_seq)

eqd_rmse_quants <- north_rmse_quants <- numeric(3)
for(i in 1:3){
  eqd_rmse_quants[i] <- rmse(eqd_est_quants[,i], true_quants[i])
}

RMSE_gauss <- data.frame(j=c(0,1,2), EQD=eqd_rmse_quants)

print(RMSE_gauss)
