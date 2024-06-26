
source("src/eqd.R")
source("src/JointMLEFunctions.R")
library(threshr)
library(Metrics)

#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)
data_matrix <- matrix(NA, nrow=1200, ncol=500)
mode <- 1
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data[data>mode], probs, names=F)
  
  #EQD method
  myres <- eqd(data, thresh=thresh)
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


eqdthr_case1 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case1 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case1 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(eqdthr_case1, "output/threshold_selection/additional_sensitivity/mode/eqd_case1_mode.csv")
write.csv(wadsthr_case1, "output/threshold_selection/additional_sensitivity/mode/wads_case1_mode.csv")
write.csv(norththr_case1, "output/threshold_selection/additional_sensitivity/mode/north_case1_mode.csv")


# Case 4: partially-observed GPD
mythresh <- wadsthresh <- norththresh <- numeric(n)
myquantile <- wadsquantile <- northquantile <- numeric(n)
myscale <- wadsscale <- northscale <- numeric(n)
myshape <- wadsshape <- northshape <- numeric(n)
mylen <- wadslen <- northlen <- numeric(n)

case_4_f_prime <- function(x, par=c(0.5,0.1), beta_par=c(1,2)){
  neg_h_prime_over_h <- ((1+par[2])/par[1])*(1 + par[2]*x/par[1])^(-1)
  f_prime <- (neg_h_prime_over_h*pbeta(x, shape1=beta_par[1], shape2=beta_par[2]) - dbeta(x, shape1=beta_par[1], shape2=beta_par[2]))^2
  return(f_prime)
}

mode_solver <- optimize(case_4_f_prime, c(0,1))
mode_case4 <- mode_solver$minimum
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-rbeta(length(data_all),1,2)
  data_above <- sample(data_all[data_all > 1], 279, replace=FALSE)
  data_below <- sample(data_all[data_all > cens_thr & data_all <= 1], 721, replace = FALSE)
  data <- c(data_below, data_above)

  thresh <- quantile(data[data>mode_case4], probs, names=F)
  
  #EQD method
  myres <- eqd(data, thresh=thresh)
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


eqdthr_case4 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
wadsthr_case4 <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)
norththr_case4 <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(eqdthr_case4, "output/threshold_selection/additional_sensitivity/mode/eqd_case4_mode.csv")
write.csv(wadsthr_case4, "output/threshold_selection/additional_sensitivity/mode/wads_case4_mode.csv")
write.csv(norththr_case4, "output/threshold_selection/additional_sensitivity/mode/north_case4_mode.csv")


# Gaussian case using thresholds which lie below the mode

#n=2000
n <- 500
probs=seq(0, 0.95, by=0.05)
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
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  myres <- eqd(data, thresh=thresh)
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num_excess
  
  # #Wadsworth 2016 method
  # wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  # if(!is.na(wadsthresh[ii])){
  #   fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
  #   optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
  #   wadsscale[ii] <- optwads$par[1]
  #   wadsshape[ii] <- optwads$par[2]
  #   wadslen[ii] <- length(fit.data_w)
  #   wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  # }
  # else{
  #   wadsscale[ii] <- NA
  #   wadsshape[ii] <- NA
  #   wadslen[ii] <- NA
  #   wadsquantile[ii] <- NA
  # }
  
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
norththr_gauss <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)

write.csv(eqdthr_gauss, "output/threshold_selection/additional_sensitivity/mode/eqd_gauss_mode.csv")
write.csv(norththr_gauss, "output/threshold_selection/additional_sensitivity/mode/north_gauss_mode.csv")

#Reading in results
eqdthr_case1_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/eqd_case1_mode.csv", row.names = 1, header=T)
wadsthr_case1_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/wads_case1_mode.csv", row.names = 1, header=T)
norththr_case1_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/north_case1_mode.csv", row.names = 1, header=T)

eqdthr_case4_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/eqd_case4_mode.csv", row.names = 1, header=T)
norththr_case4_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/north_case4_mode.csv", row.names = 1, header=T)

eqdthr_gauss_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/eqd_gauss_mode.csv", row.names = 1, header=T)
norththr_gauss_mode <- read.csv("output/threshold_selection/additional_sensitivity/mode/north_gauss_mode.csv", row.names = 1, header=T)

#RMSE, bias and var calculations
#Case 1
rmse(eqdthr_case1_mode$thr,1)
bias(eqdthr_case1_mode$thr,1)
var(eqdthr_case1_mode$thr)

err <- is.na(wadsthr_case1_mode$thr)
rmse(wadsthr_case1_mode$thr[!err],1)
bias(wadsthr_case1_mode$thr[!err],1)
var(wadsthr_case1_mode$thr[!err])

rmse(norththr_case1_mode$thr,1)
bias(norththr_case1_mode$thr,1)
var(norththr_case1_mode$thr)

#Case 4
rmse(eqdthr_case4_mode$thr,1)
bias(eqdthr_case4_mode$thr,1)
var(eqdthr_case4_mode$thr)

rmse(norththr_case4_mode$thr,1)
bias(norththr_case4_mode$thr,1)
var(norththr_case4_mode$thr)

#Gaussian
estimated_quantile <- function(df,p,n){
  lu <- df$len/n
  est_quants <- df$thr + (df$scale/df$shape)*((p/lu)^(-df$shape)-1)
  return(est_quants)
}
#---------RMSE of estimated quantiles for Gaussian case--------------------

n <- 2000
p_seq <- c(1/n,1/(10*n), 1/(100*n))

eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqdthr_gauss_mode, n=n)
north_est_quants <- sapply(p_seq, estimated_quantile, df=norththr_gauss_mode, n=n)
  
true_quants <- qnorm(1-p_seq)
  
eqd_rmse_quants <- north_rmse_quants <- numeric(3)
for(i in 1:3){
  eqd_rmse_quants[i] <- rmse(eqd_est_quants[,i], true_quants[i])
  north_rmse_quants[i] <- rmse(north_est_quants[,i], true_quants[i])
}
  
RMSE_gauss <- data.frame(j=c(0,1,2), EQD=eqd_rmse_quants, Northrop=north_rmse_quants)

print(RMSE_gauss)


