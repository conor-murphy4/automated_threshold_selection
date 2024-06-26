source('src/eqd_Varty.R')

#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
varty_thresh <- varty_quantile <- varty_scale <- varty_shape <- varty_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #varty method
  varty_res <- eqd_varty(data, thresh=thresh)
  varty_thresh[ii] <- varty_res$thresh
  varty_quantile[ii] <- probs[thresh==varty_thresh[ii]]
  varty_scale[ii] <- varty_res$par[1]
  varty_shape[ii] <- varty_res$par[2]
  varty_len[ii] <- varty_res$num_excess
}

varty_thr_case1 <- data.frame(thr=varty_thresh,quantile=varty_quantile, scale=varty_scale, shape=varty_shape, len=varty_len)

write.csv(varty_thr_case1, "output/threshold_selection/varty_case1.csv")

#Case 2: small sample size
varty_thresh <- varty_quantile <- varty_scale <- varty_shape <- varty_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(80, 0.5, 1.0)
  dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  varty_res <- eqd_varty(data, thresh=thresh)
  varty_thresh[ii] <- varty_res$thresh
  varty_quantile[ii] <- probs[thresh==varty_thresh[ii]]
  varty_scale[ii] <- varty_res$par[1]
  varty_shape[ii] <- varty_res$par[2]
  varty_len[ii] <- varty_res$num_excess
}


varty_thr_case2 <- data.frame(thr=varty_thresh,quantile=varty_quantile, scale=varty_scale, shape=varty_shape, len=varty_len)

write.csv(varty_thr_case2, "output/threshold_selection/varty_case2.csv")


#Case 3: negative shape parameter
varty_quantile <- varty_thresh <- varty_scale <- varty_shape <- varty_len <- numeric(n)
for(ii in 1:n){
  print(ii)
  set.seed(ii)
  dat1 <- runif(400, 0.5, 1.0)
  dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, probs, names=F)
  
  #EQD method
  varty_res <- eqd_varty(data, thresh=thresh)
  varty_thresh[ii] <- varty_res$thresh
  varty_quantile[ii] <- probs[thresh==varty_thresh[ii]]
  varty_scale[ii] <- varty_res$par[1]
  varty_shape[ii] <- varty_res$par[2]
  varty_len[ii] <- varty_res$num_excess
}


varty_thr_case3 <- data.frame(thr=varty_thresh,quantile=varty_quantile, scale=varty_scale, shape=varty_shape, len=varty_len)

write.csv(varty_thr_case3, "output/threshold_selection/varty_case3.csv")


# Case 4: partially-observed GPD
u <- 1.0
varty_thresh <- varty_quantile <- varty_scale <- varty_shape <- varty_len <- numeric(n)
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
  varty_res <- eqd_varty(data, thresh=thresh)
  varty_thresh[ii] <- varty_res$thresh
  varty_quantile[ii] <- probs[thresh==varty_thresh[ii]]
  varty_scale[ii] <- varty_res$par[1]
  varty_shape[ii] <- varty_res$par[2]
  varty_len[ii] <- varty_res$num_excess
  
}

varty_thr_case4 <- data.frame(thr=varty_thresh,quantile=varty_quantile, scale=varty_scale, shape=varty_shape, len=varty_len)

write.csv(varty_thr_case4, "output/threshold_selection/varty_case4.csv")

#Reading in results
eqd_case1 <- read.csv("output/threshold_selection/eqd_case1.csv", row.names = 1, header = T)
eqd_case2 <- read.csv("output/threshold_selection/eqd_case2.csv", row.names = 1, header = T)
eqd_case3 <- read.csv("output/threshold_selection/eqd_case3.csv", row.names = 1, header = T)
eqd_case4 <- read.csv("output/threshold_selection/eqd_case4.csv", row.names = 1, header = T)

varty_case1 <- read.csv("output/threshold_selection/varty_case1.csv", row.names = 1, header = T)
varty_case2 <- read.csv("output/threshold_selection/varty_case2.csv", row.names = 1, header = T)
varty_case3 <- read.csv("output/threshold_selection/varty_case3.csv", row.names = 1, header = T)
varty_case4 <- read.csv("output/threshold_selection/varty_case4.csv", row.names = 1, header = T)

# RMSE, bias, var results
library(Metrics)

rmse(varty_case4$thr, 1)
bias(varty_case4$thr, 1)
var(varty_case4$thr)
#------------------True quantiles---------------------

#Cases 1-3 , p exceedance prob, par (sig, xi), u true
case123_true_quant <- function(p,par,u){
  x_p <- (par[1]/par[2])*((6*p/5)^(-par[2])-1) + u
  return(x_p)
}
#Case 4, p exceedance prob, p_1 P(<1), par (sig, xi), u true
case4_true_quant <- function(p, par, p_1, u){
  x_p <- ((par[1]+par[2]*(u))/par[2])*((p/(1-p_1))^(-par[2])-1) + u
  return(x_p)  
}

#------------------Estimated quantiles------------------
estimated_quantile <- function(df,p,n){
  lu <- df$len/n
  est_quants <- df$thr + (df$scale/df$shape)*((p/lu)^(-df$shape)-1)
  return(est_quants)
}

#---------RMSE of estimated quantiles for Case 1-4--------------------

eqd_list <- list(eqd_case1, eqd_case2, eqd_case3, eqd_case4)
varty_list <- list(varty_case1, varty_case2, varty_case3, varty_case4)


RMSE_case1234 <- vector('list',4)
sample_sizes <- c(1200, 480, 2400, 1000)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
for(case in 1:4){
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  par <- true_pars[case,]
  eqd_df <- eqd_list[[case]]
  varty_df <- varty_list[[case]]
  
  eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_df, n=n)
  varty_est_quants <- sapply(p_seq, estimated_quantile, df=varty_df, n=n)
  
  if(case != 4){
    true_quants <- sapply(p_seq, case123_true_quant, par=par, u=1)
  }
  else{
    true_quants <- sapply(p_seq, case4_true_quant, par=par, p_1= 0.7206618, u=1)
  }
  
  eqd_rmse_quants <- varty_rmse_quants <- numeric(3)
  for(i in 1:3){
    eqd_rmse_quants[i] <- bias(eqd_est_quants[,i], true_quants[i])
    varty_rmse_quants[i] <- bias(varty_est_quants[,i], true_quants[i])
  }
  
  RMSE_case1234[[case]] <- data.frame(j=c(0,1,2), EQD=eqd_rmse_quants, Varty = varty_rmse_quants)
}

print(RMSE_case1234)

apply(eqd_est_quants, 2, var)
apply(varty_est_quants, 2, var)

