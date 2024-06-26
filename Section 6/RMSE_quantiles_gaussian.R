library(Metrics)

#Threshold choices
eqd_gauss <- read.csv("output/threshold_selection/eqd_gauss.csv", row.names = 1, header = T)
wads_gauss <- read.csv("output/threshold_selection/wads_gauss.csv", row.names = 1, header = T)
north_gauss <- read.csv("output/threshold_selection/north_gauss.csv", row.names = 1, header = T)

eqd_gauss_large <- read.csv("output/threshold_selection/eqd_gauss_large.csv", row.names = 1, header = T)
wads_gauss_large <- read.csv("output/threshold_selection/wads_gauss_large.csv", row.names = 1, header = T)
north_gauss_large <- read.csv("output/threshold_selection/north_gauss_large.csv", row.names = 1, header = T)

eqd_list <- list(eqd_gauss, eqd_gauss_large)
wads_list <- list(wads_gauss, wads_gauss_large)
north_list <- list(north_gauss, north_gauss_large)

#------------------Estimated quantiles------------------
estimated_quantile <- function(df,p,n){
  lu <- df$len/n
  est_quants <- df$thr + (df$scale/df$shape)*((p/lu)^(-df$shape)-1)
  return(est_quants)
}

#---------RMSE of estimated quantiles for Gaussian case--------------------
RMSE_gauss <- vector('list',2)
sample_sizes <- c(2000,20000)
for(case in 1:2){
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  eqd_df <- eqd_list[[case]]
  wads_df <- wads_list[[case]]
  north_df <- north_list[[case]]
  
  eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_df, n=n)
  wads_est_quants <- sapply(p_seq, estimated_quantile, df=wads_df[!is.na(wads_df$thr),], n=n)
  north_est_quants <- sapply(p_seq, estimated_quantile, df=north_df, n=n)
  
  true_quants <- qnorm(1-p_seq)
  
  eqd_rmse_quants <- wads_rmse_quants <- north_rmse_quants <- numeric(3)
  for(i in 1:3){
    eqd_rmse_quants[i] <- rmse(eqd_est_quants[,i], true_quants[i])
    wads_rmse_quants[i] <- rmse(wads_est_quants[,i], true_quants[i])
    north_rmse_quants[i] <- rmse(north_est_quants[,i], true_quants[i])
  }
  
  RMSE_gauss[[case]] <- data.frame(j=c(0,1,2), EQD=eqd_rmse_quants, Wadsworth=wads_rmse_quants, Northrop=north_rmse_quants)
}

print(RMSE_gauss)

saveRDS(RMSE_gauss, "output/tables/Table_5_rmse_quantile_estimates_gaussian.rds")
