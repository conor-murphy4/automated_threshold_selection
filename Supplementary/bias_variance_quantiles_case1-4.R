library(Metrics)

#Threshold choices
eqd_case1 <- read.csv("output/threshold_selection/eqd_case1.csv", row.names = 1, header = T)
eqd_case2 <- read.csv("output/threshold_selection/eqd_case2.csv", row.names = 1, header = T)
eqd_case3 <- read.csv("output/threshold_selection/eqd_case3.csv", row.names = 1, header = T)
eqd_case4 <- read.csv("output/threshold_selection/eqd_case4.csv", row.names = 1, header = T)

wads_case1 <- read.csv("output/threshold_selection/wads_case1.csv", row.names = 1, header = T)
wads_case2 <- read.csv("output/threshold_selection/wads_case2.csv", row.names = 1, header = T)
wads_case3 <- read.csv("output/threshold_selection/wads_case3.csv", row.names = 1, header = T)
wads_case4 <- read.csv("output/threshold_selection/wads_case4.csv", row.names = 1, header = T)

north_case1 <- read.csv("output/threshold_selection/north_case1.csv", row.names = 1, header = T)
north_case2 <- read.csv("output/threshold_selection/north_case2.csv", row.names = 1, header = T)
north_case3 <- read.csv("output/threshold_selection/north_case3.csv", row.names = 1, header = T)
north_case4 <- read.csv("output/threshold_selection/north_case4.csv", row.names = 1, header = T)

eqd_list <- list(eqd_case1, eqd_case2, eqd_case3, eqd_case4)
wads_list <- list(wads_case1, wads_case2, wads_case3, wads_case4)
north_list <- list(north_case1, north_case2, north_case3, north_case4)

#------------------True quantiles---------------------

#Ideal case , p exceedance prob, par (sig, xi), u true
case123_true_quant <- function(p,par,u){
  x_p <- (par[1]/par[2])*((6*p/5)^(-par[2])-1) + u
  return(x_p)
}

#Censored case, p exceedance prob, p_1 P(<1), par (sig, xi), u true
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

#---------bias of estimated quantiles for Case 1-4--------------------
bias_case1234 <- list(4)
sample_sizes <- c(1200, 480, 2400, 1000)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
for(case in 1:4){
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  par <- true_pars[case,]
  eqd_df <- eqd_list[[case]]
  wads_df <- wads_list[[case]]
  north_df <- north_list[[case]]
  
  eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_df, n=n)
  wads_est_quants <- sapply(p_seq, estimated_quantile, df=wads_df[!is.na(wads_df$thr),], n=n)
  north_est_quants <- sapply(p_seq, estimated_quantile, df=north_df, n=n)
  
  if(case != 4){
    true_quants <- sapply(p_seq, case123_true_quant, par=par, u=1)
  }
  else{
    true_quants <- sapply(p_seq, case4_true_quant, par=par, p_1= 0.7206618, u=1)
  }
  
  eqd_bias_quants <- wads_bias_quants <- north_bias_quants <- numeric(3)
  for(i in 1:3){
    eqd_bias_quants[i] <- bias(eqd_est_quants[,i], true_quants[i])
    wads_bias_quants[i] <- bias(wads_est_quants[,i], true_quants[i])
    north_bias_quants[i] <- bias(north_est_quants[,i], true_quants[i])
  }
  
  bias_case1234[[case]] <- data.frame(j=c(0,1,2), EQD=eqd_bias_quants, Wadsworth=wads_bias_quants, Northrop=north_bias_quants)
}

print(bias_case1234)

saveRDS(bias_case1234, "output/tables/Table_S.9_bias_quantile_estimates_case1-4.rds")


#---------var of estimated quantiles for Case 1-4--------------------
var_case1234 <- list(4)
sample_sizes <- c(1200, 480, 2400, 1000)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
for(case in 1:4){
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  par <- true_pars[case,]
  eqd_df <- eqd_list[[case]]
  wads_df <- wads_list[[case]]
  north_df <- north_list[[case]]
  
  eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_df, n=n)
  wads_est_quants <- sapply(p_seq, estimated_quantile, df=wads_df[!is.na(wads_df$thr),], n=n)
  north_est_quants <- sapply(p_seq, estimated_quantile, df=north_df, n=n)
  
  
    eqd_var_quants <- apply(eqd_est_quants, 2, var)
    wads_var_quants <- apply(wads_est_quants, 2, var)
    north_var_quants <- apply(north_est_quants, 2, var)
  
  var_case1234[[case]] <- data.frame(j=c(0,1,2), EQD=eqd_var_quants, Wadsworth=wads_var_quants, Northrop=north_var_quants)
}

print(var_case1234)

saveRDS(var_case1234, "output/tables/Table_S.9_variance_quantile_estimates_case1-4.rds")
