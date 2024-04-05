library(Metrics)

#Threshold choices
dan2001_case1 <- read.csv("output/threshold_selection/dan2001_case1.csv", row.names = 1, header = T)
dan2001_case2 <- read.csv("output/threshold_selection/dan2001_case2.csv", row.names = 1, header = T)
dan2001_case3 <- read.csv("output/threshold_selection/dan2001_case3.csv", row.names = 1, header = T)
dan2001_case4 <- read.csv("output/threshold_selection/dan2001_case4.csv", row.names = 1, header = T)

dan2019_case1 <- read.csv("output/threshold_selection/dan2019_case1.csv", row.names = 1, header = T)
dan2019_case2 <- read.csv("output/threshold_selection/dan2019_case2.csv", row.names = 1, header = T)
dan2019_case3 <- read.csv("output/threshold_selection/dan2019_case3.csv", row.names = 1, header = T)
dan2019_case4 <- read.csv("output/threshold_selection/dan2019_case4.csv", row.names = 1, header = T)

dan2001_gaussian <- read.csv("output/threshold_selection/dan2001_gauss.csv", row.names = 1, header = T)
dan2019_gaussian <- read.csv("output/threshold_selection/dan2019_gauss.csv", row.names = 1, header = T)


#Table S.8: RMSE, bias, variance of threshold choices for Cases 1-4
(Dan2001_case1234 <- data.frame(RMSE=c(rmse(dan2001_case1$thr,1), rmse(dan2001_case2$thr,1), rmse(dan2001_case3$thr,1), rmse(dan2001_case4$thr,1)),Bias=c(bias(dan2001_case1$thr,1), bias(dan2001_case2$thr,1), bias(dan2001_case3$thr,1), bias(dan2001_case4$thr,1)), Variance = c(var(dan2001_case1$thr), var(dan2001_case2$thr), var(dan2001_case3$thr), var(dan2001_case4$thr))))

write.csv(Dan2001_case1234, "output/tables/Table_S.8_dan2001_rmse_bias_variance_threshold_choice_case1-4.csv")

(Dan2019_case1234 <- data.frame(RMSE=c(rmse(dan2019_case1$thr,1), rmse(dan2019_case2$thr,1), rmse(dan2019_case3$thr,1), rmse(dan2019_case4$thr,1)),Bias=c(bias(dan2019_case1$thr,1), bias(dan2019_case2$thr,1), bias(dan2019_case3$thr,1), bias(dan2019_case4$thr,1)), Variance = c(var(dan2019_case1$thr), var(dan2019_case2$thr), var(dan2019_case3$thr), var(dan2019_case4$thr))))

write.csv(Dan2019_case1234, "output/tables/Table_S.8_dan2019_rmse_bias_variance_threshold_choice_case1-4.csv")



#Table S.11: RMSE of quantile estimates for Cases 1-4

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
dan2001_list <- list(dan2001_case1, dan2001_case2, dan2001_case3, dan2001_case4)
dan2019_list <- list(dan2019_case1, dan2019_case2, dan2019_case3, dan2019_case4)

RMSE_dan_case1234 <- list(4)
sample_sizes <- c(1200, 480, 2400, 1000)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
for(case in 1:4){
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  par <- true_pars[case,]
  dan2001_df <- dan2001_list[[case]]
  dan2019_df <- dan2019_list[[case]]
  
  dan2001_est_quants <- sapply(p_seq, estimated_quantile, df=dan2001_df, n=n)
  dan2019_est_quants <- sapply(p_seq, estimated_quantile, df=dan2019_df, n=n)
  
  if(case != 4){
    true_quants <- sapply(p_seq, case123_true_quant, par=par, u=1)
  }
  else{
    true_quants <- sapply(p_seq, case4_true_quant, par=par, p_1= 0.7206618, u=1)
  }
  
  dan2001_rmse_quants <- dan2019_rmse_quants <- numeric(3)
  for(i in 1:3){
    dan2001_rmse_quants[i] <- rmse(dan2001_est_quants[,i], true_quants[i])
    dan2019_rmse_quants[i] <- rmse(dan2019_est_quants[,i], true_quants[i])
  }
  
  RMSE_dan_case1234[[case]] <- data.frame(j=c(0,1,2), Dan2001=dan2001_rmse_quants, Dan2019=dan2019_rmse_quants)
}

print(RMSE_dan_case1234)

saveRDS(RMSE_dan_case1234, "output/tables/Table_S.11_danielsson_rmse_quantile_estimates_case1-4.rds")
