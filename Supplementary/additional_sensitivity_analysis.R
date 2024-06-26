source("src/eqd.R")
library(Metrics)

k_vals <- c(10, 20, 50, 80)
for (k in k_vals) {
  n <- 500
  probs=seq(0, 0.95, by=0.05)
  mythresh <- wadsthresh <- norththresh <- numeric(n)
  myquantile <- wadsquantile <- northquantile <- numeric(n)
  myscale <- wadsscale <- northscale <- numeric(n)
  myshape <- wadsshape <- northshape <- numeric(n)
  mylen <- wadslen <- northlen <- numeric(n)
  for(ii in 1:n){
    print(ii)
    set.seed(ii)
    dat1 <- runif(400, 0.5, 1.0)
    dat2 <- rgpd(2000, shape=-0.05, scale=0.5, mu=1.0)
    data <- c(dat1, dat2)
    thresh <- quantile(data, probs, names=F)
    
    #EQD method
    myres <- eqd(data, thresh=thresh, k=k)
    mythresh[ii] <- myres$thresh
    myquantile[ii] <- probs[thresh==mythresh[ii]]
    myscale[ii] <- myres$par[1]
    myshape[ii] <- myres$par[2]
    mylen[ii] <- myres$num
  }
  
  
  eqdthr_case3 <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)

  filename <- paste("output/threshold_selection/additional_sensitivity/eqd_case3_k",k,".csv",sep="")
  write.csv(eqdthr_case3, filename)
}

for(k in k_vals){
  var_name <- paste("eqd_case4_k",k, sep = "")
  filename <- paste("output/threshold_selection/additional_sensitivity/eqd_case4_k",k,".csv",sep="")
  assign(var_name, read.csv(filename, row.names = 1, header = T))
}

rmse(eqd_case4_k10$thr, 1)
rmse(eqd_case4_k20$thr, 1)
rmse(eqd_case4_k50$thr, 1)
rmse(eqd_case4_k80$thr, 1)

bias(eqd_case4_k10$thr, 1)
bias(eqd_case4_k20$thr, 1)
bias(eqd_case4_k50$thr, 1)
bias(eqd_case4_k80$thr, 1)

var(eqd_case4_k10$thr)
var(eqd_case4_k20$thr)
var(eqd_case4_k50$thr)
var(eqd_case4_k80$thr)

#Quantile estimates

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
eqd_list <- list(list(eqd_case1_k10, eqd_case1_k20, eqd_case1_k50, eqd_case1_k80), list(eqd_case2_k10, eqd_case2_k20, eqd_case2_k50, eqd_case2_k80), list(eqd_case3_k10, eqd_case3_k20, eqd_case3_k50, eqd_case3_k80), list(eqd_case4_k10, eqd_case4_k20, eqd_case4_k50, eqd_case4_k80))
RMSE_case1234 <- vector('list',4)
sample_sizes <- c(1200, 480, 2400, 1000)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
for(case in 1:4){
  RMSE_case1234[[case]] <- vector('list',4)
  n <- sample_sizes[case]
  p_seq <- c(1/n,1/(10*n), 1/(100*n))
  par <- true_pars[case,]
  for(k in 1:length(k_vals)){
    eqd_df <- eqd_list[[case]][[k]]
    
    eqd_est_quants <- sapply(p_seq, estimated_quantile, df=eqd_df, n=n)
    
    if(case != 4){
      true_quants <- sapply(p_seq, case123_true_quant, par=par, u=1)
    }
    else{
      true_quants <- sapply(p_seq, case4_true_quant, par=par, p_1= 0.7206618, u=1)
    }
    
    eqd_rmse_bias_var_quants <- matrix(NA, nrow=3, ncol=3)
    for(i in 1:3){
      eqd_rmse_bias_var_quants[i,1] <- rmse(eqd_est_quants[,i], true_quants[i])
      eqd_rmse_bias_var_quants[i,2] <- bias(eqd_est_quants[,i], true_quants[i])
      eqd_rmse_bias_var_quants[i,3] <- var(eqd_est_quants[,i])
    }
    RMSE_case1234[[case]][[k]] <- data.frame(j=c(0,1,2), RMSE=eqd_rmse_bias_var_quants[,1], BIAS=eqd_rmse_bias_var_quants[,2], VAR=eqd_rmse_bias_var_quants[,3])
  }
}

print(RMSE_case1234)

eqd_list[[1]][[1]]
