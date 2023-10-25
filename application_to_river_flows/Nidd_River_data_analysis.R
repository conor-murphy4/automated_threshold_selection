#Nidd River data
library(evir)
data("nidd.thresh")

#EQD method
source("thresh_qq_metric.R")

#Wadsworth method
source("simulation_study/JointMLEFunctions.r")

#Northrop method
library(threshr)


# Results for Table 9 -----------------------------------------------------

#Candidate threshold grids
(thresholds <- quantile(nidd.thresh, seq(0,0.99, by=0.01)))
(thresholds <- quantile(nidd.thresh, seq(0,0.98, by=0.01)))
(thresholds <- quantile(nidd.thresh, seq(0,0.9, by=0.01)))
(thresholds <- quantile(nidd.thresh, seq(0,0.8, by=0.01)))
(thresholds <- quantile(nidd.thresh, seq(0,0.8, by=0.2)))
(thresholds <- quantile(nidd.thresh, seq(0,0.9, by=0.3)))
(thresholds <- quantile(nidd.thresh, seq(0,0.75, by=0.25)))
(thresholds <- quantile(nidd.thresh, c(0,0.1,0.4,0.7)))

#Estimated thresholds for each method
#EQD
set.seed(11111) 
(EQD_thr <- thresh_qq_metric(nidd.thresh, thresholds, k=200))
#Wadsworth
set.seed(11111)
(wadsres <- NHPP.diag(nidd.thresh, u=thresholds, plot.out=FALSE, UseQuantiles = FALSE))
#Northrop
set.seed(11111)
(norththr <- summary(ithresh(nidd.thresh, u_vec = thresholds)))


# Results for Figure 2 -----------------------------------------------------

#QQplot with tolerance bounds incorporating parameter uncertainty
u_hat <- EQD_thr$thresh #EQD threshold estimate
excess <- nidd.thresh[nidd.thresh>u_hat] - u_hat #set of excesses of EQD threshold estimate
n_u <- length(excess) #number of excesses
probs <- c(1:n_u)/(n_u+1) #evaluation probabilities

#Sample quantiles
q_sample <- quantile(excess, probs = probs, names=F)

#Model quantile estimates
model_fit <- optim(GPD_LL, par=c(mean(excess), 0.1), z=excess, control = list(fnscale=-1))
q_model <- qgpd(probs, shape = model_fit$par[2], scale=model_fit$par[1])

m <- 200 #Number of bootstraps
boot_q_mod <- matrix(NA, nrow = m, ncol=length(probs)) 

#Obtaining bootstrap estimates of quantiles incorporating uncertainty in parameter estimates
set.seed(12345)
for(i in 1:m){
  excess_boot <- rgpd(n_u, scale = model_fit$par[1], shape= model_fit$par[2]) #simulating GPD samples of excesses of length n_u using fitted parameters
  boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) #Refitting to new GPD sample
  boot_q_mod[i,] <- qgpd(probs, shape = boot_model_fit$par[2], scale=boot_model_fit$par[1]) #Estimating quantiles at same set of evaluation probabilites
}

#95% tolerance bounds for model quantiles
q_l <- apply(boot_q_mod, 2 , quantile, probs=0.025) #lower bound
q_u <- apply(boot_q_mod, 2 , quantile, probs=0.975) #upper bound

#Return level plot with paramater uncertainty and additional threshold uncertainty

t <- c(1,2,5,10,25,50,100,500,1000) #return periods

point_estimates <- u_hat + (model_fit$par[1]/model_fit$par[2])*(((t*n_u)/35)^(model_fit$par[2])-1) #return level point estimates using fitted GPD parameters

boot_return_ests <- matrix(NA, nrow = m, ncol=length(t))

#Obtaining bootstrap estimates of return levels incorporating uncertainty in parameter estimates
set.seed(12345)
for(i in 1:m){
  excess_boot <- rgpd(n_u, scale = model_fit$par[1], shape= model_fit$par[2]) #simulating GPD samples of excesses of length n_u using fitted parameters
  boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) #Refitting to new GPD sample
  boot_return_ests[i,] <- u_hat + (boot_model_fit$par[1]/boot_model_fit$par[2])*(((t*n_u)/35)^(boot_model_fit$par[2])-1) #Estimating return levels at same set of return periods
}

n <- length(nidd.thresh)
boot_return_ests_thr_uncertainty <- matrix(NA, nrow = m*m, ncol=length(t))
boot_thr <- numeric(m)

#Obtaining bootstrap estimates of return levels incorporating uncertainty in parameter estimates and threshold
set.seed(12345)
for(i in 1:m){
  data_boot <- sample(data,n, replace=TRUE) #non-parametric bootstrap to incorporate threshold (and rate of exceedance) uncertainty
  thresholds_boot <- quantile(data_boot, seq(0,0.99, by=0.01)) #candidate threshold grid for bootstrapped sample
  EQD_thr_boot <- thresh_qq_metric(data_boot, thresh = thresholds_boot, k=200) #estimating threshold using EQD for bootstrapped sample
  boot_thr[i] <- EQD_thr_boot$thresh #recording each threshold estimate in vector
  n_excess_boot <- length(data_boot[data_boot>EQD_thr_boot$thresh]) #number of excesses of EQD threshold estimate in bootstrapped sample 
  for(jj in 1:m){
    excess_boot <- rgpd(n_excess_boot, scale = EQD_thr_boot$par[1], shape= EQD_thr_boot$par[2]) #simulating GPD samples of excesses of length n_excess_boot using fitted parameters
    boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) #Refitting to new GPD sample
    boot_return_ests_thr_uncertainty[m*(i-1) + jj,] <- boot_thr[i] + (boot_model_fit$par[1]/boot_model_fit$par[2])*(((t*n_b)/35)^(boot_model_fit$par[2])-1) #Estimating return levels at same set of return periods 
  }
}

#95% confidence intervals incorporating parameter uncertainty
return_l <- apply(boot_return_ests, 2, quantile, 0.025) #lower bound
return_u <- apply(boot_return_ests, 2, quantile, 0.975) #upper bound

#95% confidence intervals incorporating threshold and parameter uncertainty
return_thr_l <- apply(boot_return_ests_thr_uncertainty, 2, quantile, 0.025) #lower bound
return_thr_u <- apply(boot_return_ests_thr_uncertainty, 2, quantile, 0.975) #upper bound

#Plotting window
dev.new(width=9.17, height=5,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

#QQplot as in Figure 2
plot(q_sample, q_model, ylab="Model Quantiles", xlab = "Sample Quantiles", pch=19, asp = 1, cex.lab=1.1) #Plotting point estimates
polygon(c( q_l, sort(q_u, decreasing = TRUE)), c(q_model,  sort(q_model, decreasing = TRUE)), col=rgb(0, 0, 1,0.2), border=NA) #Tolerance bounds
abline(a=0, b=1) #Line of equality

#Return level plot as in Figure 2
plot(point_estimates ~ t, log='x', xlab="Return Period (Years)", type='l', lwd=2, ylab = expression(Return~Level~(m^3/s)), ylim=c(100,2500), cex.lab=1.1) #Plotting point estimates
polygon(c( t, sort(t, decreasing = TRUE)), c(return_l,  sort(return_u, decreasing = TRUE)), col=rgb(1, 0, 0,0.4), border=NA) #95% confidence intervals incorporating parameter uncertainty
polygon(c( t, sort(t, decreasing = TRUE)), c(return_thr_l,  sort(return_thr_u, decreasing = TRUE)), col=rgb(1, 0, 0,0.2), border=NA) #95% confidence intervals incorporating threshold and parameter uncertainty



