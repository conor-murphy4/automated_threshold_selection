
library(evir)
source("src/eqd.R")

data("nidd.thresh")

thresholds <- quantile(nidd.thresh, seq(0,0.93, by=0.01))

set.seed(11111) 
EQD_thr <- eqd(nidd.thresh, thresholds, k=200)

#QQplot with tolerance bounds incorporating parameter uncertainty
u_hat <- EQD_thr$thresh 
excess <- nidd.thresh[nidd.thresh>u_hat] - u_hat 
n_u <- length(excess) 
probs <- c(1:n_u)/(n_u+1) 

q_sample <- quantile(excess, probs = probs, names=F)

model_fit <- optim(GPD_LL, par=c(mean(excess), 0.1), z=excess, control = list(fnscale=-1))
q_model <- qgpd(probs, shape = EQD_thr$par[2], scale=EQD_thr$par[1])

m_boot <- 200 
boot_q_model_estimates <- matrix(NA, nrow = m_boot, ncol=length(probs)) 

set.seed(12345)
for(i in 1:m_boot){
  excess_boot <- rgpd(n_u, scale = EQD_thr$par[1], shape= EQD_thr$par[2]) 
  boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) 
  boot_q_model_estimates[i,] <- qgpd(probs, shape = boot_model_fit$par[2], scale=boot_model_fit$par[1]) 
}

#95% tolerance bounds for model quantiles
q_l <- apply(boot_q_model_estimates, 2 , quantile, probs=0.025) 
q_u <- apply(boot_q_model_estimates, 2 , quantile, probs=0.975) 

#Return level plot with parameter uncertainty and additional threshold uncertainty

t <- c(1,2,5,10,25,50,100,500,1000) #return periods

point_estimates <- u_hat + (EQD_thr$par[1]/EQD_thr$par[2])*(((t*n_u)/35)^(EQD_thr$par[2])-1) #return level point estimates using fitted GPD parameters

boot_return_ests <- matrix(NA, nrow = m_boot, ncol=length(t))

#Obtaining bootstrap estimates of return levels incorporating uncertainty in parameter estimates
set.seed(12345)
for(i in 1:m_boot){
  excess_boot <- rgpd(n_u, scale = EQD_thr$par[1], shape= EQD_thr$par[2]) 
  boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) 
  boot_return_ests[i,] <- u_hat + (boot_model_fit$par[1]/boot_model_fit$par[2])*(((t*n_u)/35)^(boot_model_fit$par[2])-1) 
}

n <- length(nidd.thresh)
boot_return_ests_thr_uncertainty <- matrix(NA, nrow = (m_boot^2), ncol=length(t))
boot_thr <- numeric(m_boot)

#Obtaining bootstrap estimates of return levels incorporating uncertainty in parameter estimates and threshold
#NOTE: This section of code will take some time to run!
set.seed(12345)
for(i in 1:m_boot){
  data_boot <- sample(nidd.thresh,n, replace=TRUE) 
  thresholds_boot <- quantile(data_boot, seq(0,0.99, by=0.01)) 
  EQD_thr_boot <- eqd(data_boot, thresh = thresholds_boot, k=200) 
  boot_thr[i] <- EQD_thr_boot$thresh 
  n_excess_boot <- length(data_boot[data_boot>EQD_thr_boot$thresh])  
  for(jj in 1:m_boot){
    excess_boot <- rgpd(n_excess_boot, scale = EQD_thr_boot$par[1], shape= EQD_thr_boot$par[2])
    boot_model_fit <- optim(GPD_LL, par=c(mean(excess_boot), 0.1), z=excess_boot, control = list(fnscale=-1)) 
    boot_return_ests_thr_uncertainty[m_boot*(i-1) + jj,] <- boot_thr[i] + (boot_model_fit$par[1]/boot_model_fit$par[2])*(((t*n_excess_boot)/35)^(boot_model_fit$par[2])-1) 
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
