EQD_thr_case4 <- readRDS("simulation_study/cases_1-4/simulation_results/estimated_thresholds/EQD_thr_case4.rds")

num_samples <- 500
num_boot <- 200
sample_length <- 1000
return_level_probs <- c(1/sample_length,1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
EQD_return_level_ests <- matrix(nrow = num_boot, ncol = length(return_level_probs))
return_level_boot_estimates <- rep(list(EQD_return_level_ests), num_samples)
for(ii in 1:num_samples){
  EQD_thresh <- EQD_thr_case4$thr[ii]
  EQD_num_excess <- EQD_thr_case4$len[ii]
  EQD_scale <- EQD_thr_case4$scale[ii]
  EQD_shape <- EQD_thr_case4$shape[ii]
  EQD_rate <- EQD_num_excess/sample_length
  set.seed(ii)
  for(jj in 1:num_boot){
    excess_b <- rgpd(EQD_num_excess, scale= EQD_scale, shape=EQD_shape)
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    EQD_return_level_ests[jj,] <- EQD_thresh + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/EQD_rate)^(-opt_b$par[2])-1)
  }
  return_level_boot_estimates[[ii]] <- EQD_return_level_ests
}


#Censored case, p exceedance prob, p_1 P(<1), par (sig, xi), u true
cens_quant <- function(p, par, p_1, u){
  x_p <- ((par[1]+par[2]*(u))/par[2])*((p/(1-p_1))^(-par[2])-1) + u
  return(x_p)  
}

true_quantiles <- cens_quant(p=return_level_probs, par=c(0.5,0.1), p_1= 0.5290265, u=1.0)

return_level_boot_estimates <- return_level_boot_ests_Alg2
covered <- matrix(0,nrow=num_samples, ncol = length(return_level_probs))
for (ii in 1:num_samples) {
  EQD_return_level_ests <- return_level_boot_estimates[[ii]]
  return_level_boot_ests_sorted <- apply(EQD_return_level_ests, 2, sort)
  lower_limits <- return_level_boot_ests_sorted[0.1*(num_boot)^2,]
  upper_limits <- return_level_boot_ests_sorted[0.9*(num_boot)^2,]
  covered[ii,(lower_limits <= true_quantiles & upper_limits >= true_quantiles)] <- 1
}

(coverage_probs <- colSums(covered)/num_samples)
  
