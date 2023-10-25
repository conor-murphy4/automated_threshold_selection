
source("thresh_qq_metric.R")

data_matrix <- readRDS("simulation_study/cases_1-4/simulation_results/data/data_sim_study_Case_4.rds")

#Note: this will take some time to run!

num_samples <- 500
num_boot <- 200
sample_length <- 1000
probs=seq(0, 0.95, by=0.05)
EQD_thresh_boot <- matrix(nrow = num_samples, ncol=num_boot)
return_level_probs <- c(1/sample_length,1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
EQD_return_level_ests <- matrix(nrow = num_boot^2, ncol = length(return_level_probs))
return_level_boot_estimates <- rep(list(EQD_return_level_ests), num_samples)
for(ii in 1:num_samples){
  data_matrix[,ii] <- data
  for(kk in 1:num_boot){
    data_b <- sample(data, sample_length, replace = TRUE)
    thresh <- quantile(data_b, probs=probs, names=FALSE)
    EQDres <- thresh_qq_metric(data_b, thresh=thresh)
    EQD_thresh_boot[ii,kk] <- EQDres$thresh
    EQD_num_excess_b <- EQDres$num
    EQD_scale_b <- EQDres$par[1]
    EQD_shape_b <- EQDres$par[2]
    EQD_rate_b <- EQD_num_excess_b/sample_length
    for(jj in 1:num_boot){
      excess_b <- rgpd(EQD_num_excess_b, scale= EQD_scale_b, shape=EQD_shape_b)
      opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
      EQD_return_level_ests[num_boot*(kk-1) + jj,] <- EQD_thresh_boot[ii,kk] + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/EQD_rate_b)^(-opt_b$par[2])-1)
    }
  }
  return_level_boot_estimates[[ii]] <- EQD_return_level_ests
}

true_quantiles <- qnorm(1-return_level_probs)
covered <- matrix(0,nrow=num_samples, ncol = length(return_level_probs))
for (ii in 1:num_samples) {
  EQD_return_level_ests <- return_level_boot_estimates[[ii]]
  return_level_boot_ests_sorted <- apply(EQD_return_level_ests, 2, sort)
  lower_limits <- return_level_boot_ests_sorted[0.025*(num_boot^2),]
  upper_limits <- return_level_boot_ests_sorted[0.975*(num_boot^2),]
  covered[ii,(lower_limits <= true_quantiles & upper_limits >= true_quantiles)] <- 1
}

(coverage_probs <- colSums(covered)/num_samples)
