mythrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Replicated with quantiles/mythrG.rds")
head(mythrG)
num_samples <- 200
num_boot <- 200
sample_length <- 2000
return_level_probs <- c(1/sample_length,1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
my_return_level_ests <- matrix(nrow = num_boot, ncol = length(return_level_probs))
return_level_boot_estimates <- rep(list(my_return_level_ests), num_samples)
for(ii in 1:num_samples){
  my_thresh <- mythrG$thr[ii]
  my_num_excess <- mythrG$len[ii]
  my_scale <- mythrG$scale[ii]
  my_shape <- mythrG$shape[ii]
  my_rate <- my_num_excess/sample_length
  set.seed(ii)
  for(jj in 1:num_boot){
    excess_b <- rgpd(my_num_excess, scale= my_scale, shape=my_shape)
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    my_return_level_ests[jj,] <- my_thresh + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/my_rate)^(-opt_b$par[2])-1)
  }
  return_level_boot_estimates[[ii]] <- my_return_level_ests
}

saveRDS(return_level_boot_estimates, "return_level_boot_ests_Alg1.rds")

true_quantiles <- qnorm(1-return_level_probs)
covered <- matrix(0,nrow=num_samples, ncol = length(return_level_probs))
for (ii in 1:num_samples) {
  my_return_level_ests <- return_level_boot_estimates[[ii]]
  return_level_boot_ests_sorted <- apply(my_return_level_ests, 2, sort)
  lower_limits <- return_level_boot_ests_sorted[0.025*num_boot,]
  upper_limits <- return_level_boot_ests_sorted[0.975*num_boot,]
  covered[ii,(lower_limits <= true_quantiles & upper_limits >= true_quantiles)] <- 1
}

(coverage_probs <- colSums(covered)/num_samples)
  
