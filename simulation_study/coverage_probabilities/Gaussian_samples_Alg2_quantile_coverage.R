
source("thr_selection_final.R")

num_samples <- 500
num_boot <- 200
sample_length <- 2000
probs=seq(0.5, 0.95, by=0.05)
my_thresh_boot <- matrix(nrow = num_samples, ncol=num_boot)
return_level_probs <- c(1/sample_length,1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
my_return_level_ests <- matrix(nrow = num_boot^2, ncol = length(return_level_probs))
return_level_boot_estimates <- rep(list(my_return_level_ests), num_samples)
for(ii in 1:num_samples){
  set.seed(ii)
  data <- rnorm(sample_length)
  for(kk in 1:num_boot){
    data_b <- sample(data, sample_length, replace = TRUE)
    thresh <- quantile(data_b, probs=probs, names=FALSE)
    myres <- thr_dist(data_b, thresh=thresh)
    my_thresh_boot[ii,kk] <- myres$thresh
    my_num_excess_b <- myres$num
    my_scale_b <- myres$par[1]
    my_shape_b <- myres$par[2]
    my_rate_b <- my_num_excess_b/sample_length
    for(jj in 1:num_boot){
      excess_b <- rgpd(my_num_excess_b, scale= my_scale_b, shape=my_shape_b)
      opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
      my_return_level_ests[num_boot*(kk-1) + jj,] <- my_thresh_boot[ii,kk] + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/my_rate_b)^(-opt_b$par[2])-1)
    }
  }
  return_level_boot_estimates[[ii]] <- my_return_level_ests
}

saveRDS(return_level_boot_estimates, "return_level_boot_ests_Alg2.rds")

true_quantiles <- qnorm(1-return_level_probs)
covered <- matrix(0,nrow=num_samples, ncol = length(return_level_probs))
for (ii in 1:num_samples) {
  my_return_level_ests <- return_level_boot_estimates[[ii]]
  return_level_boot_ests_sorted <- apply(my_return_level_ests, 2, sort)
  lower_limits <- return_level_boot_ests_sorted[0.025*(num_boot^2),]
  upper_limits <- return_level_boot_ests_sorted[0.975*(num_boot^2),]
  covered[ii,(lower_limits <= true_quantiles & upper_limits >= true_quantiles)] <- 1
}

(coverage_probs <- colSums(covered)/num_samples)
