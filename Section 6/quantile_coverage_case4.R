source("src/helper_functions.R")
source("src/eqd.R")

data_matrix <- read.csv("data/Case_4.csv", row.names = 1, header=T)

eqd_case4 <- read.csv("output/threshold_selection/eqd_case4.csv", row.names = 1, header = T)

# Algorithm 1: Parameter uncertainty with known threshold
num_samples <- 500
num_boot <- 200
sample_length <- 1000
return_level_probs <- c(1/sample_length,1/(3*sample_length),1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
return_level_ests_A1 <- matrix(nrow = num_boot, ncol = length(return_level_probs))
A1_estimates_list <- rep(list(return_level_ests_A1), num_samples)
for(ii in 1:num_samples){
  eqd_thresh <- eqd_case4$thr[ii]
  num_excess <- eqd_case4$len[ii]
  eqd_scale <- eqd_case4$scale[ii]
  eqd_shape <- eqd_case4$shape[ii]
  exceedance_rate <- num_excess/sample_length
  set.seed(ii)
  for(jj in 1:num_boot){
    excess_b <- rgpd(num_excess, scale= eqd_scale, shape=eqd_shape)
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    return_level_ests_A1[jj,] <- eqd_thresh + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/exceedance_rate)^(-opt_b$par[2])-1)
  }
  A1_estimates_list[[ii]] <- return_level_ests_A1
}

saveRDS(A1_estimates_list, "output/coverage/Alg1_quantile_estimates_case4.rds")

# Algorithm 1b: Parameter uncertainty with known threshold (including rate uncertainty)
num_samples <- 500
num_boot <- 200
sample_length <- 1000
return_level_probs <- c(1/sample_length,1/(3*sample_length),1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
return_level_ests_A1b <- matrix(nrow = num_boot, ncol = length(return_level_probs))
A1b_estimates_list <- rep(list(return_level_ests_A1b), num_samples)
for(ii in 1:num_samples){
  eqd_thresh <- eqd_case4$thr[ii]
  eqd_num_excess <- eqd_case4$len[ii]
  eqd_scale <- eqd_case4$scale[ii]
  eqd_shape <- eqd_case4$shape[ii]
  exceedance_rate <- eqd_num_excess/sample_length
  set.seed(ii)
  for(jj in 1:num_boot){
    num_excess <- rbinom(1, sample_length, exceedance_rate)
    excess_b <- rgpd(num_excess, scale= eqd_scale, shape=eqd_shape)
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    return_level_ests_A1b[jj,] <- eqd_thresh + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/exceedance_rate)^(-opt_b$par[2])-1)
  }
  A1b_estimates_list[[ii]] <- return_level_ests_A1b
}

saveRDS(A1b_estimates_list, "output/coverage/Alg1b_quantile_estimates_case4.rds")

# Algorithm 2: Parameter uncertainty for unknown threshold
probs=seq(0, 0.95, by=0.05)
eqd_thresh_boot <- matrix(nrow = num_samples, ncol=num_boot)
return_level_probs <- c(1/sample_length,1/(3*sample_length),1/(5*sample_length),1/(10*sample_length), 1/(25*sample_length), 1/(50*sample_length), 1/(100*sample_length), 1/(200*sample_length), 1/(500*sample_length))
return_level_ests_A2 <- matrix(nrow = num_boot^2, ncol = length(return_level_probs))
A2_estimates_list <- rep(list(return_level_ests_A2), num_samples)
for(ii in 1:num_samples){
  set.seed(ii)
  data <- data_matrix[,ii]
  for(kk in 1:num_boot){
    data_b <- sample(data, sample_length, replace = TRUE)
    thresh <- quantile(data_b, probs=probs, names=FALSE)
    eqd_selection <- eqd(data_b, thresh=thresh)
    eqd_thresh_boot[ii,kk] <- eqd_selection$thresh
    num_excess_b <- eqd_selection$num_excess
    scale_b <- eqd_selection$par[1]
    shape_b <- eqd_selection$par[2]
    exceedance_rate_b <- num_excess_b/sample_length
    for(jj in 1:num_boot){
      excess_b <- rgpd(num_excess_b, scale= scale_b, shape=shape_b)
      opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
      return_level_ests_A2[num_boot*(kk-1) + jj,] <- eqd_thresh_boot[ii,kk] + (opt_b$par[1]/opt_b$par[2])*((return_level_probs/exceedance_rate_b)^(-opt_b$par[2])-1)
    }
  }
  A2_estimates_list[[ii]] <- return_level_ests_A2
}

saveRDS(A2_estimates_list, "output/coverage/Alg2_quantile_estimates_case4.rds")


#Case 4, p exceedance prob, p_1 P(<1), par (sig, xi), u true
true_quant <- function(p, par, p_1, u){
  x_p <- ((par[1]+par[2]*(u))/par[2])*((p/(1-p_1))^(-par[2])-1) + u
  return(x_p)  
}

A1_estimates_list <- readRDS("output/coverage/Alg1_quantile_estimates_case4.rds")
A1b_estimates_list <- readRDS("output/coverage/Alg1b_quantile_estimates_case4.rds")
A2_estimates_list <- readRDS("output/coverage/Alg2_quantile_estimates_case4.rds")

true_quantiles <- true_quant(p=return_level_probs, par=c(0.5,0.1), p_1= 0.7206618, u=1.0)

coverage_table <- vector('list', 3)

sig_levels <- c(0.5, 0.2, 0.05) 

for(jj in 1:3){
  sig_level <- sig_levels[jj]
  
  covered_A1 <- covered_A1b <- covered_A2 <- CI_widths_A1 <- CI_widths_A2 <- matrix(0,nrow=num_samples, ncol = length(return_level_probs))
  for (ii in 1:num_samples) {
    return_level_ests_A1 <- A1_estimates_list[[ii]]
    return_level_ests_A1b <- A1b_estimates_list[[ii]]
    return_level_ests_A2 <- A2_estimates_list[[ii]]
    lower_limits_A1 <- apply(return_level_ests_A1, 2, quantile, prob = (sig_level/2))
    lower_limits_A1b <- apply(return_level_ests_A1b, 2, quantile, prob = (sig_level/2))
    lower_limits_A2 <- apply(return_level_ests_A2, 2, quantile, prob = (sig_level/2))
    upper_limits_A1 <- apply(return_level_ests_A1, 2, quantile, prob = (1-sig_level/2))
    upper_limits_A1b <- apply(return_level_ests_A1b, 2, quantile, prob = (1-sig_level/2))
    upper_limits_A2 <- apply(return_level_ests_A2, 2, quantile, prob = (1-sig_level/2))
    CI_widths_A1[ii,] <- upper_limits_A1 - lower_limits_A1
    CI_widths_A2[ii,] <- upper_limits_A2 - lower_limits_A2
    covered_A1[ii,(lower_limits_A1 <= true_quantiles & upper_limits_A1 >= true_quantiles)] <- 1
    covered_A1b[ii,(lower_limits_A1b <= true_quantiles & upper_limits_A1b >= true_quantiles)] <- 1
    covered_A2[ii,(lower_limits_A2 <= true_quantiles & upper_limits_A2 >= true_quantiles)] <- 1
  }
  
  coverage_probs_A1 <- colSums(covered_A1)/num_samples
  coverage_probs_A1b <- colSums(covered_A1b)/num_samples
  coverage_probs_A2 <- colSums(covered_A2)/num_samples
  CI_ratios <- round(colSums(CI_widths_A2/CI_widths_A1)/num_samples, 3)
  coverage_table[[jj]] <- as.data.frame(matrix(c(coverage_probs_A1, coverage_probs_A1b, coverage_probs_A2, CI_ratios), nrow=4, ncol=length(return_level_probs), byrow=T))
  rownames(coverage_table[[jj]]) <- c("Alg 1", "Alg 1b", "Alg 2", "CI ratio")
}

names(coverage_table) <- c("50%", "80%", "95%")

print(coverage_table)

saveRDS(coverage_table, "output/tables/Table_S.22_coverage_probs_estimated_quantiles_case4.rds")



