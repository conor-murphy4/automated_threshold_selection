###########################
# Function definitions ----
###########################

source("helper_functions.R") # for GPD density and quantiles

density_discontinuous <- function(x, uniform_lower, uniform_upper, proportion_uniform, gpd_scale, gpd_shape){

  is_undefined <- which(x < uniform_lower)
  is_uniform <- which((x > uniform_lower) & (x < uniform_upper))
  is_gpd <- which(x > uniform_upper)

  density <- rep(NA_real_, length(x))
  density[is_undefined] <- 0
  density[is_uniform] <- proportion_uniform / (uniform_upper - uniform_lower)
  density[is_gpd] <- dgpd(x[is_gpd], gpd_shape, gpd_scale, mu = uniform_upper)
  density[is_gpd] <- density[is_gpd] * (1 - proportion_uniform)
  return(density)
}

density_continuous <- function(x, alpha = 1, beta = 0.5){
  gpd_scale = 0.5
  gpd_shape = 0.1
  
  # probability density up to proportionality (integral less than 1)
  prop_to_density <- function(x, a = alpha, b = beta){
    dgpd(x, shape = 0.1, scale = 0.5) * pbeta(x, a, b)
  } 

  # Area correction to get valid pdf 
    AUC_0_to_1 <- integrate(prop_to_density, 0, 1)$value
    AUC_1_to_inf <- 1 - pgpd(q = 1,shape = gpd_shape, scale = gpd_scale)
    total_AUC <- AUC_0_to_1 + AUC_1_to_inf

  ## Correct density values to ensure unit integrable
    area_correction_factor <- 1 / total_AUC
    density <- prop_to_density(x) * area_correction_factor
  
  return(density)
}

density_continuous_2 <- function(x, alpha = 1, beta = 0.5){
  gpd_scale = 0.5
  gpd_shape = 0.1
  
  is_undefined <- which(x < 0)
  is_bgpd <- which((x > 0) & (x < 1))
  is_gpd <- which(x > 1)

  # probability density up to proportionality (integral less than 1)
  prop_to_density <- function(x, a = alpha, b = beta){
    dgpd(x, shape = gpd_shape, scale = gpd_scale) * pbeta(x, a, b)
  } 
  
  # Calculate area under the curve between 0 and 1
  q <- integrate(f = prop_to_density, lower = 0, upper = 1)$value

  # scale conditional distribution  
  density <- prop_to_density(x)
  density[is_gpd] <- (1 - q) * dgpd(x[is_gpd] - 1, shape = 0.1, scale = 0.6, mu = 0)

  return(density)
}


# sanity check against sample proportions
integrate(density_continuous, lower = 0, upper = 1)
integrate(density_continuous, lower = 0, upper = 100)

integrate(density_continuous_2, lower = 0, upper = 1)
integrate(density_continuous_2, lower = 0, upper = 100)

#############################################
# Evaluate probability density functions ----
#############################################

t <- seq(0, 4, length.out = 5001)

case_1_density <- density_discontinuous(
  x = t,
  uniform_lower = 0.5,
  uniform_upper = 1.0,
  proportion_uniform = 1/6,
  gpd_scale = 0.5,
  gpd_shape = 0.1
)

case_3_density <- density_discontinuous(
  x = t,
  uniform_lower = 0.5,
  uniform_upper = 1.0,
  proportion_uniform = 1/6,
  gpd_scale = 0.5,
  gpd_shape = -0.05
)

case_4_density_1 <- density_continuous(t)
case_4_density_1a <- density_continuous(t, beta = 2)
case_4_density_2 <- density_continuous_2(t)
###################
# Create Plots ----
###################

pdf("simulation_study/case_study_plots.pdf", width = 7, height = 5)
par(oma = c(0.1,0.1,0.1,0.1), mar=c(4.5, 4.5, 1, 1))
plot(
  x = t,
  y = case_1_density,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2),
  #main = "Cases 1 and 2",
  )
abline(v = 1, lty = 2, lwd = 2)

plot(
  x = t,
  y = case_3_density,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2),
  #main = "Case 3",
  )
abline(v = 1, lty = 2, lwd = 2)

plot(
  x = t,
  y = case_4_density_1,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2),
  #main = "Case 4 - ZV method beta = 0.5",
  )
abline(v = 1, lty = 2, lwd = 2)

plot(
  x = t,
  y = case_4_density_2,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2),
  #main = "Case 4 - JT method beta = 0.5",
  )
abline(v = 1, lty = 2, lwd = 2)

plot(
  x = t,
  y = case_4_density_1a,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2),
  #main = "Case 4 - ZV method beta = 2",
  )
abline(v = 1, lty = 2, lwd = 2)
dev.off()
