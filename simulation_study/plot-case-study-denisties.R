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

density_continuous <- function(x){

  is_undefined <- which(x < 0)
  is_bgpd <- which((x > 0) & (x < 1))
  is_gpd <- which(x > 1)

  density <- rep(NA_real_, length(x))
  density[is_undefined] <- 0
  density[is_gpd] <- dgpd(x[is_gpd], shape = 0.1, scale = 0.5, mu = 0)
  density[is_bgpd] <- dgpd(x[is_bgpd], shape = 0.1, scale = 0.5) *
                      pbeta(x[is_bgpd], shape1 = 1, shape2 = 1/0.5)

  # Area correction to get valid pdf
  ## Calculate area under the curve between 0 and 1
    n_eval_points = 10001
    eval_points <- ppoints(n_eval_points)
    eval_values <- dgpd(eval_points, shape = 0.1, scale = 0.5) *
      pbeta(eval_points, shape1 = 1, shape2 = 1/0.5)
    AUC_0_to_1 <- mean(eval_values)

  ## Calculate total AUC
    AUC_1_to_inf <- 1 - pgpd(q = 1,shape = 0.1, scale = 0.5)
    total_AUC <- AUC_0_to_1 + AUC_1_to_inf

  ## Correct density values to ensure unit integrabile
    area_correction_factor <- 1 / total_AUC
    density <- density * area_correction_factor

  return(density)
}

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
  proportion_uniform = 1/5,
  gpd_scale = 0.5,
  gpd_shape = -0.05
)

case_4_density <- density_continuous(t)

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
  ylim = c(0,2))
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
  ylim = c(0,2))
abline(v = 1, lty = 2, lwd = 2)

plot(
  x = t,
  y = case_4_density,
  type = "l",
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 1.75,
  xlab = "x",
  ylab = "probability density",
  bty = "n",
  ylim = c(0,2))
abline(v = 1, lty = 2, lwd = 2)
dev.off()
