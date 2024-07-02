###########################
# Function definitions ----
###########################

source("helper_functions.R") # for GPD density and quantiles

density_discontinuous <- function(x,
    uniform_lower = 0.5,
    uniform_upper = 1.0, 
    proportion_uniform = 1/6,
    gpd_scale = 0.5,
    gpd_shape = 0.1){

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

density_continuous <- function(x, alpha = 1, beta = 2){
  gpd_scale = 0.5
  gpd_shape = 0.1
  
  # probability density up to proportionality (integral less than 1)
  prop_to_density <- function(x, a = alpha, b = beta){
    dgpd(x, shape = 0.1, scale = 0.5) * pbeta(x, a, b)
  } 

  # Calculate area under prop_to_density() 
    AUC_0_to_1 <- integrate(prop_to_density, 0, 1)$value
    AUC_1_to_inf <- 1 - pgpd(q = 1,shape = gpd_shape, scale = gpd_scale)
    total_AUC <- AUC_0_to_1 + AUC_1_to_inf

  ## Correct density values to ensure unit integrable
    area_correction_factor <- 1 / total_AUC
    density <- prop_to_density(x) * area_correction_factor
  
  return(density)
}


#############################################
# Evaluate probability density functions ----
#############################################

t <- seq(0, 4, length.out = 5001)

densities_by_case <- list(
  case_1 = density_discontinuous(x = t, gpd_shape = 0.1),
  case_2 = density_discontinuous(x = t, gpd_shape = 0.1),   # same as case 1
  case_3 = density_discontinuous(x = t, gpd_shape = -0.05), 
  case_4 = density_continuous(x = t)
)

###################
# Create Plots ----
###################

pdf("simulation_study/case_study_denisty_plots.pdf", width = 7, height = 5)
par(oma = c(0.1,0.1,0.1,0.1), mar=c(4.5, 4.5, 1, 1))
for (density_values in densities_by_case){
  plot(
    x = t,
    y = density_values,
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
}
dev.off()


## Tests ----------------------------------------------------------------

# Sanity check against sample proportion below threshold, (approx 0.721)

integrate(density_continuous, lower = 0, upper = 1)
integrate(density_continuous, lower = 0, upper = 100)

# Extend over all quantiles

## Simulate large set of Xs
set.seed(4321)
n_sim <- 1e6
Y <- threshold:::rgpd(n = n_sim, shape = 0.1, scale = 0.5)
B <- rbeta(n = n_sim, shape1 = 1, shape2 = 2)
X <- Y[Y>B]
tau <- mean(X < 1)

## Define model and sample cdfs

model_cdf   <- function(to){integrate(density_continuous, 0, to)$value}
sample_cdf  <- function(to){mean(X <= to)}

## Compare CDFs directly
eval_points <- 8 * ppoints(1001)
cdf_model   <- purrr::map_dbl(eval_points, .f = model_cdf)
cdf_sample  <- purrr::map_dbl(eval_points, .f = sample_cdf)

# Overlay CDFs
pdf("simulation_study/cdf_check.pdf", width = 7, height = 5)
plot(
  x = eval_points,
  y = cdf_model, 
  xlab = "x",
  ylab = "CDF", 
  type = "l",
  lwd = 3,
  col = "darkorange")
lines(eval_points, cdf_sample, col = "black", lty = 2)
legend(
  "bottomright",
  legend = c("Model", "Empirical"),
  lwd = 2,
  lty = c(1,2),
  col = c("darkorange", "black")
)
dev.off()

## QQ plot
pdf("simulation_study/qq_plot_check.pdf", width = 5, height = 5)
plot(
  x = cdf_sample, 
  y = cdf_model, 
  col = "darkorange", 
  type = "l", 
  lwd = "3",
  xlab = "sample quantile", 
  ylab = "model quantile")
abline(a = 0 , b = 1, lwd = 1, col = "black")
legend(
  "bottomright",
  legend = c("Model", "Empirical"),
  lwd = 2,
  lty = c(1,2),
  col = c("darkorange", "black")
)
dev.off()
