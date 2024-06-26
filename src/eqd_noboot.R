# Load helper functions for working with the generalised Pareto distribution
source('src/helper_functions.R')

#-----------------------------------------------------------------------

#' Threshold selection method for univariate extremes
#'
#' 'eqd_noboot' is an adjustment of the 'eqd' function which selects a constant threshold above which the data can be most closely modelled by a Generalised Pareto distribution but this version selects a threshold based solely on the fit to the observed sample, i.e. without the use of boostrapping.
#'
#' @author Conor Murphy
#'
#' @param data A numeric vector.
#' @param thresh A numeric vector of proposed thresholds to test.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.


eqd_noboot <- function(data, thresh, m = 500){

  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.numeric(thresh)) stop("u to be tested needs to be a vector")
  # if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")

  meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
  for (i in 1:length(thresh)) {
    u <- thresh[i]
    excess <- data[data > u] - u
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      gpd.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
      xis[i] <- gpd.fit$par[[2]]
      sigmas[i] <- gpd.fit$par[[1]]
      quants <- qgpd((1:m) / (m + 1), scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
      meandistances[i] <- (1 / m) * sum(abs(quantile(excess, probs = (1:m) / (m+1)) - quants))
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold <- thresh[chosen_index]
  xi <- xis[chosen_index]
  sigma <- sigmas[chosen_index]
  len <- num_excess[chosen_index]
  result <- list(thresh = chosen_threshold, par = c(sigma,xi), num_excess = len, dists = meandistances)
  return(result)
}

