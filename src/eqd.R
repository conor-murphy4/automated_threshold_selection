# Load helper functions for working with the generalised Pareto distribution
source('src/helper_functions.R')

#-----------------------------------------------------------------------

#' Threshold selection method for univariate extremes
#'
#' 'eqd' selects a constant threshold above which the data can be most closely modelled by a Generalised Pareto distribution.
#'
#' @author Conor Murphy
#'
#' @param data A numeric vector.
#' @param thresh A numeric vector of proposed thresholds to test.
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.
#'
#' @examples
#' set.seed(12345)
#' data_test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' thresholds1 <- quantile(data_test1,seq(0,0.95,by=0.05))
#' (example1 <- eqd(data_test1, thresh = thresholds1))
#'
#' set.seed(11111)
#' 
#' 
#' test2_data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
#' u <- 1
#' cens_thr <- u*rbeta(length(test2_data_all),1,2)
#' data_above <- sample(test2_data_all[test2_data_all > u], 279, replace=FALSE)
#' data_below <- sample(test2_data_all[test2_data_all > cens_thr & test2_data_all <= u], 721, replace = FALSE) 
#' data_test2 <- c(data_below, data_above)
#' thresholds2 <- quantile(data_test2,seq(0, 0.95, by=0.05))
#' (example2 <- eqd(data_test2,thresh = thresholds2))


eqd <- function(data, thresh, k = 100, m = 500){

  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.numeric(thresh)) stop("u to be tested needs to be a vector")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")

  meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
  for (i in 1:length(thresh)) {
    u <- thresh[i]
    excess <- data[data > u] - u
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      for (j in 1:k) {
        X <- sample(excess, num_excess[i], replace = TRUE)
        mle <- mean(X)
        ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1))
        quants <- qgpd((1:m) / (m + 1), scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
        distances[j] <- (1 / m) * sum(abs(quantile(X, probs = (1:m) / (m+1)) - quants))
      }
      meandistances[i] <- mean(distances)
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

