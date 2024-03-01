#=====================================================================
# Functions for Generalised Pareto Distribution.
# Added option to use nu parameterisation
# checks that param values are valid
#=====================================================================
# pgpd
# qgpd
# dgpd
# rgpd

# GPD_LL
#=====================================================================
#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
#'
#' @author Zak Varty
#'
#' @param q vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @param skip_checks logical. Speed up evaluation by skipping checks on inputs? (Beware!)
#' @return Probability of the GPD X<=q
#' @importFrom stats pexp
#' @examples
#' pgpd(q = c(-1,1.5,3), shape = 1, scale = 1)
#' pgpd(q = 1.5, shape = c(0,-1), scale = c(0.1,1))
#' @export

pgpd <- function(q, shape, scale = NULL, nu = NULL, mu = 0, skip_checks = FALSE){

  if (!skip_checks) {
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }

    }
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }

    # Ensure q, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(q) > 1) {
      scale <- rep(scale, length(q))
    }
    if (length(shape) == 1 & length(q) > 1) {
      shape <- rep(shape, length(q))
    }
    if (length(mu) == 1 & length(q) > 1) {
      mu <- rep(mu, length(q))
    }
  } else {
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
    }
  }
  #calculate probabilities
  p <- (1 - (1 + (shape * (q - mu))/scale)^(-1/shape))
  #correct probabilities below mu or above upper end point
  p[q < mu] <- 0
  p[(shape < 0) & (q >= (mu - scale/shape))] <- 1

  #correct probabilities where xi = 0
  if (any(abs(shape) < 1e-10)) {
    #ex <- which(shape ==0)
    ex <- which(abs(shape) < 1e-10)
    p[ex] <- pexp(q = q[ex] - mu[ex], rate = 1 / scale[ex])
  }

  return(p)
}

#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
#' probability.
#'
#' @author Zak Varty
#'
#' @param p vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=x
#' @examples
#' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
#' \dontrun{ qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1) }
#' @export
qgpd <- function(p, shape, scale = NULL, nu = NULL, mu = 0){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }

  # Probabilities must all be positive
  if (!all((p >= 0) & (p <= 1))) {
    stop('Probabilities p must be in the range [0,1].')
  }

  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }

  }

  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure p, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(p) > 1) {
    scale <- rep(scale, length(p))
  }
  if (length(shape) == 1 & length(p) > 1) {
    shape <- rep(shape, length(p))
  }
  if (length(mu) == 1 & length(p) > 1) {
    mu <- rep(mu, length(p))
  }

  #calculate quantiles
  q <- mu + (scale/shape) * ((1 - p)^(-shape) - 1)

  #correct quantiles where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    q[ex] <- mu[ex] + stats::qexp(p = p[ex],rate = 1/scale[ex])
  }
  return(q)
}



#' Generalised Pareto Distribution
#'
#' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or x is
#' outside of the domain of the given distribution.
#'
#' @author Zak Varty
#'
#' @param x vector of values as which to evaluate density.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter
#' @param mu  location parameter
#' @param log  locical. Return log
#' @return density of the GPD at x
#' @examples
#' dgpd(x = c(-1,0.5,1,1.9,5),shape = -0.5, scale = 1)
#' @export
#'
dgpd <- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }

  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }

  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure x, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(x) > 1) {
    scale <- rep(scale, length(x))
  }
  if (length(shape) == 1 & length(x) > 1) {
    shape <- rep(shape, length(x))
  }
  if (length(mu) == 1 & length(x) > 1) {
    mu <- rep(mu, length(x))
  }

  if (log == FALSE) {
    out <- (scale^(-1)) * pmax((1 + shape * (x - mu)/scale),0)^((-1/shape) - 1)
    # amend values below threshold
    out[which(x < mu)] <- 0
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- 0
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape < 1e-10))) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1/scale[ex])
    }
  } else {
    out <-  -log(scale) + ((-1/shape) - 1)*log(pmax((1 + shape * (x - mu)/scale),0))
    # amend values below threshold
    out[which(x < mu)] <- -Inf
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- -Inf
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex],log = TRUE)
    }
  }
  return(out)
}

#' Generalised Pareto Distribution
#'
#' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0.
#'
#' @author Zak Varty
#'
#' @param n sample size.
#' @param shape shape parameter (xi).
#' @param scale scale parameter (sigma).
#' @param nu  alternative scale parameter.
#' @param mu  location parameter.
#' @return Random sample from generalised pareto distirbution.
#'
#' @examples
#' rgpd(n = 100, shape = 0, scale = 1:100)
#' @export
rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
  ## Input checks
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure q, scale, shape and mu are of same length.
  if ((length(scale) == 1) & (n > 1)) {
    scale <- rep(scale, n)
  }
  if ((length(shape) == 1) & (n > 1)) {
    shape <- rep(shape, n)
  }
  if ((length(mu) == 1) & (n > 1)) {
    mu <- rep(mu, n)
  }

  #simulate sample
  sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
  #correct sample values where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    sample[ex] <- mu[ex] +
      stats::rexp(n = length(ex),rate = 1/scale[ex])
  }
  return(sample)
}

#rgpd_rd <- function(n, sig, xi, mu, to_nearest, mu_latent = NULL){
#  if(is.null(mu_latent)) mu_latent = mu - 0.5 * to_nearest
#  x <- rgpd(n = n, shape = xi, scale = sig, mu = mu_latent)
#  x <- round_to_nearest(x, to_nearest)
#  return(x)
#}

#' evalue probability mass function of rounded generalised Pareto distribution
#'
#' @author Zak Varty
#'
#' @param x Vector values at which to evaluate mass function
#' @param u Vector of latent threshold values
#' @param sig_u Vector of latent scale parameters (for exceedances of u)
#' @param xi Latent shape parameter
#' @param to_nearest Level of rounding
#'
#' @return pmf evaluated at x. NOTE: does not check validity of x values.
#'
#' @examples
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = 0, sig_u = 1, xi = 0)
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = seq(0, 1.4, by = 0.1), sig_u = 1, xi = 0)
#' # CAUTION:
#' gpd_rd(x = 0.15, to_nearest = 0.1,  u = 0, sig_u = 1, xi = 0)
dgpd_rd <- function(x, u, sig_u, xi, to_nearest){

  # If (Y_i - u_i | Y_i > u_i) ~ GPD(sig_i, xi)
  # then Z_i  = [(Y_i - u_i)/sig_i | Y_i - u_i > 0] ~ GPD(1, xi)

  # range of z values that lead to observing x
  x_low <- pmax(x - to_nearest/2, u)
  z_low <- (x_low - u) / sig_u
  z_high <- (x - u + to_nearest/2)/sig_u

  # calculate probability of z in that range
  p_high <- pgpd(q = z_high, scale = 1, shape = xi, mu = 0)
  p_low  <- pgpd(q = z_low,  scale = 1, shape = xi, mu = 0)
  p <- p_high - p_low

  return(p)
}


#' Generalised Pareto log-likelihood
#'
#' @author Conor Murphy
#'
#' @param par A numeric vector of parameter values of length 2.
#' @param z A numeric vector of excesses of some threshold.
#'
#' @returns A numeric value of the log-likeihood.
#'
#' @examples
#' test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' excess <- test1[test1>1.5] - 1.5
#' GPD_LL(par=c(1,0.4), z=excess)


GPD_LL <- function(par, z){
  sigma <- par[1]
  xi <- par[2]
  if (sigma > 0) {
    if (abs(xi) < 1e-10) {
      return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
    }
    else {
      if (all(1 + (xi * z) / sigma > 0)) {
        return(-(length(z) * log(sigma)) - ((1 + 1 / xi)*(sum(log(1 + (xi * z) / sigma)))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}


transform_to_exp <- function (y,sig, xi){
  std_exp <- (1 / xi) * log( 1 + xi * (y/sig))  
  return(std_exp)
}

