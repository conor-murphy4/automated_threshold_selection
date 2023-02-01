source('_gpd.R') #rgpd and qgpd function

#GPD likelihood
GPD_LL<-function(par,z){
  sigma<-par[1]
  xi<-par[2]
  if(sigma>0){
    if(abs(xi)<1e-10){
      return(-length(z)*log(sigma)-((1/sigma)*sum(z)))
    }
    else{
      if( all(1+(xi*z)/sigma >0)){
        return(-(length(z)*log(sigma))-((1+1/xi)*(sum(log(1+(xi*z)/sigma)))))
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

#threshold selection method

library(roxygen2)
vignette("rd")

#' Threshold selection for univariate extremes
#' 
#' 'thr_dist' selects a constant threshold above which the data can be most closely modelled by a Generalised Pareto distribution.
#' 
#' Detailed description
#' 
#' @param data A numeric vector.
#' @param thresh A numeric vector of proposed thresholds to test.
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#' 
#' @returns A list.
#' 
#' @examples 
#' test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' Q1=seq(0, 2, by=0.1)
#' example <- thr_dist(test1, thresh=Q1)
#' 
#' test2 <- rgpd(10000, shape = 0.1, scale=0.5)
#' u=1
#' cens_thr<-u*rbeta(length(test2),1,0.5)
#' keep <- test>cens_thr
#' data_test2 <- test[keep]
#' Q2=seq(0, 3, by=0.1)
#' example2 <- thr_dist(data_test2,thresh=Q2)


thr_dist <- function(data,thresh,k=100,m=500){
  if(length(dim(data)) != 0){
    stop("Data must be a vector")
  }
  if(length(dim(thresh)) != 0){
    stop("u to be tested needs to be a vector")
  }
  if(k <= 0 | k%%1 != 0){
    stop("Number of bootstrapped samples must be a positive integer")
  }
  if(m <= 0 | m%%1 != 0){
    stop("Number of equally spaced probabilities must be a positive integer")
  }
  meandistances <- xis <- sigmas <- lens <- numeric(length(thresh))
  for(i in 1:length(thresh)){
    u <- thresh[i]
    data_ex <- data - u #quantile(data, q, names=FALSE)
    excess <- data_ex[data_ex >0]
    lens[i] <- length(excess)
    if(lens[i]>10){
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL,z=excess,par=c(mle0,0.1),control=list(fnscale=-1))
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      for(j in 1:k){
        X <- sample(excess, lens[i], replace=TRUE)
        mle <- mean(X)
        ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit<-optim(GPD_LL,z=X,par=pars_init,control=list(fnscale=-1))
        quants<-qgpd(c(1:m)/(m+1),shape=gpd.fit$par[[2]],scale = gpd.fit$par[[1]])
        distances[j] <- (1/m)*sum(abs(quantile(X,probs=c(1:m)/(m+1))-quants))
      }
      meandistances[i] <- mean(distances)
    }
    else{
      meandistances[i] <- NA
    }
  }
  choice <- thresh[which.min(meandistances)]
  xi <- xis[which.min(meandistances)]
  sigma <- sigmas[which.min(meandistances)]
  len <- lens[which.min(meandistances)]
  result <- list(thresh=choice, par=c(sigma,xi), num=len, dists=meandistances)
  return(result)
}
