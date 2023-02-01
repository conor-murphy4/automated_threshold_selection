source('_gpd.R') #rgpd and qgpd function
required_pckgs <- c( "dplyr", "purrr")
#install.packages(required_pckgs, dependencies = TRUE)
lapply(required_pckgs, require, character.only = TRUE)

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

thrs_distance <- function(data, thresh, B = 100, m = 501, verbose = FALSE){
  #checks on the inputs
  if(length(dim(data)) != 0){
    stop("Data must be a vector")
  }
  if(length(dim(thresh)) != 0){
    stop("u to be tested needs to be a vector")
  }
  if(B <= 0 | B%%1 != 0){
    stop("Number of bootstrapped samples must be a positive integer")
  }
  if(m <= 0 | m%%1 != 0){
    stop("Number of equally spaced probabilities must be a positive integer")
  }
  
  #wrangle data
  data <- data.frame(y = data)
  n <- NROW(data)
  thresholds <- map(.x = thresh, .f = function(x){rep(x,n)})
  
  #Now will create a list for each threshold with a data frame of exceedances in each list
  data_filtered_list <- map(
    .x = thresholds,
    .f = filter_data,
    data = data)
  
  #converts data frame to a vector of excesses
  data_list_vecs <- map(data_filtered_list, function(df){df$y - df$v})
  thresholds_list_vecs <- map(data_filtered_list, function(df){df$v})
  
  #obtain the length of each vector
  nu <- map_dbl(.x = map(.x = data_filtered_list, NROW), .f = function(vec){vec[1]})
  
  #remove those thresholds where fewer than 10 exceedances are present
  if(any(nu < 10)){
    index <- which(nu < 10)
    data_list_vecs <- data_list_vecs[-index]
    thresholds_list_vecs <- thresholds_list_vecs[-index]
    nu <- nu[-index]
    thresh <- thresh[-index]
  }
  
  #fit the models to get the point estimates
  init_mles <- pmap(.l = list(data = data_list_vecs),
                    .f = get_mle,
                    llh_val = FALSE)
  
  #extract mles as vectors
  sig_u_mles <- map_dbl(.x = init_mles, .f = function(vec){vec[[1]]})
  xi_mles <- map_dbl(.x = init_mles, .f = function(vec){vec[[2]]})
  
  #setting up the equally spaced probabilities
  p <- ppoints(m)
  
  #empty vector to save the output
  dists <- c()
  
  for(i in seq_along(thresholds_list_vecs)){
    #create a list of non-parametrically bootstrapped samples
    boot_samples <- replicate(n = B, 
                              expr = sample(x = data_list_vecs[[i]], size = nu[i], replace = TRUE),
                              simplify = FALSE)
    
    #initialise starting parameters for bootstrapping
    if(xi_mles[i] <= 0){
      par <- rbind(sapply(X = boot_samples, FUN = mean), rep(0.1, times = B))
      #conver the 
      par <- lapply(X = seq_len(ncol(par)), FUN = function(j){par[,j]})
    }
    else{
      par <- rep(list(c(sig_u_mles[i], xi_mles[i])), B)
    }
    
    #get the mles of the bootstrapped samples
    boot_mles <- pmap(.f = get_mle,
                      .l = list(data = boot_samples, par=par),
                      llh_val = FALSE)
    
    #extract mles as vectors
    sig_boot_mles <- map_dbl(.x = boot_mles, .f = function(vec){vec[[1]]})
    xi_boot_mles <- map_dbl(.x = boot_mles, .f = function(vec){vec[[2]]})
    
    #calculate the theoretical quantiles using the bootstrapped MLES
    tq <- pmap(.f = qgpd,
               .l = list(scale = sig_boot_mles, shape = xi_boot_mles),
               p = p)
    
    #calculate the empirical quantiles of the bootstrapped samples
    eq <- pmap(.f = quantile,
               .l = list(x = boot_samples),
               p = p)
    
    #obtain the distance metrics
    dq1 <- pmap(.f = dq1_metric, .l = list(x = eq, y = tq))
    dists[i] <- mean(map_dbl(.x = dq1, .f = function(vec){vec[1]}))
    
    if(verbose){
      print(paste(i, "thresholds of", length(thresholds_list_vecs), "completed"))
    }
  }
  
  #obtain the choice
  ind <- which.min(dists)
  choice <- thresh[ind]
  xi <- xi_mles[ind]
  sig <- sig_u_mles[ind]
  len <- nu[ind]
  out <- list(u = choice,
              par = c(sig, xi),
              nu = len,
              dists = dists)
  return(out)
}

#test
dat1 <- runif(200, 0.5,1)
dat2 <- rgpd(1000, shape=0.2, scale=0.5, mu=1)
data <- c(dat1,dat2)
thresholds <- quantile(data,seq(0,0.95,by=0.01), names = F)

start_n <- Sys.time()
(new_thr <- thrs_distance(data, thresh = thresholds))
end_n <- Sys.time()

(Time_n <- end_n - start_n)

