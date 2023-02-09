#Beirlant 2017 Truncated GPD
trunc_GPDLL <- function(par, z){
  sig <- par[1]
  xi <- par[2]
  if(sig>0){
    if(all(1+(xi*z)/sig >0)){
      return(-(length(z)*log(sig))-((1+1/xi)*(sum(log(1+(xi*z)/sig)))) -(length(z)*log(1-(1+(xi*max(z))/sig)^(-1/xi))))
    }
    else{
      return(-1e6)
    }
  }
  else{
    return(-1e7)
  }
  
}

trunc_gpd_endpoint <- function(data, par, u){
  k <- length(data[data>u]) + 1
  scale <- par[1]
  shape <- par[2]
  tau <- shape/scale
  p <- (1+tau*(max(data)-u))^(-1/shape)
  if(p > 1/k){
    end <- u + (1/tau)*(((1-1/k)/(p-1/k))^shape-1)
    return(list(end=end, shape=shape, Dtk = 1))
  }
  else{
    if(shape < 0){
      end <- u - 1/tau
    } 
    else{
      end <- 1e10
    }
    return(list(end=end, shape=shape, Dtk = 0))
  }
}

dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape = 0.2, scale=0.5, mu=1.0)
data <- c(dat1,dat2)


thresh <- quantile(data, seq(0,0.95,by=0.01), names=F)
(thr <- thr_dist(data,thresh = thresh))

excess <- data[data > thr$thresh] - thr$thresh
(opt_trunc <- optim(trunc_GPDLL, par=c(0.1,0.1), z=excess, control = list(fnscale=-1)))
(endpoint_trunc <- trunc_gpd_endpoint(data, par=opt_trunc$par, u = thr$thresh))


#Simulation study
#Taking selected thresholds from previous simulations and comparing endpoint estimators

mythrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/C_fixed n/mythrC.rds")
mythr3A <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythr3A.rds")
mythrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/mythrG.rds")

#GPD with negative shape
#truncation point
xik <- trunc_point <- Dtk <- numeric(500)
set.seed(11111)
for(i in 1:500){
  dat1 <- runif(200, 0.5,1)
  dat2 <- rgpd(1000, shape=-0.2, scale=0.5, mu=1.0)
  data <- c(dat1,dat2)
  thr <- mythr3A$thr[i]
  excess <- data[data > thr] - thr
  opt_trunc <- optim(trunc_GPDLL, par=c(0.1,0.1), z=excess, control = list(fnscale=-1))
  beirlant <- trunc_gpd_endpoint(data, par=opt_trunc$par, u=thr)
  trunc_point[i] <- beirlant$end 
  xik[i] <- beirlant$shape
  Dtk[i] <- beirlant$Dtk
}
#endpoint
neg <- mythr3A$shape < 0
endpoint <- mythr3A$thr[neg] - mythr3A$scale[neg]/mythr3A$shape[neg] 
library(Metrics)
rmse(endpoint, 3.5)
rmse(trunc_point, 3.5)
#plots
par(mfrow=c(1,2))
plot(endpoint, trunc_point, xlab = "GPD", ylab = "TGPD", main="")
abline(a=0, b=1)
plot(mythr3A$shape, xik, xlab = "Shape", ylab = "Truncated Shape", main="")
length(Dtk[Dtk == 1])

length(mythrC$shape[mythrC$shape < 0])
length(xik_C[xik_C < 0])

#Censored case
#truncation point
xik_C <- trunc_point_C <- Dtk_C <- numeric(500)
set.seed(11111)
for(i in 1:500){
  data_all <- rgpd(4000, shape=0.1, scale=0.5)
  cens_thr <- rbeta(4000, 1, 0.5)
  keep <- data_all>cens_thr
  data_kept <- data_all[keep]
  data <- sample(data_kept, 1000, replace = FALSE)
  thr <- mythrC$thr[i]
  excess <- data[data > thr] - thr
  opt_trunc <- optim(trunc_GPDLL, par=c(0.1,0.1), z=excess, control = list(fnscale=-1))
  beirlant <- trunc_gpd_endpoint(data, par=opt_trunc$par, u=thr)
  trunc_point_C[i] <- beirlant$end 
  xik_C[i] <- beirlant$shape
  Dtk_C[i] <- beirlant$Dtk
}

#endpoint
neg <- mythrC$shape < 0
endpoint_C <- mythrC$thr[neg] - mythrC$scale[neg]/mythrC$shape[neg] 
endpoint_C[36:500] <- 4000
#trunc_point_C[trunc_point_C>1000] <- 100
#plots
length(trunc_point_C[trunc_point_C == 1e10])
plot(endpoint_C, trunc_point_C, xlab = "GPD", ylab = "TGPD", main="")
abline(a=0, b=1)
plot(mythrC$shape, xik_C, xlab = "Shape", ylab = "Truncated Shape", main="")

#Gaussian
#truncation point
xik_G <- trunc_point_G <- Dtk_G <- numeric(500)
set.seed(12345)
for(i in 1:500){
  data <- rnorm(2000)
  thr <- mythrG$thr[i]
  excess <- data[data > thr] - thr
  opt_trunc <- optim(trunc_GPDLL, par=c(0.1,0.1), z=excess, control = list(fnscale=-1))
  beirlant <- trunc_gpd_endpoint(data, par=opt_trunc$par, u=thr)
  trunc_point_G[i] <- beirlant$end 
  xik_G[i] <- beirlant$shape
  Dtk_G[i] <- beirlant$Dtk
}

#endpoint
neg <- mythrG$shape < 0
endpoint_G <- mythrG$thr[neg] - mythrG$scale[neg]/mythrG$shape[neg] 
endpoint_G[498:500] <- 40
#plots
plot(endpoint_G, trunc_point_G, xlab = "GPD", ylab = "TGPD", main="")
abline(a=0, b=1)
plot(mythrG$shape, xik_G, xlab = "Shape", ylab = "Truncated Shape", main="")
