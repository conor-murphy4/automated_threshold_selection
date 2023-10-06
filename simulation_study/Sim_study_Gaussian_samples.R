#setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
#source("STORM_input/thr_selection_final.R")
source("thr_selection_final.R")
# source("STORM_input/Wadsworth2016.R")
source("Wadsworth2016.R")
library(threshr)

#Gaussian data
n4 <- 500
# Time_mythr <- Time_wads <- Time_north <- Time_north_x <- numeric(n4)
mythresh <- wadsthresh <- norththresh <- numeric(n4)
myquantile <- wadsquantile <- northquantile <- numeric(n4)
myscale <- wadsscale <- northscale <- numeric(n4)
myshape <- wadsshape <- northshape <- numeric(n4)
mylen <- wadslen <- northlen <- numeric(n4)
for(ii in 1:n4){
  set.seed(ii)
  data <- rnorm(20000)
  probs=seq(0.5, 0.95, by=0.05)
  thresh <- quantile(data, probs=probs, names=FALSE)
  #My method
  # start_mythr <- Sys.time()
  myres <- thr_dist(data, thresh=thresh)
  # end_mythr <- Sys.time()
  mythresh[ii] <- myres$thresh
  myquantile[ii] <- probs[thresh==mythresh[ii]]
  myscale[ii] <- myres$par[1]
  myshape[ii] <- myres$par[2]
  mylen[ii] <- myres$num
  # Time_mythr[ii] <- end_mythr-start_mythr
  #Wadsworth 2016 method
  # start_wads <- Sys.time()
  wadsthresh[ii] <- NHPP.diag(data, u=thresh, plot.out = FALSE, UseQuantiles = FALSE)$thresh[[1]]
  #end_wads <- Sys.time()
  if(!is.na(wadsthresh[ii])){
    fit.data_w <- data[data>wadsthresh[ii]] - wadsthresh[ii]
    optwads <- optim(GPD_LL, z=fit.data_w, par=c(mean(fit.data_w), 0.1), control=list(fnscale=-1))
    wadsscale[ii] <- optwads$par[1]
    wadsshape[ii] <- optwads$par[2]
    wadslen[ii] <- length(fit.data_w)
    wadsquantile[ii] <- probs[thresh==wadsthresh[ii]]
  }
  else{
    wadsscale[ii] <- NA
    wadsshape[ii] <- NA
    wadslen[ii] <- NA
    wadsquantile[ii] <- NA
  }
  # Time_wads[ii] <- end_wads-start_wads
  #Northrop 2017 method
  # start_north <- Sys.time()
  norththresh[ii] <- summary(ithresh(data, u_vec = thresh))[3]
  northquantile[ii] <- probs[thresh==norththresh[ii]]
  fit.data_n <- data[data>norththresh[ii]] - norththresh[ii]
  optnorth <- optim(GPD_LL, z=fit.data_n, par=c(mean(fit.data_n), 0.1), control=list(fnscale=-1))
  northscale[ii] <- optnorth$par[1]
  northshape[ii] <- optnorth$par[2]
  northlen[ii] <- length(fit.data_n)
  # Time_north[ii] <- end_north-start_north
  # Time_north_x[ii] <- end_north1 - start_north
}

mythrG <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)#, time=Time_mythr)
wadsthrG <- data.frame(thr=wadsthresh,quantile=wadsquantile, scale=wadsscale, shape=wadsshape, len=wadslen)#, time=Time_wads)
norththrG <- data.frame(thr=norththresh, quantile=northquantile, scale=northscale, shape=northshape, len=northlen)#, time=Time_north, timeX = Time_north_x)

# saveRDS(mythrG, "STORM_output/mythrG.rds")
# saveRDS(wadsthrG, "STORM_output/wadsthrG.rds" )
# saveRDS(norththrG, "STORM_output/norththrG.rds")

saveRDS(mythrG, "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/mythrG_large_sample.rds")
saveRDS(wadsthrG, "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/wadsthrG_large_sample.rds" )
saveRDS(norththrG, "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/norththrG_large_sample.rds")

#Orginal samples summaries of chosen thresholds
mythrG_sort <- mythrG[order(mythrG$thr),]
wadsthrG_sort <- wadsthrG[order(wadsthrG$thr),]
norththrG_sort <- norththrG[order(norththrG$thr),]

(medians <- c(mythrG_sort$thr[0.5*n4], wadsthrG_sort$thr[0.5*n4], norththrG_sort$thr[0.5*n4]))
(means <- c(mean(mythrG_sort$thr), mean(wadsthrG_sort$thr, na.rm = TRUE), mean(norththrG_sort$thr)))
(CI_lower <- c(mythrG_sort$thr[floor(0.025*n4)], wadsthrG_sort$thr[floor(0.025*n4)], norththrG_sort$thr[floor(0.025*n4)]))
(CI_upper <- c(mythrG_sort$thr[ceiling(0.975*n4)], wadsthrG_sort$thr[ceiling(0.975*n4)], norththrG_sort$thr[ceiling(0.975*n4)]))

(medians <- c(mythrG_sort$quantile[0.5*n4], wadsthrG_sort$quantile[0.5*n4], norththrG_sort$quantile[0.5*n4]))
(means <- c(mean(mythrG_sort$quantile), mean(wadsthrG_sort$quantile, na.rm = TRUE), mean(norththrG_sort$quantile)))
(CI_lower <- c(mythrG_sort$quantile[floor(0.025*n4)], wadsthrG_sort$quantile[floor(0.025*n4)], norththrG_sort$quantile[floor(0.025*n4)]))
(CI_upper <- c(mythrG_sort$quantile[ceiling(0.975*n4)], wadsthrG_sort$quantile[ceiling(0.975*n4)], norththrG_sort$quantile[ceiling(0.975*n4)]))

#Larger samples summaries of chosen thresholds
mythrG_large_sort <- mythrG_large_sample[order(mythrG_large_sample$thr),]
wadsthrG_large_sort <- wadsthrG_large_sample[order(wadsthrG_large_sample$thr),]
norththrG_large_sort <- norththrG_large_sample[order(norththrG_large_sample$thr),]

(medians <- c(mythrG_large_sort$thr[0.5*n4], wadsthrG_large_sort$thr[0.5*n4], norththrG_large_sort$thr[0.5*n4]))
(means <- c(mean(mythrG_large_sort$thr), mean(wadsthrG_large_sort$thr, na.rm = TRUE), mean(norththrG_large_sort$thr)))
(CI_lower <- c(mythrG_large_sort$thr[floor(0.025*n4)], wadsthrG_large_sort$thr[floor(0.025*n4)], norththrG_large_sort$thr[floor(0.025*n4)]))
(CI_upper <- c(mythrG_large_sort$thr[ceiling(0.975*n4)], wadsthrG_large_sort$thr[ceiling(0.975*n4)], norththrG_large_sort$thr[ceiling(0.975*n4)]))

(medians <- c(mythrG_large_sort$quantile[0.5*n4], wadsthrG_large_sort$quantile[0.5*n4], norththrG_large_sort$quantile[0.5*n4]))
(means <- c(mean(mythrG_large_sort$quantile), mean(wadsthrG_large_sort$quantile, na.rm = TRUE), mean(norththrG_large_sort$quantile)))
(CI_lower <- c(mythrG_large_sort$quantile[floor(0.025*n4)], wadsthrG_large_sort$quantile[floor(0.025*n4)], norththrG_large_sort$quantile[floor(0.025*n4)]))
(CI_upper <- c(mythrG_large_sort$quantile[ceiling(0.975*n4)], wadsthrG_large_sort$quantile[ceiling(0.975*n4)], norththrG_large_sort$quantile[ceiling(0.975*n4)]))
