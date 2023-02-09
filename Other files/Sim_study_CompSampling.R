setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
source("STORM_input/thr_selection_final.R")
#source("STORM_input/Wadsworth2016.R")
library(threshr)

#Idealistic with small sample size
n1 <- 500
Time_mythr <- Time_poi <- numeric(n1)
mythresh <- poithresh <- numeric(n1)
for(ii in 1:n1){
  dat1 <- runif(200, 0.5, 1.0)
  dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
  data <- c(dat1, dat2)
  thresh <- quantile(data, seq(0,0.95,by=0.05))
  #My method
  start_mythr <- Sys.time()
  mythr <- thr_dist(data, thresh=thresh)$thresh
  end_mythr <- Sys.time()
  mythresh[ii] <- mythr 
  Time_mythr[ii] <- end_mythr-start_mythr
  #My method with Poi
  start_mythr_poi <- Sys.time()
  mythrpoi <- thr_dist_poi(data, thresh=thresh)$thresh
  end_mythr_poi <- Sys.time()
  poithresh[ii] <- mythrpoi 
  Time_poi[ii] <- end_mythr_poi-start_mythr_poi
  
}

mythrI<- data.frame(thr=mythresh, time=Time_mythr)
poithrI <- data.frame(thr=poithresh, time=Time_poi)

saveRDS(mythrI, "STORM_output/Sampcomp_mythr.rds")
saveRDS(poithrI, "STORM_output/Sampcomp_poi.rds")

