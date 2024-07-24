source("src/eqd_quantile_definition.R")
library(Metrics)

#Case 1
n <- 500
probs=seq(0, 0.95, by=0.05)
types <- c(6,7,6)
probs_defs <- c(1,2,2)
for(jj in 1:3){
  mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
  type <- types[jj]
  prob_def <- probs_defs[jj]
  for(ii in 1:n){
    print(ii)
    set.seed(ii)
    dat1 <- runif(200, 0.5, 1.0)
    dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
    data <- c(dat1, dat2)
    thresh <- quantile(data, probs, names=F)
    
    #EQD method
    myres <- eqd_quant_def(data, thresh, type=type, probs_def = prob_def)
    mythresh[ii] <- myres$thresh
    myquantile[ii] <- probs[thresh==mythresh[ii]]
    myscale[ii] <- myres$par[1]
    myshape[ii] <- myres$par[2]
  }

  eqdthr <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
  
  filename <- paste("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef",jj,".csv",sep="")
  write.csv(eqdthr, filename)

}


eqd_case1_probdef1 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef1.csv", row.names = 1, header = T)
eqd_case1_probdef2 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef2.csv", row.names = 1, header = T)
eqd_case1_probdef3 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef3.csv", row.names = 1, header = T)


#Case 4
n <- 500
probs=seq(0, 0.95, by=0.05)
types <- c(6,7,6)
probs_defs <- c(1,2,2)
for(jj in 1:3){
  mythresh <- myquantile <- myscale <- myshape <- mylen <- numeric(n)
  type <- types[jj]
  prob_def <- probs_defs[jj]
  for(ii in 1:n){
    print(ii)
    set.seed(ii)
    data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
    cens_thr<-rbeta(length(data_all),1,2)
    data_above <- sample(data_all[data_all > 1], 279, replace=FALSE)
    data_below <- sample(data_all[data_all > cens_thr & data_all <= 1], 721, replace = FALSE)
    data <- c(data_below, data_above)
    thresh <- quantile(data, probs, names=F)
    
    #EQD method
    myres <- eqd_quant_def(data, thresh, type=type, probs_def = prob_def)
    mythresh[ii] <- myres$thresh
    myquantile[ii] <- probs[thresh==mythresh[ii]]
    myscale[ii] <- myres$par[1]
    myshape[ii] <- myres$par[2]
  }
  
  eqdthr <- data.frame(thr=mythresh,quantile=myquantile, scale=myscale, shape=myshape, len=mylen)
  
  filename <- paste("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case4_probdef",jj,".csv",sep="")
  write.csv(eqdthr, filename)
  
}

eqd_case1_probdef1 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef1.csv", row.names = 1, header = T)
eqd_case1_probdef2 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef2.csv", row.names = 1, header = T)
eqd_case1_probdef3 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case1_probdef3.csv", row.names = 1, header = T)


rmse(eqd_case1_probdef1$thr,1)
rmse(eqd_case1_probdef2$thr,1)
rmse(eqd_case1_probdef3$thr,1)

bias(eqd_case1_probdef1$thr,1)
bias(eqd_case1_probdef2$thr,1)
bias(eqd_case1_probdef3$thr,1)

var(eqd_case1_probdef1$thr)
var(eqd_case1_probdef2$thr)
var(eqd_case1_probdef3$thr)

eqd_case4_probdef1 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case4_probdef1.csv", row.names = 1, header = T)
eqd_case4_probdef2 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case4_probdef2.csv", row.names = 1, header = T)
eqd_case4_probdef3 <- read.csv("output/threshold_selection/additional_sensitivity/prob_definition/eqd_case4_probdef3.csv", row.names = 1, header = T)

rmse(eqd_case4_probdef1$thr,1)
rmse(eqd_case4_probdef2$thr,1)
rmse(eqd_case4_probdef3$thr,1)

bias(eqd_case4_probdef1$thr,1)
bias(eqd_case4_probdef2$thr,1)
bias(eqd_case4_probdef3$thr,1)

var(eqd_case4_probdef1$thr)
var(eqd_case4_probdef2$thr)
var(eqd_case4_probdef3$thr)
