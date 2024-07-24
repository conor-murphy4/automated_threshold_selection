
library(evir)
source("src/eqd.R")
source("src/JointMLEFunctions.r")
library(threshr)

data("nidd.thresh")

# Table 7:  -----------------------------------------------------

#Candidate threshold grids
thresholds_list <- list(quantile(nidd.thresh, seq(0,0.93, by=0.01)),quantile(nidd.thresh, seq(0,0.9, by=0.01)), quantile(nidd.thresh, seq(0,0.8, by=0.01)), 
                        quantile(nidd.thresh, seq(0,0.8, by=0.2)), quantile(nidd.thresh, seq(0,0.9, by=0.3)), quantile(nidd.thresh, seq(0,0.75, by=0.25)),
                        quantile(nidd.thresh, c(0,0.1,0.4,0.7)))

EQD_thr <- wads_thr <- north_thr <- numeric(7)
for(i in 1:7){
  thresholds <- thresholds_list[[i]]
  #EQD
  set.seed(11111) 
  EQD_thr[i] <- eqd(nidd.thresh, thresholds, k=200)$thresh
  #Wadsworth
  set.seed(11111)
  wads_thr[i] <- NHPP.diag(nidd.thresh, u=thresholds, plot.out=FALSE, UseQuantiles = FALSE)$thresh
  #Northrop
  set.seed(11111)
  north_thr[i] <- summary(ithresh(nidd.thresh, u_vec = thresholds))[3]
}

(table_7 <- data.frame(EQD=EQD_thr, Wadsworth=wads_thr, Northrop=north_thr))

write.csv(table_7, "output/tables/Table_7_River_Nidd_dataset_selected_thresholds.csv")
