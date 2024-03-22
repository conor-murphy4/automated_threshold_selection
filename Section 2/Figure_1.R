source("src/parameter_stability.R")

Case_4 <- read.csv("data/Case_4.csv", row.names = 1)
case_4_sample_1 <- Case_4[,1]

data("nidd.thresh")

Q <-  seq(0,0.95,by=0.05)

#Candidate threshold grids
thresh_4 <- quantile(case_4_sample_1, Q, names=FALSE)
thresh_Nidd <- quantile(nidd.thresh,Q, names=F)

#Plotting window
dev.new(width=9.17, height=3.7,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

#Parameter stability plots
shapestabboot(case_4_sample_1, thresholds = thresh_4, Q=Q, reverse=F, boot=TRUE,m.boot=200)
abline(v=1.0, col="green")

shapestabboot(nidd.thresh, thresholds = thresh_Nidd, Q=Q, reverse=F, boot=TRUE,m.boot=200)