source("src/parameter_stability.R")

Case_4 <- read.csv("data/Case_4.csv", row.names = 1)
case_4_sample_1 <- Case_4[,1]

data("nidd.thresh")

#Candidate threshold grids
Q_4 <- seq(0,0.95,by=0.05)
thresh_4 <- quantile(case_4_sample_1, Q_4, names=FALSE)

Q_Nidd <- seq(0,0.94,by=0.04)
thresh_Nidd <- quantile(nidd.thresh,Q_Nidd, names=F)

#Plotting window
dev.new(width=9.17, height=3.7,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

#Parameter stability plots
shapestabboot(case_4_sample_1, thresholds = thresh_4, Q=Q_4, reverse=F, boot=TRUE,m.boot=200)
abline(v=1.0, col="green")

shapestabboot(nidd.thresh, thresholds = thresh_Nidd, Q=Q_Nidd, reverse=F, boot=TRUE,m.boot=200)