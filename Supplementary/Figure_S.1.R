source("src/parameter_stability.R")

Case_1 <- read.csv("data/Case_1.csv", row.names = 1)
Case_2 <- read.csv("data/Case_2.csv", row.names = 1)
Case_3 <- read.csv("data/Case_3.csv", row.names = 1)

case_1_sample_1 <- Case_1[,1]
case_2_sample_1 <- Case_2[,1]
case_3_sample_1 <- Case_3[,1]

Q <-  seq(0,0.95,by=0.05)

#Candidate threshold grids
thresholds_1 <- quantile(case_1_sample_1,Q, names=F)
thresholds_2 <- quantile(case_2_sample_1,Q, names=F)
thresholds_3 <- quantile(case_3_sample_1,Q, names=F)

#Plotting window
dev.new(width=9.17, height=3.7,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')

#Parameter stability plot
shapestabboot(case_1_sample_1, thresholds = thresholds_1, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")
shapestabboot(case_2_sample_1, thresholds = thresholds_2, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")
shapestabboot(case_3_sample_1, thresholds = thresholds_3, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")