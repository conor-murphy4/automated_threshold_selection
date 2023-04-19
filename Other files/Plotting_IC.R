library("Metrics")

#Data

mythrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI.rds")
mythrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI1.rds")
mythrI1.1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/mythrI1.1.rds")
mythrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI2.rds")
mythrI2.2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/mythrI2.2.rds")
mythrI2.3 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/mythrI2.3.rds")
mythrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/C_fixed n/mythrC.rds")
mythrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/mythrG.rds")


norththrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/norththrI.rds")
norththrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/norththrI1.rds")
norththrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/norththrI2.rds")
norththrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/C_fixed n/norththrC.rds")
norththrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/norththrG.rds")

wadsthrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/wadsthrI.rds")
wadsthrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/wadsthrI1.rds")
wadsthrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/wadsthrI2.rds")
wadsthrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/C_fixed n/wadsthrC.rds")
wadsthrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/wadsthrG.rds")

plot(wadsthrI2$time, ylim=c(0,15))
points(mythrI2$time, col="red")

#Ideal case

# plot(mythrI$thr, ylab = "Threshold Choice", xlab = "Iteration", ylim=c(0.7,3))
# points(wadsthrI$thr, col="red", pch=6)
# points(norththrI$thr, col="blue", pch=0)
# legend("topleft", c("Mythr", "Wadsworth", "Northrop"), col=c("black", "red", "blue"), pch=c(1,6,0))


#Density plots
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')

err <- is.na(wadsthrI$thr)
d_my <- density(mythrI$thr)
d_wads <- density(wadsthrI$thr[!err])
d_north <- density(norththrI$thr)
plot(d_wads, ylab="Density", xlab="Thresholds", main="Case 1", lwd=2, xlim=c(0.5,3.0), col="grey")
lines(d_my, lwd=2, lty=2, col="green")
lines(d_north, lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Wadsworth", "Northrop"), col=c("green", "grey", "blue"), lty=c(2,1,3),lwd=2, cex=0.8)

err <- is.na(wadsthrI1$thr)
d_my <- density(mythrI1$thr)
d_wads <- density(wadsthrI1$thr[!err])
d_north <- density(norththrI1$thr)
plot(d_wads, ylab="Density", xlab="Thresholds", main="Case 2",ylim=c(0,50), lwd=2, xlim=c(0.5,3.0), col="grey")
lines(d_my, lwd=2, lty=2, col="green")
lines(d_north, lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Wadsworth", "Northrop"), col=c("green", "grey", "blue"), lty=c(2,1,3),lwd=2, cex=0.8)

err <- is.na(wadsthrI2$thr)
d_my <- density(mythrI2$thr)
d_wads <- density(wadsthrI2$thr[!err])
d_north <- density(norththrI2$thr)
plot(d_wads, ylab="Density", xlab="Thresholds", main="Case 3", lwd=2, xlim=c(0.5,3.0), col="grey")
lines(d_my, lwd=2, lty=2, col="green")
lines(d_north, lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Wadsworth", "Northrop"), col=c("green", "grey", "blue"), lty=c(2,1,3),lwd=2, cex=0.8)
length(mythrI2$thr<1.05 & mythrI2$thr>0.95)
length(wadsthrI2$thr<1.05 & wadsthrI2$thr>0.95)
summary(mythrI2$thr)
summary(wadsthrI2$thr)

err <- is.na(wadsthrC$thr)
length(wadsthrC$thr[err])
d_my <- density(mythrC$thr)
d_wads <- density(wadsthrC$thr[!err])
d_north <- density(norththrC$thr)
plot(d_wads,ylab="Density", xlab="Thresholds", main="Case 4", lwd=2, xlim=c(0,3.0),col="grey" )
lines(d_my, lwd=2, col="green", lty=2)
lines(d_north, lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Wadsworth", "Northrop"), col=c("green", "grey", "blue"), lty=c(2,1,3),lwd=2, cex=0.8)

#G
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
err <- is.na(wadsthrG$thr)
length(wadsthrG$thr[err])
d_my <- density(mythrG$thr)
d_wads <- density(wadsthrG$thr[!err])
d_north <- density(norththrG$thr)
plot(d_wads,ylab="Density", xlab="Thresholds", main="Gaussian Case", lwd=2,col="grey" )
lines(d_my, lwd=2, col="green", lty=2)
lines(d_north, lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Wadsworth", "Northrop"), col=c("green", "grey", "blue"), lty=c(2,1,3),lwd=2, cex=0.8)
wadsthrG$
#I1.1
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
hist(mythrI1.1$thr, xlab="Threshold", breaks = 30, main="")

plot(density(mythrI1.1$thr))

dev.new(width=8.00, height=10.00,noRStudioGD = TRUE)
par(mfrow=c(3,1),bg='transparent')
hist(danthrI$thr, xlab="Threshold", breaks = 30, main="Case 1")
hist(danthrI1$thr, xlab="Threshold", breaks = 30, main="Case 2")
hist(danthrC$thr, xlab="Threshold", breaks = 30, main = "Case 4")


#Censored case

# plot(mythrC$thr, ylab = "Threshold Choice", xlab = "Iteration", ylim=c(0,3))
# points(wadsthrC$thr, col="red", pch=6)
# points(norththrC$thr, col="blue", pch=0)
# legend("topleft", c("Mythr", "Wadsworth", "Northrop"), col=c("black", "red", "blue"), pch=c(1,6,0))
# 







#Overlaying histograms
# par(mfrow=c(2,1))
# c1 <- rgb(250,50,130,max = 255, alpha = 80, names = "blue")
# c2 <- rgb(130,250,150, max = 255, alpha = 80, names = "red")
# c3 <- rgb(50,130,250, max = 255, alpha = 255, names = "black")
# A <- mythrI$thr
# B <- wadsthrI$thr[!is.na(wadsthrI$thr)]
# C <- norththrI$thr
# A1 <- mythrC$thr
# B1 <- wadsthrC$thr[!is.na(wadsthrC$thr)]
# C1 <- norththrC$thr
# ax=seq(0.8,3.0,by=0.05)
# hgA <- hist(A, breaks = ax , plot = FALSE) # Save first histogram data
# hgB <- hist(B, breaks = ax, plot = FALSE) # Save 2nd histogram data
# hgC <- hist(C, breaks = ax, plot = FALSE)
# hgA1 <- hist(A1, breaks = ax, plot = FALSE) # Save first histogram data
# hgB1 <- hist(B1, breaks = ax, plot = FALSE) # Save 2nd histogram data
# hgC1 <- hist(C1, breaks = ax, plot = FALSE)
# 
# plot(hgA, col = c3, main="", xlab="Thresholds") # Plot 1st histogram using a transparent color
# plot(hgC, col = c1, add = TRUE) # Add 2nd histogram using different color
# plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
# legend("topright", c("Ours", "Wadsworth", "Northrop"), col=c(c3, c2, c1), pch=c(15,15,15), cex=0.8)
# 
# plot(hgA1, col = c3, main="", xlab="Thresholds") # Plot 1st histogram using a transparent color
# plot(hgC1, col = c1, add = TRUE) # Add 2nd histogram using different color
# plot(hgB1, col = c2, add = TRUE) # Add 2nd histogram using different color
# legend("topright", c("Ours", "Wadsworth", "Northrop"), col=c(c3, c2, c1), pch=c(15,15,15), cex=0.8)
# 

#--------------------------------------

mythrI1.1 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/SimIC/mythrI1.1.rds")
mythrI2.1 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/SimIC/mythrI2.1.rds")

hist(mythrI1.1$thr, breaks = 30)
hist(mythrI2.1$thr, breaks = 30)

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(density(mythr2A$thr),  lwd=2, col="green", lty=2, xlim=c(0.5,2.5), main="Case 2A", ylab="Density", xlab="Thresholds")
lines(density(norththr2A$thr), lwd=2, col="blue", lty=3)
abline(v=1.0, col="red", lwd=1.5)
legend("topright", c("Our method", "Northrop"), col=c("green", "blue"), lty=c(2,3),lwd=2, cex=0.8)
plot(density(mythr3A$thr),  lwd=2, col="green", lty=2, xlim=c(0.5,2.5), main="Case 3A", ylab="Density", xlab="Thresholds")
lines(density(norththr3A$thr), lwd=2, col="blue", lty=3)
abline(v=1.0, col="red", lwd=1.5)
legend("topright", c("Our method", "Northrop"), col=c("green", "blue"), lty=c(2,3),lwd=2, cex=0.8)
plot(density(mythr3B$thr),  lwd=2, col="green", lty=2, xlim=c(0.5,2.5), main="Case 3B", ylab="Density", xlab="Thresholds")
lines(density(norththr3B$thr), lwd=2, col="blue", lty=3)
abline(v=1.0, col="red", lwd=1.5)
legend("topright", c("Our method", "Northrop"), col=c("green", "blue"), lty=c(2,3),lwd=2, cex=0.8)

##--------Showing what simulated data looks like------------
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')
set.seed(12345)
dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data, breaks=seq(0.5,7, by=0.2), main="Case 1", xlab='')
abline(v=1.0, col="red", lwd=2)
set.seed(12345)
dat1 <- runif(80, 0.5, 1.0)
dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data,breaks=seq(0.5,7, by=0.2), main="Case 2", xlab='')
abline(v=1.0, col="red", lwd=2)
set.seed(12345)
dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=-0.05, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
hist(data,breaks=seq(0.5,7, by=0.2), main="Case 3",xlab='')
abline(v=1.0, col="red", lwd=2)
set.seed(12345)
data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
cens_thr<-1*rbeta(length(data_all),1,0.5)
keep <- data_all>cens_thr
data_keep <- data_all[keep]
data <- sample(data_keep, 1000, replace = FALSE)
hist(data,breaks=seq(0,6.5, by=0.2), main="Case 4", xlab='')
abline(v=1.0, col="red", lwd=2)
head(mythrI2.2)


library(Metrics)

variance <- function(x){
  n <- length(x)
  x_bar <- sum(x)/n
  diff <- x - x_bar
  sdiff <- diff^2
  vary<- sum(sdiff)/n
  return(vary)
}

#Bias and variance of methods by case
#Case 1
pred=dan2001C$thr
(rmse <- rmse(1.0, pred))
(bias <- bias(1.0, pred))
(var <- variance(pred))

#Case 1
diff_mythrI <- mythrI$thr-1
best80 <- diff_mythrI[abs(diff_mythrI) <= quantile(abs(diff_mythrI), 0.8)]
rmse(0,best80)

err <- is.na(wadsthrI$thr)
wadsthrI$thr[err] <- 100
diff_wadsthrI <- wadsthrI$thr-1
best80 <- diff_wadsthrI[abs(diff_wadsthrI) <= quantile(abs(diff_wadsthrI), 0.8)]
rmse(0,best80)

diff_norththrI <- norththrI$thr-1
best80 <- diff_norththrI[abs(diff_norththrI) <= quantile(abs(diff_norththrI), 0.8)]
rmse(0,best80)

#Case 2
diff_mythrI1 <- mythrI1$thr-1
best80 <- diff_mythrI1[abs(diff_mythrI1) <= quantile(abs(diff_mythrI1), 0.8)]
rmse(0,best80)

err <- is.na(wadsthrI1$thr) #This is only 71.6% of threshold choices so no need to adjust to best 80% 
rmse(1,wadsthrI1$thr[!err])

diff_norththrI1 <- norththrI1$thr-1
best80 <- diff_norththrI1[abs(diff_norththrI1) <= quantile(abs(diff_norththrI1), 0.8)]
rmse(0,best80)

#Case 3
diff_mythrI2 <- mythrI2$thr-1
best80 <- diff_mythrI2[abs(diff_mythrI2) <= quantile(abs(diff_mythrI2), 0.8)]
rmse(0,best80)

err <- is.na(wadsthrI2$thr)
wadsthrI2$thr[err] <- 100
diff_wadsthrI2 <- wadsthrI2$thr-1
best80 <- diff_wadsthrI2[abs(diff_wadsthrI2) <= quantile(abs(diff_wadsthrI2), 0.8)]
rmse(0,best80)

diff_norththrI2 <- norththrI2$thr-1
best80 <- diff_norththrI2[abs(diff_norththrI2) <= quantile(abs(diff_norththrI2), 0.8)]
rmse(0,best80)

#Case 4
diff_mythrC <- mythrC$thr-1
best80 <- diff_mythrC[abs(diff_mythrC) <= quantile(abs(diff_mythrC), 0.8)]
rmse(0,best80)

err <- is.na(wadsthrC$thr)
wadsthrC$thr[err] <- 100
diff_wadsthrC <- wadsthrC$thr-1
best80 <- diff_wadsthrC[abs(diff_wadsthrC) <= quantile(abs(diff_wadsthrC), 0.8)]
rmse(0,best80)

diff_norththrC <- norththrC$thr-1
best80 <- diff_norththrC[abs(diff_norththrC) <= quantile(abs(diff_norththrC), 0.8)]
rmse(0,best80)
