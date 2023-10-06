library(Metrics)

#Data

mythrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI.rds")
mythrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI1.rds")
mythrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/mythrI2.rds")
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


mythrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/mythrG_large_sample.rds")
norththrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/norththrG_large_sample.rds")
wadsthrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/wadsthrG_large_sample.rds")
#------------------True quantiles---------------------

#Ideal case , p exceedance prob, par (sig, xi), u true
ideal_quant <- function(p,par,u){
  x_p <- (par[1]/par[2])*((6*p/5)^(-par[2])-1) + u
  return(x_p)
}

#Censored case, p exceedance prob, p_1 P(<1), par (sig, xi), u true
cens_quant <- function(p, par, p_1, u){
  x_p <- ((par[1]+par[2]*(u))/par[2])*((p/(1-p_1))^(-par[2])-1) + u
  return(x_p)  
}

#------------------Estimated quantiles------------------
estimated_quantile <- function(df,p,n){
  lu <- df$len/n
  est_quants <- df$thr + (df$scale/df$shape)*((p/lu)^(-df$shape)-1)
  return(est_quants)
}

#Ideal case
n <- 1200
t <- c(1,10,100)
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- wads_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythrI <- estimated_quantile(df=mythrI, p=p, n=n)
  est_northI <- estimated_quantile(df=norththrI, p=p, n=n)
  err<-is.na(wadsthrI$thr)
  est_wadsI <- estimated_quantile(df=wadsthrI[!err,], p=p, n=n)
  true_I <- ideal_quant(p=p, par=c(0.5,0.1), u=1.0)
  
  my_RMSE[i] <- rmse(true_I, est_mythrI)
  north_RMSE[i] <- rmse(true_I, est_northI)
  wads_RMSE[i] <- rmse(true_I, est_wadsI)
}

RMSE_I <- data.frame(p=p_seq, my=my_RMSE, wads=wads_RMSE, north=north_RMSE,  t=t)

n <- 1200
t <- c(1,10,100)
p_seq <- c(1/n,1/(10*n), 1/(100*n))
dan01_RMSE <- dan19_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythrI <- estimated_quantile(df=dan2001I, p=p, n=n)
  est_northI <- estimated_quantile(df=danthrI, p=p, n=n)
  true_I <- ideal_quant(p=p, par=c(0.5,0.1), u=1.0)
  
  dan01_RMSE[i] <- rmse(true_I, est_mythrI)
  dan19_RMSE[i] <- rmse(true_I, est_northI)
}

RMSE_I <- data.frame(p=p_seq, dan01=dan01_RMSE, dan19=dan19_RMSE,  t=t)

#I1 (I with small sample size)
n <- 480
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_dan2001I1 <- estimated_quantile(df=dan2001I1, p=p, n=n)
  est_northI1 <- estimated_quantile(df=danthrI1, p=p, n=n)
  true_I1 <- ideal_quant(p=p, par=c(0.5,0.1), u=1.0)

  my_RMSE[i] <- rmse(true_I1, est_dan2001I1)
  north_RMSE[i] <- rmse(true_I1, est_northI1)
}

RMSE_I1 <- data.frame(p=p_seq, dan01=my_RMSE, dan19=north_RMSE, t=t)


#I2 (I with, n=2400 shape=-0.05)
n <- 2400
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_dan2001I2 <- estimated_quantile(df=dan2001I2, p=p, n=n)
  est_northI2 <- estimated_quantile(df=danthrI2, p=p, n=n)
  true_I2 <- ideal_quant(p=p, par=c(0.5,-0.05), u=1.0)
  
  my_RMSE[i] <- rmse(true_I2, est_dan2001I2)
  north_RMSE[i] <- rmse(true_I2, est_northI2)
}

RMSE_I2 <- data.frame(p=p_seq, dan01=my_RMSE,  dan19=north_RMSE, t=t)

#C 
n <- 1000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_dan2001C <- estimated_quantile(df=dan2001C, p=p, n=n)
  est_northC <- estimated_quantile(df=danthrC, p=p, n=n)
  true_C <- cens_quant(p=p, par=c(0.5,0.1), p_1= 0.5290265, u=1.0)

  my_RMSE[i] <- rmse(true_C, est_dan2001C)
  north_RMSE[i] <- rmse(true_C, est_northC)
}

RMSE_C <- data.frame(p=p_seq, dan01=my_RMSE, dan19=north_RMSE,  t=t)


#G
n <- 20000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythrG <- estimated_quantile(df=norththrG_large_sample, p=p, n=n)
  true_G <- qnorm(1-p)

  my_RMSE[i] <- rmse(true_G, est_mythrG) 
}
my_RMSE

#Bias and variance breakdown
variance <- function(x){
  n <- length(x)
  x_bar <- sum(x)/n
  diff <- x - x_bar
  sdiff <- diff^2
  vary<- sum(sdiff)/n
  return(vary)
}

p <- 1/(100*n)
est_G1 <- estimated_quantile(df=mythrG, p=p, n=n)
err<-is.na(wadsthrG$thr)
est_G2 <- estimated_quantile(df=wadsthrG[!err,], p=p, n=n)
est_G3 <- estimated_quantile(df=norththrG, p=p, n=n)
est_G4 <- estimated_quantile(df=dan2001G, p=p, n=n)
est_G5 <- estimated_quantile(df=danthrG, p=p, n=n)
true <- qnorm(1-p)

bias(true,est_G1)
bias(true, est_G2)
bias(true, est_G3)
bias(true, est_G4)
bias(true, est_G5)

variance(est_G1)
variance(est_G2)
variance(est_G3)
variance(est_G4)
variance(est_G5)

RMSE_G <- data.frame(p=p_seq, my=my_RMSE, wads=wads_RMSE, north=north_RMSE,  t=t)
RMSE_C

err<-is.na(wadsthrG$thr)
length(wadsthrG$thr[err])

#Endpoint estimation
estimated_end <- function(df){
  end <- df$thr - df$scale/df$shape
  return(end)
}

end_mythr <- estimated_end(mythr3A)
end_north <- estimated_end(norththr3A)
end_mythr
end_north
length(norththr3A$thr[norththr3A$shape>0])


dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
plot(density(end_mythr),  lwd=2, col="green", lty=2, main="Case 3A", ylab="Density", xlab="Thresholds")
lines(density(end_north), lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Northrop"), col=c("green", "blue"), lty=c(2,3),lwd=2, cex=0.8)
abline(v=3.5, col="red", lwd=1.5)
plot(density(mythr3B$thr),  lwd=2, col="green", lty=2, xlim=c(0.5,2.5), main="Case 3B", ylab="Density", xlab="Thresholds")
lines(density(norththr3B$thr), lwd=2, col="blue", lty=3)
legend("topright", c("Our method", "Northrop"), col=c("green", "blue"), lty=c(2,3),lwd=2, cex=0.8)

RMSE_I
