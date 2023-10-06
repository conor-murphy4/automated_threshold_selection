source('threshold_paper_code-main/00_src/_gpd.R')

#install.packages("evir")

library(evir)

reps <- 30

estimated_quantile <- function(p,xis,sigmas,lens,mus,n){
  est_quants <- NULL
  #est_quants_2 <- NULL
  for(i in 1:length(xis)){
    lu <- lens[i]/n[i]
    q_est <- qgpd(1-(p/lu),shape = xis[i],scale = sigmas[i],mu=mus[i])
    #q_est_2 <- (betas[i]/xis[i])*((lu/p)^(xis[i])-1) + mus[i]
    est_quants <- c(est_quants,q_est)
    #est_quants_2 <- c(est_quants_2,q_est_2)
  }
  #est <- data.frame(x=est_quants,y=est_quants_2)
  return(est_quants)
}

load("C:/Users/murphyc4/OneDrive - Lancaster University/Results.Rda")
p=1/100000
Est_quantiles_o <- estimated_quantile(p,xis=Results.df$XI_norm,sigmas = Results.df$SIG_norm, lens = Results.df$LEN_norm, mus=Results.df$GPD_norm, n=Results.df$N)
Est_quantiles_e <- estimated_quantile(p,xis=Results.df$XI_enorm,sigmas = Results.df$SIG_enorm, lens = Results.df$LEN_enorm, mus=Results.df$EXP_norm, n=Results.df$N)
True_quantiles <- qnorm(1-p)

library(Metrics)#RMSEs
RMSEs_o <- tapply(Results.df$GPD, rep(1:(length(Results.df$GPD)/reps), each=reps), rmse, actual=rep(1.0,reps))
RMSEs_e <- tapply(Results.df$EXP, rep(1:(length(Results.df$EXP)/reps), each=reps), rmse, actual=rep(1.0,reps))
RMSEs_o_norm <- tapply(Est_quantiles_o, rep(1:(length(Results.df$GPD_norm)/reps), each=reps), rmse, True_quantiles)
RMSEs_e_norm <- tapply(Est_quantiles_e, rep(1:(length(Results.df$EXP_norm)/reps), each=reps), rmse, True_quantiles)

#Absolute error of estimates

sample_sizes <- c(seq(1000,10000,by=1000),seq(12000,20000, by=2000), 30000, 40000, 50000)
load("C:/Users/murphyc4/OneDrive - Lancaster University/Results.Rda")

error <- abs(Results.df$GPD-1.0) - abs(Results.df$EXP-1.0)
length(error[error<0])/length(error)
length(error[error>0])/length(error)
length(error[error==0])/length(error)

DF <- data.frame(sample_sizes,error)
dev.new(width=10.5, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')
sizeplot(sample_sizes, error, xlab='', ylab='',lwd=2,cex.axis=1.5)
points(DF$sample_sizes[which(DF$error<0)],DF$error[which(DF$error<0)],col='green',pch=19,cex=1.3)
points(DF$sample_sizes[which(DF$error>0)],DF$error[which(DF$error>0)],col='red',pch=19,cex=1.3)
points(5000,-0.05,col='green',pch=19,cex=2.3)
points(7000,-0.05,col='green',pch=19,cex=2.0)
points(9000,-0.2,col='green',pch=19,cex=1.7)
points(14000,-0.05,col='green',pch=19,cex=1.7)
points(14000,-0.1,col='green',pch=19,cex=1.7)
points(14000,-0.15,col='green',pch=19,cex=1.7)
points(50000,-0.05,col='green',pch=19,cex=1.7)
abline(h=0, lty="dashed", lwd=2)
mtext("EQD better", at=25000, padj=15,cex=1.5,col='green')
mtext("Varty better", at=25000, padj=2,cex=1.5,col='red')
mtext("No difference", at=29000, padj=7, cex=1.2, col="black")
mtext('Sample Sizes',side=1, cex=1.5,line=2.5)
mtext("Difference in Absolute Error", side=2, cex=1.5, line=2.5)

plot(sample_sizes,RMSEs_o_norm, pch=19,type='l', xlab= "", ylab="",lwd=2,cex.axis=1.5)
lines(sample_sizes, RMSEs_e_norm, col="red", lwd=2)
mtext('Sample Sizes',side=1, cex=1.5,line=2.5)
mtext("RMSE", side=2, cex=1.5, line=2.5)
legend("topright", c("EQD","Varty method"),col=c("black", "red"),lwd=c(2,2),lty=c(1,1) , cex=1.3)

###-------------------------
par(mfrow=c(2,1))
plot(Results.df$GPD_norm)
points(Results.df$EXP_norm, col="red")

library(Metrics)

variance <- function(x){
  n <- length(x)
  x_bar <- sum(x)/n
  diff <- x - x_bar
  sdiff <- diff^2
  vary<- sum(sdiff)/n
  return(vary)
}
#Bias
bias_GPD <- tapply(Results.df$GPD, rep(1:(length(Results.df$GPD)/reps), each=reps), bias, actual=1.0)
bias_EXP <- tapply(Results.df$EXP, rep(1:(length(Results.df$EXP)/reps), each=reps), bias, actual=1.0)
bias_GPDnorm <- tapply(Est_quantiles_o,rep(1:(length(Est_quantiles_o)/reps), each=reps), bias, actual=True_quantiles )
bias_EXPnorm <- tapply(Est_quantiles_e,rep(1:(length(Est_quantiles_e)/reps), each=reps), bias, actual=True_quantiles )

#Variance
var_GPD <- tapply(Results.df$GPD, rep(1:(length(Results.df$GPD)/reps), each=reps), variance)
var_EXP <- tapply(Results.df$EXP, rep(1:(length(Results.df$EXP)/reps), each=reps), variance)
var_GPDnorm <- tapply(Est_quantiles_o,rep(1:(length(Est_quantiles_o)/reps), each=reps), variance)
var_EXPnorm <- tapply(Est_quantiles_e,rep(1:(length(Est_quantiles_e)/reps), each=reps), variance)


#MSE
mse_GPD <- tapply(Results.df$GPD, rep(1:(length(Results.df$GPD)/reps), each=reps), mse, actual=1.0)
mse_EXP <- tapply(Results.df$EXP, rep(1:(length(Results.df$EXP)/reps), each=reps), mse, actual=1.0)
mse_GPDnorm <- tapply(Est_quantiles_o,rep(1:(length(Est_quantiles_o)/reps), each=reps), mse, actual=True_quantiles )
mse_EXPnorm <- tapply(Est_quantiles_e,rep(1:(length(Est_quantiles_e)/reps), each=reps), mse, actual=True_quantiles )

sqrt(mse_GPD)
#Check
MSE <- bias_GPDnorm^2 + var_GPDnorm
mse_GPDnorm


var(Results.df$GPD[300:330])
var(Results.df$EXP[300:330])

