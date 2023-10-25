library(Metrics)

#Data

#Case 1
mythrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/mythrI.rds")
wadsthrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/wadsthrI.rds")
norththrI <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/norththrI.rds")

mythrI_large <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/mythrI_large.rds")
wadsthrI_large <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/wadsthrI_large.rds")
norththrI_large <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/norththrI_large.rds")

wadsthrI_large_halfpercent <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/wadsthrI_large_halfpercent.rds")
norththrI_large_halfpercent <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/norththrI_large_halfpercent.rds")
mythrI_large_halfpercent <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/mythrI_large_halfpercent.rds")

#Case 2
wadsthrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/wadsthrI1.rds")
norththrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/norththrI1.rds")
mythrI1 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/mythrI1.rds")

#Case 3
mythrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/mythrI2.rds")
wadsthrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/wadsthrI2.rds")
norththrI2 <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/norththrI2.rds")

#Case 4
mythrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/mythrC.rds")
wadsthrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/wadsthrC.rds")
norththrC <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimIC/Seeded_final/Rerun with quantiles/norththrC.rds")

#Gaussian case 
mythrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Original/mythrG.rds")
wadsthrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Original/wadsthrG.rds")
norththrG <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Original/norththrG.rds")

#with large sample size (20000)
mythrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/mythrG_large_sample_halfpercentgrid.rds")
norththrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/norththrG_large_sample_halfpercentgrid.rds")
wadsthrG_large_sample <- readRDS("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Constant Threshold Selection/SimGauss/Seeded_final/Larger sample size/wadsthrG_large_sample_halfpercentgrid.rds")

#RMSE of threshold choice
#Case 1
rmse(mythrI$thr,1)
rmse(wadsthrI$thr[!is.na(wadsthrI$thr)],1)
rmse(norththrI$thr, 1)

#Case 1 large sample
rmse(mythrI_large_sample$thr,1)
rmse(wadsthrI_large_sample$thr[!is.na(wadsthrI_large_sample$thr)],1)
rmse(norththrI_large_sample$thr, 1)

#Case 2
rmse(mythrI1$thr,1)
rmse(wadsthrI1$thr[!is.na(wadsthrI1$thr)],1)
rmse(norththrI1$thr, 1)

#Case 3
rmse(mythrI2$thr,1)
rmse(wadsthrI2$thr[!is.na(wadsthrI2$thr)],1)
rmse(norththrI2$thr, 1)

#Case 4
rmse(mythrC$thr,1)
rmse(wadsthrC$thr[!is.na(wadsthrC$thr)],1)
rmse(norththrC$thr, 1)

#Percentage fail of Wadsworth in Case 1-4
length(wadsthrI$thr[is.na(wadsthrI$thr)])/500
length(wadsthrI1$thr[is.na(wadsthrI1$thr)])/500
length(wadsthrI2$thr[is.na(wadsthrI2$thr)])/500
length(wadsthrC$thr[is.na(wadsthrC$thr)])/500

length(wadsthrI_large_sample$thr[is.na(wadsthrI_large_sample$thr)])/500

#RMSE on samples where Wadsworth was successful only
err <- is.na(wadsthrC$thr)
rmse(mythrC$thr[!err],1)
rmse(wadsthrC$thr[!err],1)
rmse(norththrC$thr[!err], 1)

#Bias of threshold choice
#Case 1
bias(mythrI$thr,1)
bias(wadsthrI$thr[!is.na(wadsthrI$thr)],1)
bias(norththrI$thr, 1)

#Case 2
bias(mythrI1$thr,1)
bias(wadsthrI1$thr[!is.na(wadsthrI1$thr)],1)
bias(norththrI1$thr, 1)

#Case 3
bias(mythrI2$thr,1)
bias(wadsthrI2$thr[!is.na(wadsthrI2$thr)],1)
bias(norththrI2$thr, 1)

#Case 4
bias(mythrC$thr,1)
bias(wadsthrC$thr[!is.na(wadsthrC$thr)],1)
bias(norththrC$thr, 1)

#Case 1
variance(mythrI$thr )
variance(wadsthrI$thr[!is.na(wadsthrI$thr)] )
variance(norththrI$thr )

#Case 2
variance(mythrI1$thr )
variance(wadsthrI1$thr[!is.na(wadsthrI1$thr)] )
variance(norththrI1$thr )

#Case 3
variance(mythrI2$thr )
variance(wadsthrI2$thr[!is.na(wadsthrI2$thr)] )
variance(norththrI2$thr )

#Case 4
variance(mythrC$thr )
variance(wadsthrC$thr[!is.na(wadsthrC$thr)] )
variance(norththrC$thr )

#Variance of threshold choice
variance <- function(x){
  n <- length(x)
  x_bar <- sum(x)/n
  diff <- x - x_bar
  sdiff <- diff^2
  vary<- sum(sdiff)/n
  return(vary)
}

#Danielsson methods

rmse(dan2001thrC$thr, 1)
bias(dan2001thrC$thr, 1)
variance(dan2001thrC$thr)

rmse(dan2019thrC$thr, 1)
bias(dan2019thrC$thr, 1)
variance(dan2019thrC$thr)

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

#---------RMSE of estimated quantiles for Case 1-4--------------------
#Case 1-4
#Edit sample size, estimated threshold vectors and true parameters for case 1-4
n <- 20000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- wads_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythr <- estimated_quantile(df=mythrI_large_halfpercent,p=p, n=n)
  est_north <- estimated_quantile(df=norththrI_large_halfpercent, p=p, n=n)
  err<-is.na(wadsthrI_large_halfpercent$thr)
  est_wads <- estimated_quantile(df=wadsthrI_large_halfpercent[!err,], p=p, n=n)
  true <- ideal_quant(p, par=c(0.5,0.1), u=1.0)
  
  my_RMSE[i] <- rmse(est_mythr, true)
  north_RMSE[i] <- rmse(est_north, true)
  wads_RMSE[i] <- rmse(est_wads, true)
}

(RMSE <- data.frame(p=p_seq, my=my_RMSE, wads=wads_RMSE, north=north_RMSE))#,  t=t))

#Danielsson methods Case 1-4
n <- 1000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
dan01_RMSE <- dan19_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythrI <- estimated_quantile(df=dan2001thrI, p=p, n=n)
  est_northI <- estimated_quantile(df=dan2019thrI, p=p, n=n)
  true_I <- cens_quant(p=p, par=c(0.5,0.1),  p_1= 0.5290265, u=1.0 )
  
  dan01_RMSE[i] <- rmse(est_mythrI,true_I)
  dan19_RMSE[i] <- rmse(est_northI,true_I)
}

(RMSE_I <- data.frame(p=p_seq, dan01=dan01_RMSE, dan19=dan19_RMSE))

#C 
n <- 1000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- north_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_dan2001C <- estimated_quantile(df=dan2001C, p=p, n=n)
  est_northC <- estimated_quantile(df=danthrC, p=p, n=n)
  true_C <- cens_quant(p=p, par=c(0.5,0.1), p_1= 0.5290265, u=1.0)

  my_RMSE[i] <- rmse(est_dan2001C , true_C)
  north_RMSE[i] <- rmse(est_northC, true_C)
}

RMSE_C <- data.frame(p=p_seq, dan01=my_RMSE, dan19=north_RMSE,  t=t)

#Gaussian samples
n <- 2000
p_seq <- c(1/n,1/(10*n), 1/(100*n))
my_RMSE <- numeric(3)
for(i in 1:3){
  p <- p_seq[i]
  est_mythrG <- estimated_quantile(df=dan2019thrG, p=p, n=n)
  true_G <- qnorm(1-p)

  my_RMSE[i] <- variance(est_mythrG)#, true_G) 
}
my_RMSE

