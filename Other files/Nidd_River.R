#Nidd River data
library(evir)
data("nidd.thresh")
plot(nidd.thresh)
plot(niddnew)
plot(nidd1966)

#Parstab
# shapestab <- function (data, thresholds, reverse = TRUE,#models = 30, start = 15, end = 500,, 
#                        ci = 0.95, auto.scale = TRUE, labels = TRUE) 
# {
#   data <- as.numeric(data)
#   qq <- 0
#   if (ci) 
#     qq <- qnorm(1 - (1 - ci)/2)
#   x <- thresholds #trunc(seq(from = length(data[threshold, to = start, 
#   #length = models))
#   gpd.dummy <- function(thr, data) {
#     out <- gpd(data = data, threshold = thr, information = "expected")
#     c(out$n.exceed, out$par.ests[1], out$par.ses[1])
#   }
#   mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
#   mat <- rbind(mat, x)
#   dimnames(mat) <- list(c("exceedances", "shape", 
#                           "se", "thresholds"), NULL)
#   thresh <- mat[4, ]
#   exceed <- mat[1, ]
#   y <- mat[2, ]
#   yrange <- range(y)
#   if (ci) {
#     u <- y + mat[3, ] * qq
#     l <- y - mat[3, ] * qq
#     yrange <- range(y, u, l)
#   }
#   index <- x
#   if (reverse) 
#     index <- -x
#   if (auto.scale) 
#     plot(index, y, ylim = yrange, type = "l", xlab = "", 
#          ylab = "", axes = FALSE)
#   else plot(index, y, type = "l", xlab = "", ylab = "",axes = FALSE)
#   abline(v=1.0, lwd=1.5, lty=2, col="red")
#   axis(1, at = index, labels = paste(round(thresholds,2)), tick = T)
#   axis(2)
#   axis(3, at = index, labels = paste(seq(0,0.95,by=0.05)), tick = T)
#   box()
#   if (ci) {
#     lines(index, u, lty = 2, col = 2)
#     lines(index, l, lty = 2, col = 2)
#   }
#   if (labels) {
#     labely <- expression(paste("",hat(xi)))
#     # if (ci) 
#     #   labely <- paste(labely, " (CI, p = ", ci, ")", 
#     #                   sep = "")
#     title(xlab = "Threshold", ylab = labely)
#     mtext("Quantile", side = 3, line = 3)
#   }
#   invisible(mat)
# }

#My method
source("thr_selection_final.R")

#Wadsworth
source("JointMLEFunctions.r")

#Northrop
library(threshr)

# thresholds <- quantile(nidd.thresh, c(0, 0.3, 0.6, 0.9), names=F)
# thresholds
# wadsres <- NHPP.diag(nidd.thresh, u=thresholds, plot.out=FALSE, UseQuantiles = FALSE)
# wadsres
# 
# thresholds <- quantile(nidd.thresh, seq(0,0.95, by=0.05), names=F)
# mythr <- thr_dist(nidd.thresh, thresholds)
# mythr
# thresholds
# norththr <- summary(ithresh(nidd.thresh, u_vec = thresholds))
# norththr
# length(nidd.thresh[nidd.thresh>109.076])
# opt_north <- optim(GPD_LL, par=c(mean(excess), 0.1), z=excess, control = list(fnscale=-1))
# dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
# par(mfrow=c(1,3),bg='transparent')
# thr <- 149.09
# excess <- nidd.thresh[nidd.thresh>thr] - thr 
# opt_north <- optim(GPD_LL, par=c(mean(excess), 0.1), z=excess, control = list(fnscale=-1))
# n <- length(excess)
# probs=c(1:n)/(n+1)
# samp_q <- quantile(excess, probs, names = F )
# model_q <- qgpd(probs, shape=opt_north$par[2], scale=opt_north$par[1])
# plot(samp_q, model_q, ylab="Model Quantiles", xlab = "Sample Quantiles", main="Wadsworth")
# abline(a=0,b=1)

# 
# #Algorithm 2: threshold uncertainty
# 
# data <- nidd.thresh
# thresholds <- quantile(data, seq(0,0.95, by=0.05), names=F)
# (mythr <- thr_dist(data, thresh = thresholds))
# lambda_1 <- length(data[data>mythr$thresh])/n
# point_est <- mythr$thresh + (mythr$par[1]/mythr$par[2])*((p/lambda_1)^(-mythr$par[2])-1)
# n <- length(data)
# m <- 200 #no. of bootstraps
# p <- 1/(10*n)
# boot_thr_ests <- boot_thr <- numeric(m)
# for(i in 1:m){
#   data_b <- sample(data,n, replace=TRUE)
#   mythr_b <- thr_dist(data_b, thresh = thresholds)
#   lambda <- length(data_b[data_b>mythr_b$thresh])/n
#   boot_thr[i] <- mythr_b$thresh
#   boot_thr_ests[i] <- mythr_b$thresh + (mythr_b$par[1]/mythr_b$par[2])*((p/lambda)^(-mythr_b$par[2])-1)
# }
# 
# #Confidence interval on threshold
# boot_thr_sorted <- sort(boot_thr)
# thr_l <- boot_thr_sorted[0.025*m]
# thr_u <- boot_thr_sorted[0.975*m]
# (ci_thr <- c(thr_l, thr_u))
# 
# hist(boot_thr_sorted)
# boot_q_sorted <- sort(boot_thr_ests)
# q_l <- boot_q_sorted[0.025*m]
# q_u <- boot_q_sorted[0.975*m]
# (ci_q <- c(q_l, q_u))
# point_est
# max(nidd.thresh)
# sort(nidd.thresh)
# hist(boot_q_sorted)

#Algorithm 1: par uncertainty
data <- nidd.thresh
u <- 69.74
n <- length(data)
m <- 400 #no. of bootstraps
t <- c(1,2,5,10,25,50,100)
excess <- nidd.thresh[nidd.thresh>u] - u
n_u <- length(excess)
opt <- optim(GPD_LL, par=c(0.1, 0.1), z=excess, control = list(fnscale=-1) )
point_est_old <- u + (opt$par[1]/opt$par[2])*(((t*n_u)/35)^(opt$par[2])-1)
boot_ests <- matrix(NA, nrow = m, ncol=t)
for(i in 1:m){
    excess_b <- rgpd(n_u, scale = opt$par[1], shape= opt$par[2]) #parametric bootstrap
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    boot_ests[i,] <- u + (opt_b$par[1]/opt_b$par[2])*(((t*n_u)/35)^(opt_b$par[2])-1)
}

#Algorithm 3: parameter and threshold uncertianty
data <- nidd.thresh
thresholds <- quantile(data, seq(0,0.95, by=0.05), names=F)
(mythr <- thr_dist(data, thresh = thresholds))
n <- length(data)
m <- 50 #no. of bootstraps
p_seq <- c(1/n, 1/(5*n), 1/(10*n), 1/(20*n), 1/(50*n), 1/(100*n))
p_len <- length(p_seq)
point_est <- numeric(p_len)
lambda_1 <- length(data[data>mythr$thresh])/n
boot_thr_ests <- matrix(NA, nrow = m, ncol=p_len)
boot_thr <- numeric(m)
for(j in 1:p_len){
  p <- p_seq[j]
  point_est[j] <- mythr$thresh + (mythr$par[1]/mythr$par[2])*((p/lambda_1)^(-mythr$par[2])-1)
  for(i in 1:m){
    data_b <- sample(data,n, replace=TRUE) #non-parametric bootstrap
    mythr_b <- thr_dist(data_b, thresh = thresholds)
    boot_thr[i] <- mythr_b$thresh
    n_b <- length(data_b[data_b>mythr_b$thresh])
    lambda <- n_b/n
    excess_b <- rgpd(n_b, scale = mythr_b$par[1], shape= mythr_b$par[2]) #parametric bootstrap
    opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
    boot_thr_ests[i,j] <- boot_thr[i] + (opt_b$par[1]/opt_b$par[2])*((p/lambda)^(-opt_b$par[2])-1)
  }
}
#Confidence interval on threshold
boot_thr_sorted <- sort(boot_thr)
thr_l <- boot_thr_sorted[floor(0.025*m)]
thr_u <- boot_thr_sorted[ceiling(0.975*m)]
(ci_thr <- c(thr_l, thr_u))

#Confidence bounds for return levels
m <- 400
boot_q_sorted <- apply(boot_ests_1, 2, sort)
q_l <- boot_q_sorted[floor(0.025*m),]
q_u <- boot_q_sorted[ceiling(0.975*m),]
plot((1-p_seq), point_est, ylim=c(100, 1500), type='l', ylab="Return level", xlab="Non-exceedance probability")
lines((1-p_seq), q_l, col="blue", lty="dashed")
lines((1-p_seq), q_u, col="blue", lty="dashed")

saveRDS(nidd.thresh, "Nidddata.rds")


#Analysis of uncertainty sim

#Read in files from Nidd river -> return levels -> 2

b1_sort <- apply(boot_ests_1, 2, sort)
b3_sort <- apply(boot_ests_3, 2, sort )
     
#Plot of new data return levels using all the data 
ind <- c(1,2,5,10,25,50,100)
plot(ind, point_est,ylim = c(0,800), ylab="Return Levels", xlab="Return period", pch=19)
ql1 <- b1_sort[400*0.025, ]
lines(ind, ql1[1:7])
qu1 <- b1_sort[400*0.975, ]
lines(ind, qu1[1:7])
ql2 <- b3_sort[160000*0.025,]
lines(ind, ql2, col="black")
qu2 <- b3_sort[160000*0.975,]
lines(ind, qu2, col="black")

#------------------------------------------------------------------------------------------------------------

#Fine grid of thresholds for Old Nidd data

thresholds <- unique(sort(nidd.thresh)[-c(144:length(nidd.thresh))])
thresholds <- quantile(nidd.thresh, seq(0,0.95,by=0.01), names=F)
(thr_fine <- thr_dist(nidd.thresh, thresh = thresholds, k=40))
#Comparison with new River Nidd data

niddnew <- read.delim("Nidd River Case study/NRFAPeakFlow_v11/suitable-for-pooling/027001new.pt.txt", header=F, sep = ',')
colnames(niddnew) <- c("Date", "Flow", "V3")
niddnew <- niddnew[-length(niddnew$Date),]

plot(niddnew$Flow, ylab="Flow rate")
#CIs on old point estimates (Algorithm 1) using new data from 1934-69
data <- niddnew$Flow[1:136]
data <- nidd.thresh
plot(nidd.thresh, ylab="Flow rate")
thresh <- quantile(data, seq(0,0.95, by=0.05), names=F)
(mythr <- thr_dist(data, thresh=thresh))
u <-mythr$thresh
u <- 69.74
m <- 400 #no. of bootstraps
t <- c(1,2,5,10,25,50,100)
excess <- data[data>u] - u
n_u <- length(excess)
opt <- optim(GPD_LL, par=c(0.1, 0.1), z=excess, control = list(fnscale=-1))
point_est_old <- u + (opt$par[1]/opt$par[2])*(((t*n_u)/35)^(opt$par[2])-1)
points(t, point_est_old,ylim=c(0,1000), pch=4, col="red")
boot_ests_old <- matrix(NA, nrow = m, ncol=length(t))
for(i in 1:m){
  excess_b <- rgpd(n_u, scale = opt$par[1], shape= opt$par[2]) #parametric bootstrap
  opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
  boot_ests_old[i,] <- u + (opt_b$par[1]/opt_b$par[2])*(((t*n_u)/35)^(opt_b$par[2])-1)
}

boot_old_sort <- apply(boot_ests_old, 2, sort)
q_l_old <-  boot_old_sort[400*0.025, ]
lines(ind, q_l_old, col="red")
q_u_old <- boot_old_sort[400*0.975, ]
lines(ind, q_u_old, col="red")

u <-109.08
(norththr <- summary(ithresh(data, thresh))[3])
u <- 57.31
m <- 400 #no. of bootstraps
t <- c(1,2,5,10,25,50,100)
excess <- data[data>u] - u
n_u <- length(excess)
opt <- optim(GPD_LL, par=c(0.1, 0.1), z=excess, control = list(fnscale=-1) )
point_est_old <- u + (opt$par[1]/opt$par[2])*(((t*n_u)/35)^(opt$par[2])-1)
points(t, point_est_old, pch=15, col="blue")
boot_ests_old <- matrix(NA, nrow = m, ncol=length(t))
for(i in 1:m){
  excess_b <- rgpd(n_u, scale = opt$par[1], shape= opt$par[2]) #parametric bootstrap
  opt_b <- optim(GPD_LL, par=c(mean(excess_b), 0.1), z=excess_b, control = list(fnscale=-1))
  boot_ests_old[i,] <- u + (opt_b$par[1]/opt_b$par[2])*(((t*n_u)/35)^(opt_b$par[2])-1)
}

boot_old_sort <- apply(boot_ests_old, 2, sort)
q_l_old <-  boot_old_sort[400*0.025, ]
lines(ind, q_l_old, col="blue")
q_u_old <- boot_old_sort[400*0.975, ]
lines(ind, q_u_old, col="blue")
#Check for seasonality in new data to apply to Model selection

#Sorting dates

library(lubridate)
niddnew$Date <- dmy(niddnew$Date)

nidd1966 <- niddnew$Flow[109:length(niddnew$Flow)]

thresholds <- quantile(nidd1966, seq(0,0.95, by=0.01), names = F)
(thr66 <- thr_dist(nidd1966, thresh = thresholds)) 
t <- c(1,2,5,10,25,50,100)
u <- thr66$thresh
n_u <- length(nidd6698[nidd6698>u])
point_est_1966 <- u + (thr66$par[1]/thr66$par[2])*(((t*n_u)/55)^(thr66$par[2])-1)
points(t, point_est_1966, pch=4, col="red")

#Splitting data into new block of 32 years from 1966 and comparing to using all data beyond 1966 

nidd6698 <- niddnew$Flow[109:252]
thresh <- quantile(nidd6698, seq(0, 0.95, by=0.01), names=F)
(mythr <- thr_dist(nidd6698, thresh = thresh))
t <- c(1,2,5,10,25,50,100)
u <- mythr$thresh
n_u <- length(nidd6698[nidd6698>u])
point_est_6698 <- u + (mythr$par[1]/mythr$par[2])*(((t*n_u)/32)^(mythr$par[2])-1)
plot(t, point_est_6698,ylim = c(90,270), main="", xlab= "Return period", ylab = "Return level")
