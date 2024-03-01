library(evir)
source("helper_functions.R")
library(mev)

#Code from evir::shape function adjusted to include bootstrapped confidence intervals in parameter stability plot

#Function to produce plot showing how the estimate of the shape parameter varies with threshold/quantile

#Arguments:

# data - numeric vector of data
# thresholds - numeric vector of thresholds
# Q - numeric vector of quantiles corresponding to thresholds
# reverse - should plot be by increasing threshold (TRUE) or number of extremes (FALSE)
# ci - probability for asymptotic confidence band; for no confidence band set to zero
# auto.scale - whether or not plot should be automatically scaled; if not, xlim and ylim graphical parameters may be entered
# labels - whether or not axes should be labelled
# boot - whether or not bootstrapped confidence intervals should be included
# m.boot - number of bootstraps

shapestabboot <- function (data, thresholds,Q, reverse = TRUE, ci = 0.95, auto.scale = TRUE, labels = TRUE, boot=FALSE, m.boot=200){
  data <- as.numeric(data)
  n <- length(data)
  qq <- 0
  if (ci) 
    qq <- qnorm(1 - (1 - ci)/2)
  x <- thresholds 
  gpd.dummy <- function(thr, data) {
    out <- gpd(data = data, threshold = thr, information = "expected")
    c(out$n.exceed, out$par.ests[1],out$par.ests[2], out$par.ses[1])
  }
  mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("exceedances", "shape","scale", "se", "thresholds"), NULL)
  thresh <- mat[5, ]
  exceed <- mat[1, ]
  y <- mat[2, ]
  sc <- mat[3, ]
  yrange <- range(y)
  if (ci) {
    u <- y + mat[4, ] * qq
    l <- y - mat[4, ] * qq
    yrange <- range(y, u, l)
  }
  if(boot){
    u_boot <- l_boot <- numeric(length(thresh))
    for(i in 1:length(thresh)){
      n_boot <- length(data[data>thresh[i]])
      shapes <- numeric(m.boot)
      for(j in 1:m.boot){
        ex <- rgpd(n_boot, shape=y[i], scale=sc[i], mu=0)
        bootfit <- optim(GPD_LL, z=ex, par=c(mean(ex), 0.1), control=list(fnscale=-1))
        shapes[j] <- bootfit$par[2]
      }
      u_boot[i] <- quantile(shapes, 0.975)
      l_boot[i] <- quantile(shapes, 0.025)
    }
    yrange <- range(y, u, l, u_boot, l_boot)
  }
 
  index <- x
  if (reverse) 
    index <- -x
  if (auto.scale) 
    plot(index, y, ylim = yrange, type = "l", xlab = "", 
         ylab = "", axes = FALSE)
  else plot(index, y, type = "l", xlab = "", ylab = "",axes = FALSE)
  axis(1, at = index[seq(1,length(index), by=2)], labels = paste(round(thresholds[seq(1,length(index), by=2)],1)), tick = T)
  axis(2)
  axis(3, at = index[seq(1,length(index), by=2)], labels = paste(Q[seq(1,length(index), by=2)]), tick = T)
  box()
  if (ci) {
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if(boot){ #added by CM to add bootsrapped intervals to plot
    lines(index, u_boot, lty = 3, col = "blue")
    lines(index, l_boot, lty = 3, col = "blue")
  }
  if (labels) {
    labely <- expression(hat(xi))
    title(xlab = "Threshold", ylab = labely, cex.lab=1.1)
    mtext("Quantile", side = 3, line = 3, cex=1.1)
  }
  invisible(mat)
}


# Results for Figure 1 in main text----------------------------------------------------

#Plotting window
dev.new(width=9.17, height=3.7,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

#Generating Case 4 data with sample size of 200 
# set.seed(12345)
# data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
# cens_thr<-1.0*rbeta(length(data_all),1,0.5)
# keep <- data_all>cens_thr
# data_keep <- data_all[keep]
# data <- sample(data_keep, 200, replace = FALSE)

#Candidate threshold grid
Q <-  seq(0,0.95,by=0.05)
thresh <- quantile(data_4, Q, names=FALSE)

#Parameter stability plot
shapestabboot(data_4, thresholds = thresh, Q=Q, reverse=F, boot=TRUE,m.boot=200)
abline(v=1.0, col="green")

#Comparison with profile intervals
# tstab.gpd(xdat=data, thresh = thresh, method = "profile", which = "shape")
# abline(v=1.0, col="green")

####---------Nidd river data----------------------------
data(nidd.thresh)

#Candidate threshold grid
Q <- seq(0,0.95,by=0.05)
thresholds <- quantile(nidd.thresh,Q, names=F)

#Parameter stability plot
shapestabboot(nidd.thresh, thresholds = thresholds[-20], Q=Q, reverse=F, boot=TRUE,m.boot=200)

#Comparison with profile intervals
tstab.gpd(xdat=nidd.thresh, thresh = thresholds[-20], method = "profile", which = "shape")



#Figure S.1: Parameter stability plots on extra simulated examples
par(mfrow=c(1,3))
data_sim_study_Case_1 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/data_sim_study_Case_1.rds")
data_sim_study_Case_2 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/data_sim_study_Case_2.rds")
data_sim_study_Case_3 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/data_sim_study_Case_3.rds")
data_sim_study_Case_4 <- readRDS("//luna/FST/MA/Stor-i/murphyc4/Constant Threshold Selection/Rerun with quantiles/data_sim_study_Case_4.rds")

data_1 <- data_sim_study_Case_1[,1]
data_2 <- data_sim_study_Case_2[,1]
data_3 <- data_sim_study_Case_3[,1]
data_4 <- data_sim_study_Case_4[,1]

#Candidate threshold grid
Q <- seq(0,0.95,by=0.05)
thresholds_1 <- quantile(data_1,Q, names=F)
thresholds_2 <- quantile(data_2,Q, names=F)
thresholds_3 <- quantile(data_3,Q, names=F)
thresholds_4 <- quantile(data_4,Q, names=F)

#Parameter stability plot
shapestabboot(data_1, thresholds = thresholds_1, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")
shapestabboot(data_2, thresholds = thresholds_2, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")
shapestabboot(data_3, thresholds = thresholds_3, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")
shapestabboot(data_4, thresholds = thresholds_4, Q=Q, reverse=F, boot=FALSE,m.boot=200)
abline(v=1, lty="dashed", col="green")


