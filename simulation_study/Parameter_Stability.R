library(evir)
source("helper_functions.R")
source("_gpd.R")

shapestab <- function (data, thresholds, reverse = TRUE, ci = 0.95, auto.scale = TRUE, labels = TRUE) 
{
  data <- as.numeric(data)
  qq <- 0
  if (ci) 
    qq <- qnorm(1 - (1 - ci)/2)
  x <- thresholds
  gpd.dummy <- function(thr, data) {
    out <- gpd(data = data, threshold = thr, information = "expected")
    c(out$n.exceed, out$par.ests[1], out$par.ses[1])
  }
  mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("exceedances", "shape", 
                          "se", "thresholds"), NULL)
  thresh <- mat[4, ]
  exceed <- mat[1, ]
  y <- mat[2, ]
  yrange <- range(y)
  if (ci) {
    u <- y + mat[3, ] * qq
    l <- y - mat[3, ] * qq
    yrange <- range(y, u, l)
  }
  index <- x
  if (reverse) 
    index <- -x
  if (auto.scale) 
    plot(index, y, ylim = yrange, type = "l", xlab = "", 
         ylab = "", axes = FALSE)
  else plot(index, y, type = "l", xlab = "", ylab = "",axes = FALSE)
  abline(v=1.0, lwd=1.5, lty=2, col="red")
  axis(1, at = index, labels = paste(round(thresholds,2)), tick = T)
  axis(2)
  axis(3, at = index, labels = paste(seq(0,0.95,by=0.05)), tick = T)
  box()
  if (ci) {
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if (labels) {
    labely <- expression(paste("",hat(xi)))
    title(xlab = "Threshold", ylab = labely)
    mtext("Quantile", side = 3, line = 3)
  }
  invisible(mat)
}

shapestabboot <- function (data, thresholds,Q, reverse = TRUE,#models = 30, start = 15, end = 500,, 
                       ci = 0.95, auto.scale = TRUE, labels = TRUE, boot=FALSE, m.boot=200) 
{
  data <- as.numeric(data)
  n <- length(data)
  qq <- 0
  if (ci) 
    qq <- qnorm(1 - (1 - ci)/2)
  x <- thresholds #trunc(seq(from = length(data[threshold, to = start, 
  #length = models))
  gpd.dummy <- function(thr, data) {
    out <- gpd(data = data, threshold = thr, information = "expected")
    c(out$n.exceed, out$par.ests[1],out$par.ests[2], out$par.ses[1])
  }
  mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("exceedances", "shape","scale", 
                          "se", "thresholds"), NULL)
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
      shapes_sorted <- sort(shapes)
      u_boot[i] <- shapes_sorted[0.975*m.boot]
      l_boot[i] <- shapes_sorted[0.025*m.boot]
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
  #abline(v=1.0, lwd=1.5, lty=2, col="red")
  axis(1, at = index, labels = paste(round(thresholds,2)), tick = T)
  axis(2)
  axis(3, at = index, labels = paste(Q), tick = T)
  box()
  if (ci) {
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if(boot){
    lines(index, u_boot, lty = 3, col = "blue")
    lines(index, l_boot, lty = 3, col = "blue")
  }
  if (labels) {
    labely <- expression(paste("",hat(xi)))
    # if (ci) 
    #   labely <- paste(labely, " (CI, p = ", ci, ")", 
    #                   sep = "")
    title(xlab = "Threshold", ylab = labely)
    mtext("Quantile", side = 3, line = 3)
  }
  invisible(mat)
}


scalestab <- function (data, thresholds, reverse = TRUE,#models = 30, start = 15, end = 500,, 
                       ci = 0.95, auto.scale = TRUE, labels = TRUE) 
{
  data <- as.numeric(data)
  qq <- 0
  if (ci) 
    qq <- qnorm(1 - (1 - ci)/2)
  x <- thresholds #trunc(seq(from = length(data[threshold, to = start, 
  #length = models))
  gpd.dummy <- function(thr, data) {
    out <- gpd(data = data, threshold = thr, information = "expected")
    mod_scale <- out$par.ests[2] - out$par.ests[1]*thr
    mod_se <- out$par.ses[2] + (thr^2)*out$par.ses[1]
    c(out$n.exceed, mod_scale,mod_se)
  }
  mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("exceedances", "scale", 
                          "se", "thresholds"), NULL)
  thresh <- mat[4, ]
  exceed <- mat[1, ]
  y <- mat[2, ]
  yrange <- range(y)
  if (ci) {
    u <- y + mat[3, ] * qq
    l <- y - mat[3, ] * qq
    yrange <- range(y, u, l)
  }
  index <- x
  if (reverse) 
    index <- -x
  if (auto.scale) 
    plot(index, y, ylim = yrange, type = "l", xlab = "", 
         ylab = "", axes = FALSE)
  else plot(index, y, type = "l", xlab = "", ylab = "",axes = FALSE)
  axis(1, at = index, labels = paste(x), tick = T)
  axis(2)
  axis(3, at = index, labels = paste(format(signif(exceed, 3))), tick = T)
  box()
  if (ci) {
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if (labels) {
    labely <- "Modified Scale"
    if (ci) 
      labely <- paste(labely, " (CI, p = ", ci, ")", 
                      sep = "")
    title(xlab = "Threshold", ylab = labely)
    mtext("Number of Exceedances", side = 3, line = 3)
  }
  invisible(mat)
}


#----------------------------------Plots------------------------------------
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
set.seed(12345)
dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
thresh <- quantile(data, seq(0,0.95,by=0.05), names=FALSE)
shapestab(data, thresholds = thresh, reverse = F)
set.seed(12345)
dat1 <- runif(80, 0.5, 1.0)
dat2 <- rgpd(400, shape=0.1, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
length(data[data < 1])/480
thresh <- quantile(data, seq(0,0.95,by=0.05), names=FALSE)
shapestab(data, thresholds = thresh, reverse = F)
set.seed(12345)
dat1 <- runif(200, 0.5, 1.0)
dat2 <- rgpd(1000, shape=-0.05, scale=0.5, mu=1.0)
data <- c(dat1, dat2)
thresh <- quantile(data, seq(0,0.95,by=0.05), names=FALSE)
shapestab(data, thresholds = thresh, reverse = F)
set.seed(12345)
data_all <- rgpd(4000, shape=0.1, scale=0.5, mu=0)
cens_thr<-1.0*rbeta(length(data_all),1,0.5)
keep <- data_all>cens_thr
data_keep <- data_all[keep]
data <- sample(data_keep, 200, replace = FALSE)
Q <-  seq(0,0.95,by=0.05)
thresh <- quantile(data, seq(0,0.95,by=0.05), names=FALSE)
shapestab(data, thresholds = thresh, reverse=F)
shapestabboot(data, thresholds = thresh, Q=Q, reverse=F, boot=TRUE,m.boot=200)
abline(h=-0.34, col="green")
abline(h=-0.38, col="green")
install.packages("mev")
library(mev)
NC.diag
NC.diag(data, u=thresh,do.LRT = FALSE,size = NULL,plot = TRUE,axes=FALSE,col.lab="white",xi.tol = 0.001 )

axis(1, at = thresh, labels = paste(round(thresh,2)), tick = T)
axis(2)
axis(3, at = thresh, pos=0.92,labels = paste(seq(0,0.95,by=0.05)), tick = T, col.ticks = "white", col="white")

title(xlab = "Threshold", ylab = "p-value")
mtext("Number of exceedances/Quantile", side = 3, line = 3)
abline(v=1.0, col="red", lty=2)
####---------Nidd river data----------------------------
library(evir)
data(nidd.thresh)
Q <- seq(0,0.95,by=0.01)
thresholds <- quantile(nidd.thresh,Q, names=F)
dev.new(width=9.17, height=4,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')
shapestab(nidd.thresh, thresholds = thresholds, reverse = F)
shapestabboot(nidd.thresh, thresholds = thresholds, Q=Q, reverse=F, boot=TRUE,m.boot=200)
length(thresholds)

quantile(nidd.thresh, 0.9)
thr <- 109.08
data <- nidd.thresh[nidd.thresh > thr] - thr
(opt <- optim(GPD_LL, par=c(mean(data), 0.1), z=data, control = list(fnscale=-1), hessian=T)) 
sqrt(diag(solve(-1*opt$hessian)))
