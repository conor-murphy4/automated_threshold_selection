library(evir)
source("src/helper_functions.R")

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
    title(xlab = "Threshold", ylab = labely, cex.lab=1)
    mtext("Quantile", side = 3, line = 2.5, cex=1)
  }
  invisible(mat)
}



