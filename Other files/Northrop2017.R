library(threshr)

X  <- rgpd(1000, shape=0.1, scale=0.5)
u <- quantile(X, seq(0.5,0.9, by=0.05), names=FALSE)

ith1 <- ithresh(X, u)

par(mfrow=c(1,1))
plot(ith1)
