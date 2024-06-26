
source('src/eqd_boot.R')

data_matrix <- read.csv("data/Case_1.csv", row.names = 1, header=T) 
data_matrix_gauss <- read.csv("data/Gaussian.csv", row.names = 1, header = T)

data <- data_matrix[,1]
data_gauss <- data_matrix_gauss[,1]

thresholds <- quantile(data,seq(0,0.95, by=0.05))
thresholds_gauss <- quantile(data_gauss, seq(0.5,0.95, by=0.05))

k <- 100
set.seed(12345)
eqd_output_500 <- eqd_boot(data,thresholds, k=k)
eqd_output_500_gauss <- eqd_boot(data_gauss,thresholds_gauss, k=k)

set.seed(12345)
eqd_output_nu <- eqd_mnu_boot(data, thresholds, k=k)
eqd_output_nu_gauss <- eqd_mnu_boot(data_gauss, thresholds_gauss, k=k)

#---------------Plot 1: Boxplots of d_b(u) values for m=500 vs m=nu--------------------
min_boot_dist <- min(c(eqd_output_500$boot_distances, eqd_output_nu$boot_distances))
max_boot_dist <- max(c(eqd_output_500$boot_distances, eqd_output_nu$boot_distances))
min_boot_dist_g <- min(c(eqd_output_500_gauss$boot_distances, eqd_output_nu_gauss$boot_distances))
max_boot_dist_g <- max(c(eqd_output_500_gauss$boot_distances, eqd_output_nu_gauss$boot_distances))

dev.new(width=9.17, height=6,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

boxplot(eqd_output_500$boot_distances~round(thresholds,3), col='red', ylim=c(min_boot_dist, max_boot_dist), xlab="Thresholds", ylab=expression(d[b](u)), main="m=500")
points(eqd_output_500$dists, pch=19)
boxplot(eqd_output_nu$boot_distances~round(thresholds,3), col='blue', xlab="Thresholds", ylim=c(min_boot_dist, max_boot_dist), ylab=expression(d[b](u)), main="m=nu")
points(eqd_output_nu$dists, pch=19)

boxplot(eqd_output_500_gauss$boot_distances~round(thresholds_gauss,3), col='red', ylim=c(min_boot_dist_g, max_boot_dist_g), xlab="Thresholds", ylab=expression(d[b](u)), main="m=500")
points(eqd_output_500_gauss$dists, pch=19)
boxplot(eqd_output_nu_gauss$boot_distances~round(thresholds_gauss,3), col='blue', xlab="Thresholds", ylim=c(min_boot_dist_g, max_boot_dist_g), ylab=expression(d[b](u)), main="m=nu")
points(eqd_output_nu_gauss$dists, pch=19)
#-----------------Plot 2: Difference relative to EQD value-------------------------
dbu_diff<- (eqd_output_500$boot_distances-eqd_output_nu$boot_distances)
rel_diff_to_eqd <- matrix(NA, nrow = length(thresholds), ncol=k)
for(i in 1:length(thresholds)){
  rel_diff_to_eqd[i,] <- dbu_diff[i,]/eqd_output_500$dists[i] 
}

dbu_diff<- (eqd_output_500_gauss$boot_distances-eqd_output_nu_gauss$boot_distances)
rel_diff_to_eqd <- matrix(NA, nrow = length(thresholds_gauss), ncol=k)
for(i in 1:length(thresholds_gauss)){
  rel_diff_to_eqd[i,] <- dbu_diff[i,]/eqd_output_500_gauss$dists[i] 
}

dev.new(width=9.17, height=6,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
boxplot(rel_diff_to_eqd~round(thresholds_gauss,3), col='red', xlab="Thresholds", ylab="Difference in db(u) relative to dE(u)")
points(apply(rel_diff_to_eqd, 1, mean), pch=19)
abline(h=0)


