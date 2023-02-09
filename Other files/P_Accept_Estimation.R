source("_gpd.R")

# Distribution of P_accept estimation

n <- 100000
u <- 1.0
P_Accept <- numeric(n)
for(ii in 1:n){
  saveRDS(ii, "STORM_output/iter.rds")
  data_all <- rgpd(1000, shape=0.1, scale=0.5, mu=0)
  cens_thr<-u*rbeta(length(data_all),1,0.5)
  keep <- data_all>cens_thr
  data <- data_all[keep]
  p_accept <- length(data[data<1])/length(data)
  P_Accept[ii] <- p_accept 
}


#True P_accept estimation
m1 <- 10000
u <- 1.0
data_all <- rgpd(m1, shape=0.1, scale=0.5, mu=0)
cens_thr<-u*rbeta(length(data_all),1,0.5)
keep <- data_all>cens_thr
data <- data_all[keep]
p_accept_true <- length(data[data<1])/length(data)

(error <- 2*(p_accept_true)*(1- p_accept_true)/length(data))
p_accept_true


