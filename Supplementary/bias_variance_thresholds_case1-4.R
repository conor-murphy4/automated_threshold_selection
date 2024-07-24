library(Metrics)

#Threshold choices
eqd_case1 <- read.csv("output/threshold_selection/eqd_case1.csv", row.names = 1, header = T)
eqd_case2 <- read.csv("output/threshold_selection/eqd_case2.csv", row.names = 1, header = T)
eqd_case3 <- read.csv("output/threshold_selection/eqd_case3.csv", row.names = 1, header = T)
eqd_case4 <- read.csv("output/threshold_selection/eqd_case4.csv", row.names = 1, header = T)

wads_case1 <- read.csv("output/threshold_selection/wads_case1.csv", row.names = 1, header = T)
wads_case2 <- read.csv("output/threshold_selection/wads_case2.csv", row.names = 1, header = T)
wads_case3 <- read.csv("output/threshold_selection/wads_case3.csv", row.names = 1, header = T)
wads_case4 <- read.csv("output/threshold_selection/wads_case4.csv", row.names = 1, header = T)

north_case1 <- read.csv("output/threshold_selection/north_case1.csv", row.names = 1, header = T)
north_case2 <- read.csv("output/threshold_selection/north_case2.csv", row.names = 1, header = T)
north_case3 <- read.csv("output/threshold_selection/north_case3.csv", row.names = 1, header = T)
north_case4 <- read.csv("output/threshold_selection/north_case4.csv", row.names = 1, header = T)


#Table S.7: Bias and variance of threshold choice for Cases 1-4
(bias_thr_cases1234 <- data.frame(EQD=c(bias(eqd_case1$thr,1), bias(eqd_case2$thr,1), bias(eqd_case3$thr,1), bias(eqd_case4$thr,1)), Wadsworth=c(bias(wads_case1$thr[!is.na(wads_case1$thr)],1), bias(wads_case2$thr[!is.na(wads_case2$thr)],1), bias(wads_case3$thr[!is.na(wads_case3$thr)],1), bias(wads_case4$thr[!is.na(wads_case4$thr)],1)), Northrop=c(bias(north_case1$thr,1), bias(north_case2$thr,1), bias(north_case3$thr,1), bias(north_case4$thr,1))))
write.csv(bias_thr_cases1234, "output/tables/Table_S.7_bias_threshold_choice_case1-4.csv")

(variance_thr_cases1234 <- data.frame(EQD=c(var(eqd_case1$thr), var(eqd_case2$thr), var(eqd_case3$thr), var(eqd_case4$thr)), Wadsworth=c(var(wads_case1$thr[!is.na(wads_case1$thr)]), var(wads_case2$thr[!is.na(wads_case2$thr)]), var(wads_case3$thr[!is.na(wads_case3$thr)]), var(wads_case4$thr[!is.na(wads_case4$thr)])), Northrop=c(var(north_case1$thr), var(north_case2$thr), var(north_case3$thr), var(north_case4$thr))))
write.csv(variance_thr_cases1234, "output/tables/Table_S.7_variance_threshold_choice_case1-4.csv")
