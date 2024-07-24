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

#Table 2: RMSE of threshold choices for Cases 1-4
(rmse_thr_cases1234 <- data.frame(EQD=c(rmse(eqd_case1$thr,1), rmse(eqd_case2$thr,1), rmse(eqd_case3$thr,1), rmse(eqd_case4$thr,1)), Wadsworth=c(rmse(wads_case1$thr[!is.na(wads_case1$thr)],1), rmse(wads_case2$thr[!is.na(wads_case2$thr)],1), rmse(wads_case3$thr[!is.na(wads_case3$thr)],1), rmse(wads_case4$thr[!is.na(wads_case4$thr)],1)), Northrop=c(rmse(north_case1$thr,1), rmse(north_case2$thr,1), rmse(north_case3$thr,1), rmse(north_case4$thr,1))))

write.csv(rmse_thr_cases1234, "output/tables/Table_2_rmse_threshold_choice_case1-4.csv")



