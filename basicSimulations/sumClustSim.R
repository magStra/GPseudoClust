#simulated data set 2 without dropout
library(R.matlab)
library(nonparametricSummaryPSM)
library(mcclust)

true_clusters <- c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))
A <- readMat("PSMSim2NoDropout.mat")
PSMS <- A$PSM
npSummarySim1 <- processPSMs(PSMS)
Am <- apply(PSMS,c(1,2),mean)
m <- maxpear(Am)
GPseudoClust_PY <- npSummarySim1$sumClustPEAR
GPseudoClust_DP <- npSummarySim1$sumClustPEAR_DP
GPseudoClust_mean <- m$cl
save(file="npSumClustSim2.rda", GPseudoClust_PY,GPseudoClust_DP,GPseudoClust_mean)


true_clusters <- c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))
A <- readMat("PSMSimNoDropout.mat")
PSMS <- A$PSM
npSummarySim1 <- processPSMs(PSMS)
Am <- apply(PSMS,c(1,2),mean)
Am <- (Am+t(Am))/2
Am[Am>1]<- 1
Am[Am<0]<-0
Am <- (Am+t(Am))/2
m <- maxpear(Am)
GPseudoClust_PY <- npSummarySim1$sumClustPEAR
GPseudoClust_DP <- npSummarySim1$sumClustPEAR_DP
GPseudoClust_mean <- m$cl
save(file="npSumClustSim.rda", GPseudoClust_PY,GPseudoClust_DP,GPseudoClust_mean)
