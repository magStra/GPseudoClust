library(nonparametricSummaryPSM)
library(R.matlab)
library(PReMiuM)
library(mcclust)
library(ClusterR)
summaryResults <- list()
ARI_PY <- rep(0,24)
NMI_PY <- rep(0,24)
FMI_PY <- rep(0,24)
ARI_DPM <- rep(0,24)
NMI_DPM <- rep(0,24)
FMI_DPM <- rep(0,24)

for (k in 1:24){
  print(k)
  simData <- readMat(sprintf("simLinear%d.mat",k))
  PSMs = array(0,dim=c(simData$nGenes,simData$nGenes,24))
  for (j in 1:24)
  {
    A <- as.matrix(read.table(sprintf("simLinear%d_Results_PSM_Chain%d.csv",k,j),sep=","))
    PSMs[,,j] <-  A
  }
  summaryResults[[k]] <- processPSMs(PSMs)
  trueClust <- as.vector(simData$z)
  ARI_PY[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR,
                                  method = "adjusted_rand_index")
  ARI_DPM[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR_DP,
                                   method = "adjusted_rand_index")
  FMI_PY[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR,
                                  method = "fowlkes_mallows_index")
  FMI_DPM[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR_DP,
                                   method = "fowlkes_mallows_index")
  NMI_PY[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR,
                                  method = "nmi")
  NMI_DPM[k] = external_validation(trueClust,summaryResults[[k]]$sumClustPEAR_DP,
                                   method = "nmi")
  
}
save(summaryResults,ARI_PY,ARI_DPM,FMI_PY,FMI_DPM,NMI_PY,NMI_DPM,file="LinearNonParametricSummaryPSM.rda")


