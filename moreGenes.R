library(nonparametricSummaryPSM)
library(R.matlab)
library(ClusterR)

#simulation with 5000 genes
PSMs <- array(dim=c(5000,5000,12))
for (j in 1:10){
  print(j)
  A <- readMat(sprintf("PSMsMoreGenes%d.mat",j))
    PSMs[,,j] <- A$PSM
}

summaryResults <- processPSMs(PSMs)
save(summaryResults,file="Clust5000Genes.rda")

#compare to true clustering
trueClust <- as.vector(readMat("simMoreGenes111.mat")$z)
ARI_PY = external_validation(trueClust,summaryResults$sumClustPEAR,
                                   method = "adjusted_rand_index")
ARI_DPM = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                    method = "adjusted_rand_index")
FMI_PY = external_validation(trueClust,summaryResults$sumClustPEAR,
                                   method = "fowlkes_mallows_index")
FMI_DPM = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                    method = "fowlkes_mallows_index")
NMI_PY = external_validation(trueClust,summaryResults$sumClustPEAR,
                                   method = "nmi")
NMI_DPM = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                    method = "nmi")

