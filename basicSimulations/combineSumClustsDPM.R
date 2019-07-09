#This file contains R functions to compute summary PSMs from multiple PSMs, obtained eg. by subsampling.

library(PReMiuM)
library(mcclust)
library(R.matlab)
library(ClusterR)
#first compute a summary clustering from each PSM
computeSumClustPEAR = function(PSM,maxCl=10)
{
  return(norm.label(maxpear(PSM,method="comp",max.k=maxCl)$cl))
}

#now compute weights for the samples and overall summary clustering
computeWeightsSumClust = function(allocs)
{
  colN <- c()
  for (j in 1:dim(allocs)[2])
  {
    colN <- c(colN,sprintf("var%d",j))
  }
  colnames(allocs) <- colN
  allocs <- as.data.frame(allocs)
  output <- c()
  clust <- profRegr(covNames=colN, data=allocs,nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, 
           nClusInit=20,  xModel="Discrete",  sampler="Truncated",  excludeY=T, 
           varSelectType="BinaryCluster",dPitmanYor =0.8,alpha=5)

   output$clusObj <- clust
   ng <- dim(allocs)[1]
   dissimObj <- calcDissimilarityMatrix(clust)
   PSM <- matrix(0,ng,ng)
   PSM[lower.tri(PSM, diag=F)] <- 1-dissimObj$disSimMat
   diag(PSM) <- 0.5
   PSM <- t(PSM) + PSM
   output$PSM <- PSM
   clusObj <- calcOptimalClustering(dissimObj)
   output$sumCl <- clusObj
   #get the weights from the output file
   output$rho <- read.table("output_rho.txt",header=F)
   output
}

computeWeightsSumClustDPM = function(allocs)
{
  colN <- c()
  for (j in 1:dim(allocs)[2])
  {
    colN <- c(colN,sprintf("var%d",j))
  }
  colnames(allocs) <- colN
  allocs <- as.data.frame(allocs)
  output <- c()
  clust <- profRegr(covNames=colN, data=allocs,nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, 
                    nClusInit=20,  xModel="Discrete",   excludeY=T, 
                    varSelectType="BinaryCluster")
  
  output$clusObj <- clust
  ng <- dim(allocs)[1]
  dissimObj <- calcDissimilarityMatrix(clust)
  PSM <- matrix(0,ng,ng)
  PSM[lower.tri(PSM, diag=F)] <- 1-dissimObj$disSimMat
  diag(PSM) <- 0.5
  PSM <- t(PSM) + PSM
  output$PSM <- PSM
  clusObj <- calcOptimalClustering(dissimObj)
  output$sumCl <- clusObj
  #get the weights from the output file
  output$rho <- read.table("output_rho.txt",header=F)
  output
}




processPSMs <- function(PSMs,fileName)
{
  PSMs[PSMs>1]=1
  PSMs[PSMs<0]=0
  d = dim(PSMs)
  for (k in 1:d[3])
  {diag(PSMs[,,k]) <- 1
  PSMs[,,k] = 0.5*(PSMs[,,k]+t(PSMs[,,k]))}
  
  sumClusts <- apply(PSMs,3,computeSumClustPEAR)
  clustOb <- computeWeightsSumClust(sumClusts)
  clustObDPM <- computeWeightsSumClustDPM(sumClusts)
  meanWeights <- apply(clustOb$rho,2,mean)
  meanWeightsDPM <- apply(clustObDPM$rho,2,mean)
  #summary PSM obtained from the PY process clustering of the summary clusterings
  
  weights <-meanWeights/sum(meanWeights)
  weightsDPM <-meanWeightsDPM/sum(meanWeightsDPM)
  PSM_PY <- clustOb$PSM
  PSM_DP <- clustObDPM$PSM
  sumClustPSM <- clustOb$sumCl$clustering
  sumClustPSMDP <- clustObDPM$sumCl$clustering
  weightedPSM <- matrix(0,d[1],d[1])
  for (j in 1:d[3])
  {
    weightedPSM <- weightedPSM + PSMs[,,j]*weights[j]
  }
  weightedPSM_DP <- matrix(0,d[1],d[1])
  for (j in 1:d[3])
  {
    weightedPSM_DP<- weightedPSM_DP + PSMs[,,j]*weightsDPM[j]
  }
  weightedPSM[weightedPSM>1] <- 1
  weightedPSM[weightedPSM<0] <- 0
  diag(weightedPSM) <- 1
  weightedPSM = 0.5*(weightedPSM+t(weightedPSM))
  weightedPSM_DP[weightedPSM_DP>1] <- 1
  weightedPSM_DP[weightedPSM_DP<0] <- 0
  diag(weightedPSM_DP) <- 1
  weightedPSM_DP = 0.5*(weightedPSM_DP+t(weightedPSM_DP))
  #summary clustering using PEAR from the weighted PSM
  sumClustPEAR <- computeSumClustPEAR(weightedPSM)
  sumClustPEAR_DP <- computeSumClustPEAR(weightedPSM_DP)
  writeMat(fileName,weightedPSM=weightedPSM,sumClustPEAR=sumClustPEAR,
           sumClustPEAR_DP = sumClustPEAR_DP,sumClustDPM=sumClustPSMDP, sumClustPY = sumClustPSM,
           weights=weights, weightedPSM_DP=weightedPSM_DP,weights_DP=weightsDPM)
  
}

#example: Stumpf data
PSMs_Stumpf <- array(dim=c(94,94,192))
PSMs_StumpfE <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_StumpfE.mat")
PSMs_StumpfR <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_StumpfR.mat")
PSMs_Stumpf[,,1:96] <- PSMs_StumpfE$PSMs
PSMs_Stumpf[,,97:192] <- PSMs_StumpfR$PSMs
processPSMs(PSMs_Stumpf,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/StumpfPY.mat")



####example Shalek data
library(R.matlab)
PSMs_Shalek <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMsShalek.mat")
PSMs_Shalek <- PSMs_Shalek$PSMsShalek
processPSMs(PSMs_Shalek,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/ShalekPY.mat")
#Shalek data, first 24 chains
PSMs_Shalek24 <- PSMs_Shalek[,,1:24]
processPSMs(PSMs_Shalek24,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/ShalekPY24.mat")
#Shalek data, first 4 chains
PSMs_Shalek4 <- PSMs_Shalek[,,1:4]
processPSMs(PSMs_Shalek4,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/ShalekPY4.mat")



# #example 3: Goettgens data
# branch 2
PSMs_ESC2 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC_2.mat")
PSMs_ESC2 <- PSMs_ESC2$PSMs.ESC.1
processPSMs(PSMs_ESC2,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC2_PY.mat")
#branch 1
PSMs_ESC1 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC_1.mat")
PSMs_ESC1 <- PSMs_ESC1$PSMs.ESC.1
processPSMs(PSMs_ESC1,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC1_PY.mat")
#branch 3
PSMs_ESC3 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC_3.mat")
PSMs_ESC3 <- PSMs_ESC3$PSMs.ESC.3
processPSMs(PSMs_ESC3,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_ESC3_PY.mat")


#Sasagawa
PSMs_mes1_12 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_mes1_12.mat")
PSMs_mes13_24 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_mes13_24.mat")
PSMs_mes25_36 <- readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_mes25_36.mat")
PSMs_mes <- array(0,dim=c(600,600,36))
PSMs_mes[,,1:12] <- PSMs_mes1_12$PSMS
PSMs_mes[,,13:24] <- PSMs_mes13_24$PSMS
PSMs_mes[,,25:36] <- PSMs_mes25_36$PSMS
processPSMs(PSMs_mes,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/mesPY.mat")
# Klein
PSMs_Klein1 = array(data=0,dim=c(2047,2047,28))
for (j in 1:28)
{
  print(j)
  A <- unlist(read.table(sprintf("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_Klein_1_%d.csv",j),
                  header=F,sep = ","))
  PSMs_Klein1[,,j] <- A
}
processPSMs(PSMs_Klein1,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Klein1PY.mat")

PSMs_Klein21 = array(data=0,dim=c(2047,2047,28))
for (j in 1:28)
{
  print(j)
  A <- unlist(read.table(sprintf("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_Klein_2_1_%d.csv",j),
                         header=F,sep = ","))
  PSMs_Klein21[,,j] <- A
}
processPSMs(PSMs_Klein21,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Klein21PY.mat")

PSMs_Klein22 = array(data=0,dim=c(2047,2047,28))
for (j in 1:28)
{
  print(j)
  A <- unlist(read.table(sprintf("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_Klein_2_2_%d.csv",j),
                         header=F,sep = ","))
  PSMs_Klein22[,,j] <- A
}
processPSMs(PSMs_Klein22,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Klein22PY.mat")

PSMs_Klein23 = array(data=0,dim=c(2047,2047,28))
for (j in 1:28)
{
  print(j)
  A <- unlist(read.table(sprintf("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_Klein_2_3_%d.csv",j),
                         header=F,sep = ","))
  PSMs_Klein23[,,j] <- A
}
processPSMs(PSMs_Klein23,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Klein23PY.mat")

PSMs_Klein3 = array(data=0,dim=c(2047,2047,28))
for (j in 1:28)
{
  print(j)
  A <- unlist(read.table(sprintf("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/PSMs_Klein_3_%d.csv",j),
                         header=F,sep = ","))
  PSMs_Klein3[,,j] <- A
}
processPSMs(PSMs_Klein3,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Klein3PY.mat")

#simulated data set 2 without dropout
A = readMat("~/Documents/GPseudoClustCode/PSMSim2NoDropout.mat")
PSMS = A$PSM
processPSMs(PSMS,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Sim2NoDP.mat")
Am = apply(PSMS,c(1,2),mean)
m = maxpear(Am)
adjustedRandIndex(m$cl,true_clusters)
external_validation(true_clusters, m$cl,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, m$cl,
                    method = "nmi")
A = readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/Sim2NoDP.mat")
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))
m = maxpear(A$weightedPSM)
adjustedRandIndex(m$cl,true_clusters)
external_validation(true_clusters, m$cl,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, m$cl,
                    method = "nmi")
m = maxpear(A$weightedPSM.DP)
adjustedRandIndex(m$cl,true_clusters)
#lmkk
A = read.table("~/Documents/GPseudoClustCode/sim2lmkk.csv",header=F)
adjustedRandIndex(A$V1,true_clusters)

external_validation(true_clusters, A$V1,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, A$V1,
                       method = "nmi")

#simulated data set 1 without dropout
A = readMat("~/Documents/GPseudoClustCode/PSMSimNoDropout.mat")
PSMS = A$PSM
processPSMs(PSMS,"~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/SimNoDP.mat")
Am = apply(PSMS,c(1,2),mean)
Am = Am/max(Am)
diag(Am) <- 1
m = maxpear(Am)
adjustedRandIndex(m$cl,true_clusters)
external_validation(true_clusters, m$cl,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, m$cl,
                    method = "nmi")
A = readMat("~/Documents/GPseudoClustCode/GPseudoClust/postProcPSM/SimNoDP.mat")
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))
m = maxpear(A$weightedPSM)
adjustedRandIndex(m$cl,true_clusters)
external_validation(true_clusters, m$cl,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, m$cl,
                    method = "nmi")
m = maxpear(A$weightedPSM.DP)
adjustedRandIndex(m$cl,true_clusters)
#lmkk
A = read.table("~/Documents/GPseudoClustCode/simlmkk.csv",header=F)
adjustedRandIndex(A$V1,true_clusters)

external_validation(true_clusters, A$V1,
                    method = "fowlkes_mallows_index")
external_validation(true_clusters, A$V1,
                    method = "nmi")
