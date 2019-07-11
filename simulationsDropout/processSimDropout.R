library(ClusterR)
library(R.matlab)



sumClusts_lkkm <- as.matrix(read.table("sumClust_simDropout_lmkk.csv",header=F,sep=","))
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))


ARI_lkkm = c()
NMI_lkkm = c()
FMI_lkkm = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "adjusted_rand_index")
  ARI_lkkm = c(ARI_lkkm,b)
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "fowlkes_mallows_index")
  FMI_lkkm = c(FMI_lkkm,b)
  b =external_validation(true_clusters, sumClusts_lkkm[,j],
                         method = "nmi")
  NMI_lkkm = c(NMI_lkkm,b)
}


write.table(rbind(ARI_lkkm,FMI_lkkm,NMI_lkkm),"scores_dropout_lmkk.csv",row.names=F,col.names=F,sep=",")

# consensus clustering
sumClusts_dropout <- matrix(nrow=52,ncol=100)
library(mcclust)
library(ClusterR)
for (j in 1:100)
{PSM_dropout = as.matrix(read.table(sprintf("PSM_simDropout_consensus_%d.csv",j) ,sep=",",header=F))
sumClusts_dropout[,j] = norm.label(maxpear(PSM_dropout,method="comp",max.k=10)$cl)
}



ARI_Sub = c()
NMI_Sub = c()
FMI_Sub = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "adjusted_rand_index")
  ARI_Sub = c(ARI_Sub,b)
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "fowlkes_mallows_index")
  FMI_Sub = c(FMI_Sub,b)
  b =external_validation(true_clusters, sumClusts_dropout[,j],
                         method = "nmi")
  NMI_Sub = c(NMI_Sub,b)
}


write.table(rbind(ARI_Sub,FMI_Sub,NMI_Sub),"scores_dropout_consensus.csv",row.names=F,col.names=F,sep=",")


##higher levels of dropout

sumClusts_lkkm <- as.matrix(read.table("sumClust_simDropout30_lmkk.csv",header=F,sep=","))
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))


ARI_lkkm = c()
NMI_lkkm = c()
FMI_lkkm = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "adjusted_rand_index")
  ARI_lkkm = c(ARI_lkkm,b)
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "fowlkes_mallows_index")
  FMI_lkkm = c(FMI_lkkm,b)
  b =external_validation(true_clusters, sumClusts_lkkm[,j],
                         method = "nmi")
  NMI_lkkm = c(NMI_lkkm,b)
}


write.table(rbind(ARI_lkkm,FMI_lkkm,NMI_lkkm),"scores_dropout30_lmkk.csv",row.names=F,col.names=F,sep=",")

#lmkk, more dropout, assume number of clusters known
sumClusts_lkkm <- as.matrix(read.table("sumClust_simDropout30_lmkkNumberKnown.csv",header=F,sep=","))
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))


ARI_lkkm = c()
NMI_lkkm = c()
FMI_lkkm = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "adjusted_rand_index")
  ARI_lkkm = c(ARI_lkkm,b)
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "fowlkes_mallows_index")
  FMI_lkkm = c(FMI_lkkm,b)
  b =external_validation(true_clusters, sumClusts_lkkm[,j],
                         method = "nmi")
  NMI_lkkm = c(NMI_lkkm,b)
}


write.table(rbind(ARI_lkkm,FMI_lkkm,NMI_lkkm),"scores_dropout30_lmkkNumberKnown.csv",row.names=F,col.names=F,sep=",")


#consensus clustering
sumClusts_dropout <- matrix(nrow=52,ncol=100)
library(mcclust)
library(ClusterR)
for (j in 1:100)
{PSM_dropout = as.matrix(read.table(sprintf("PSM_simDropout_consensus30_%d.csv",j) ,sep=",",header=F))
sumClusts_dropout[,j] = norm.label(maxpear(PSM_dropout,method="comp",max.k=10)$cl)
}



ARI_Sub = c()
NMI_Sub = c()
FMI_Sub = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "adjusted_rand_index")
  ARI_Sub = c(ARI_Sub,b)
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "fowlkes_mallows_index")
  FMI_Sub = c(FMI_Sub,b)
  b =external_validation(true_clusters, sumClusts_dropout[,j],
                         method = "nmi")
  NMI_Sub = c(NMI_Sub,b)
}


write.table(rbind(ARI_Sub,FMI_Sub,NMI_Sub),"scores_dropout30_consensus.csv",row.names=F,col.names=F,sep=",")




sumClusts_lkkm <- as.matrix(read.table("sumClust_simDropoutMixed_lmkk.csv",header=F,sep=","))
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))


ARI_lkkm = c()
NMI_lkkm = c()
FMI_lkkm = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "adjusted_rand_index")
  ARI_lkkm = c(ARI_lkkm,b)
  b = external_validation(true_clusters, sumClusts_lkkm[,j],
                          method = "fowlkes_mallows_index")
  FMI_lkkm = c(FMI_lkkm,b)
  b =external_validation(true_clusters, sumClusts_lkkm[,j],
                         method = "nmi")
  NMI_lkkm = c(NMI_lkkm,b)
}


write.table(rbind(ARI_lkkm,FMI_lkkm,NMI_lkkm),"scores_dropoutMixed_lmkk.csv",row.names=F,col.names=F,sep=",")

#consensus clustering
sumClusts_dropout <- matrix(nrow=52,ncol=100)
library(mcclust)
library(ClusterR)
for (j in 1:100)
{PSM_dropout = as.matrix(read.table(sprintf("PSM_simDropoutMixed_consensus_%d.csv",j) ,sep=",",header=F))
sumClusts_dropout[,j] = norm.label(maxpear(PSM_dropout,method="comp",max.k=10)$cl)
}



ARI_Sub = c()
NMI_Sub = c()
FMI_Sub = c()

for (j in 1:100){
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "adjusted_rand_index")
  ARI_Sub = c(ARI_Sub,b)
  b = external_validation(true_clusters, sumClusts_dropout[,j],
                          method = "fowlkes_mallows_index")
  FMI_Sub = c(FMI_Sub,b)
  b =external_validation(true_clusters, sumClusts_dropout[,j],
                         method = "nmi")
  NMI_Sub = c(NMI_Sub,b)
}


write.table(rbind(ARI_Sub,FMI_Sub,NMI_Sub),"scores_dropoutMixed_consensus.csv",row.names=F,col.names=F,sep=",")

#####
#use DPM and PY
library(nonparametricSummaryPSM)
#less dropout
for (j in 1:100)
{
  print(j)
  PSMs <- readMat(sprintf("PSM_simDropout_Subsample_%d.mat",j))
  PSMS <- PSMs$PSMs
  #some entries are a small trace larger than one by numeral error
  PSMS[PSMS>1] = 1
  for (k in 1:24)
  {
    PSMS[,,k] = 0.5*( PSMS[,,k]+t(PSMS[,,k]))#after correcting above for entries
    #larger than one, make sure the matrix is perfectly symmetric
  }
  A <- processPSMs(PSMS)
 save(A,file=sprintf("PSM_simDropout_DPM_%d.mat",j))
}

#more dropout
for (j in 1:100)
{
  print(j)
  PSMs <- readMat(sprintf("PSM_simDropout30_Subsample_%d.mat",j))
  PSMS <- PSMs$PSMs
  #some entries are a small trace larger than one by numeral error
  PSMS[PSMS>1] = 1
  for (k in 1:24)
  {
    PSMS[,,k] = 0.5*( PSMS[,,k]+t(PSMS[,,k]))#after correcting above for entries
    #larger than one, make sure the matrix is perfectly symmetric
  }
  A <- processPSMs(PSMS)
  save(A,file=sprintf("PSM_simDropout_DPM_%d.mat",j))
}

#mixed dropout
for (j in 1:100)
{
  print(j)
  PSMs <- readMat(sprintf("dropout/PSM_simDropout_Mixed_%d.mat",j))
  PSMS <- PSMs$PSMs
  #some entries are a small trace larger than one by numeral error
  PSMS[PSMS>1] = 1
  for (k in 1:24)
  {
    PSMS[,,k] = 0.5*( PSMS[,,k]+t(PSMS[,,k]))#after correcting above for entries
    #larger than one, make sure the matrix is perfectly symmetric
  }
  A <- processPSMs(PSMS)
  save(A,file=sprintf("PSM_simDropout_DPM_%d.mat",j))
}



#for PY and DPM methods, compute ARI with true clustering
true_clusters = c(rep(1,8),rep(2,4),rep(3,12),rep(4,16),rep(5,12))
ARI_PY = c()
FMI_PY = c()
NMI_PY = c()
ARI_DPM = c()
FMI_DPM = c()
NMI_DPM = c()
for (j in 1:100)
{
  print(j)
  clusObj = readMat(sprintf("PSM_simDropout_DPM_%d.mat",j))

  ARI_PY <- c(ARI_PY, external_validation(true_clusters,clusObj$sumClustPEAR,
                                                      method = "adjusted_rand_index"))
  ARI_DPM <- c(ARI_DPM, external_validation(true_clusters,clusObj$sumClustPEAR.DP,
                                          method = "adjusted_rand_index"))
  b = external_validation(true_clusters, clusObj$sumClustPEAR,
                          method = "fowlkes_mallows_index")
  FMI_PY = c(FMI_PY,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR,
                         method = "nmi")
  NMI_PY = c(NMI_PY,b)
  b = external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                          method = "fowlkes_mallows_index")
  FMI_DPM = c(FMI_DPM,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                         method = "nmi")
  NMI_DPM = c(NMI_DPM,b)
}
write.table(rbind(ARI_PY,FMI_PY,NMI_PY),"scores_dropout_PitmanYor.csv",row.names=F,col.names=F,sep=",")
write.table(rbind(ARI_DPM,FMI_DPM,NMI_DPM),"scores_dropout_DPM.csv",row.names=F,col.names=F,sep=",")
#for more dropout
ARI_PY = c()
FMI_PY = c()
NMI_PY = c()
ARI_DPM = c()
FMI_DPM = c()
NMI_DPM = c()
for (j in 1:100)
{
  print(j)
  clusObj = readMat(sprintf("PSM_simDropout30_DPM_%d.mat",j))

  ARI_PY <- c(ARI_PY, external_validation(true_clusters,clusObj$sumClustPEAR,
                                            method = "adjusted_rand_index"))
  b = external_validation(true_clusters, clusObj$sumClustPEAR,
                          method = "fowlkes_mallows_index")
  FMI_PY = c(FMI_PY,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR,
                         method = "nmi")
  NMI_PY = c(NMI_PY,b)


  ARI_DPM <- c(ARI_DPM, external_validation(true_clusters,clusObj$sumClustPEAR.DP,
                                          method = "adjusted_rand_index"))
  b = external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                          method = "fowlkes_mallows_index")
  FMI_DPM = c(FMI_DPM,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                         method = "nmi")
  NMI_DPM = c(NMI_DPM,b)
}
write.table(rbind(ARI_PY,FMI_PY,NMI_PY),"scores_dropout30_PitmanYor.csv",row.names=F,col.names=F,sep=",")
write.table(rbind(ARI_DPM,FMI_DPM,NMI_DPM),"scores_dropout30_DPM.csv",row.names=F,col.names=F,sep=",")

#for mixed dropout
ARI_PY = c()
FMI_PY = c()
NMI_PY = c()
ARI_DPM = c()
FMI_DPM = c()
NMI_DPM = c()
for (j in 1:100)
{
  print(j)
  clusObj = readMat(sprintf("PSM_simDropout_Mixed_DPM_%d.mat",j))

  ARI_PY <- c(ARI_PY, external_validation(true_clusters,clusObj$sumClustPEAR,
                                            method = "adjusted_rand_index"))
  b = external_validation(true_clusters, clusObj$sumClustPEAR,
                          method = "fowlkes_mallows_index")
  FMI_PY = c(FMI_PY,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR,
                         method = "nmi")
  NMI_PY = c(NMI_PY,b)

  ARI_DPM <- c(ARI_DPM, external_validation(true_clusters,clusObj$sumClustPEAR.DP,
                                          method = "adjusted_rand_index"))
  b = external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                          method = "fowlkes_mallows_index")
  FMI_DPM = c(FMI_DPM,b)
  b =external_validation(true_clusters, clusObj$sumClustPEAR.DP,
                         method = "nmi")
  NMI_DPM = c(NMI_DPM,b)
}
write.table(rbind(ARI_PY,FMI_PY,NMI_PY),"scores_dropoutMixed_PitmanYor.csv",row.names=F,col.names=F,sep=",")
write.table(rbind(ARI_DPM,FMI_DPM,NMI_DPM),"scores_dropoutMixed_DPM.csv",row.names=F,col.names=F,sep=",")
