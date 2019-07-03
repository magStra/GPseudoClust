library(nonparametricSummaryPSM)
library(R.matlab)
library(ClusterR)
library(xtable)

ARI_PY <- matrix(0,3,10)
NMI_PY <- matrix(0,3,10)
FMI_PY <- matrix(0,3,10)
ARI_DPM <- matrix(0,3,10)
NMI_DPM <- matrix(0,3,10)
FMI_DPM <- matrix(0,3,10)
PSMs_PY_10 <- list()
PSMs_PY_20 <- list()
PSMs_PY_30 <- list()
#computing PSMs for simulation with large number of cells
#with 10, 20, and finally 30 subsampled cells per chain and capture time
for (j in 1:10){
  A <- readMat(sprintf("simLarge%d.mat",j))
  trueClust <- A$z[1,]
  for (kk in 1:3){
    k = kk*10
    PSMs <- readMat(sprintf("PSMsLarge%d_%d.mat",k,j))$PSMs

    summaryResults <- processPSMs(PSMs)
    ARI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                    method = "adjusted_rand_index")
    ARI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                     method = "adjusted_rand_index")
    FMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                    method = "fowlkes_mallows_index")
    FMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                     method = "fowlkes_mallows_index")
    NMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                    method = "nmi")
    NMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                     method = "nmi")
    #also save the PSM
    if (kk==1)
      PSMs_PY_10[[j]] <- summaryResults$weightedPSM
    else if (kk==2)
      PSMs_PY_20[[j]] <- summaryResults$weightedPSM
    else
      PSMs_PY_30[[j]] <- summaryResults$weightedPSM
  }
}
save(ARI_PY,ARI_DPM,FMI_PY,FMI_DPM,NMI_PY,NMI_DPM,PSMs_PY_10,PSMs_PY_20,PSMs_PY_30,file="ClustLarge.rda")
load("ClustLarge.rda")

xtable(ARI_PY,digits=2,align=rep("c",11))
xtable(ARI_DPM,digits=2,align=rep("c",11))

for (j in 1:10)
{
  #save PSMs ordered by true cluster labels
  A <- readMat(sprintf("simLarge%d.mat",j))
  trueClust <- A$z[1,]
  xx <- order(trueClust)
  write.table(PSMs_PY_10[[j]][xx,xx],file=sprintf("PSMs_PY_10_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY_20[[j]][xx,xx],file=sprintf("PSMs_PY_20_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY_30[[j]][xx,xx],file=sprintf("PSMs_PY_30_%d.csv",j),sep=",",col.names = F,row.names = F)
}


#computing PSMs for 10 subsampled chains per capture time, but different numbers of chains
ARI_PY <- matrix(0,6,10)
NMI_PY <- matrix(0,6,10)
FMI_PY <- matrix(0,6,10)
ARI_DPM <- matrix(0,6,10)
NMI_DPM <- matrix(0,6,10)
FMI_DPM <- matrix(0,6,10)

xx = c(2,4,8,12,48,96)
PSMs_PY2 <- list()
PSMs_PY4 <- list()
PSMs_PY8 <- list()
PSMs_PY12 <- list()
PSMs_PY48 <- list()
PSMs_PY96 <- list()
for (j in 1:10){
  A <- readMat(sprintf("simLarge%d.mat",j))
  trueClust <- A$z[1,]
  PSMS <- readMat(sprintf("PSMsLarge%d_%d.mat",10,j))$PSMs
  PSMs_PY <- list()
  for (kk in 1:6){
    k = xx[kk]
    PSMs <- PSMS[,,1:k]
    summaryResults <- processPSMs(PSMs)
    ARI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "adjusted_rand_index")
    ARI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "adjusted_rand_index")
    FMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "fowlkes_mallows_index")
    FMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "fowlkes_mallows_index")
    NMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "nmi")
    NMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "nmi")
    if (kk==1)
      PSMs_PY2[[j]] <- summaryResults$weightedPSM
    else if (kk==2)
      PSMs_PY4[[j]] <- summaryResults$weightedPSM
    else if (kk==3)
      PSMs_PY8[[j]] <- summaryResults$weightedPSM
    else if (kk==4)
      PSMs_PY12[[j]] <- summaryResults$weightedPSM
    else if (kk==5)
      PSMs_PY48[[j]] <- summaryResults$weightedPSM
    else
      PSMs_PY96[[j]] <- summaryResults$weightedPSM
  }
}
save(ARI_PY,ARI_DPM,FMI_PY,FMI_DPM,NMI_PY,NMI_DPM,PSMs_PY2,PSMs_PY4,PSMs_PY8,PSMs_PY12,
     PSMs_PY48,PSMs_PY96,file="ClustLarge10.rda")
load("ClustLarge10.rda")
xtable(ARI_PY,digits=2,align=rep("c",11))
xtable(ARI_DPM,digits=2,align=rep("c",11))


for (j in 1:10)
{
  #save PSMs ordered by true cluster labels
  A <- readMat(sprintf("simLarge%d.mat",j))
  trueClust <- A$z[1,]
  xx <- order(trueClust)
  write.table(PSMs_PY2[[j]][xx,xx],file=sprintf("PSMs_PY2_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY4[[j]][xx,xx],file=sprintf("PSMs_PY4_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY8[[j]][xx,xx],file=sprintf("PSMs_PY8_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY12[[j]][xx,xx],file=sprintf("PSMs_PY12_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY48[[j]][xx,xx],file=sprintf("PSMs_PY48_%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY96[[j]][xx,xx],file=sprintf("PSMs_PY96_%d.csv",j),sep=",",col.names = F,row.names = F)
}

#larger data sets with more noise
ARI_PY <- matrix(0,3,5)
NMI_PY <- matrix(0,3,5)
FMI_PY <- matrix(0,3,5)
ARI_DPM <- matrix(0,3,5)
NMI_DPM <- matrix(0,3,5)
FMI_DPM <- matrix(0,3,5)
PSMs_PY_10 <- list()
PSMs_PY_20 <- list()
PSMs_PY_30 <- list()
#computing PSMs for simulation with large number of cells
#with 10, 20, and finally 30 subsampled cells per chain and capture time
for (j in 1:5){
  A <- readMat(sprintf("simLargeMoreNoise%d.mat",j))
  trueClust <- A$z[1,]
  for (kk in 1:3){
    k = kk*10
    PSMs <- readMat(sprintf("PSMsLargeMoreNoise%d_%d.mat",k,j))$PSMs

    summaryResults <- processPSMs(PSMs)
    ARI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "adjusted_rand_index")
    ARI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "adjusted_rand_index")
    FMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "fowlkes_mallows_index")
    FMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "fowlkes_mallows_index")
    NMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "nmi")
    NMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "nmi")
    #also save the PSM
    if (kk==1)
      PSMs_PY_10[[j]] <- summaryResults$weightedPSM
    else if (kk==2)
      PSMs_PY_20[[j]] <- summaryResults$weightedPSM
    else
      PSMs_PY_30[[j]] <- summaryResults$weightedPSM
  }
}
save(ARI_PY,ARI_DPM,FMI_PY,FMI_DPM,NMI_PY,NMI_DPM,PSMs_PY_10,PSMs_PY_20,PSMs_PY_30,file="ClustLargeMoreNoise.rda")
load("ClustLargeMoreNoise.rda")

xtable(ARI_PY,digits=2,align=rep("c",6))
xtable(ARI_DPM,digits=2,align=rep("c",6))

for (j in 1:5)
{
  #save PSMs ordered by true cluster labels
  A <- readMat(sprintf("simLargeMoreNoise%d.mat",j))
  trueClust <- A$z[1,]
  xx <- order(trueClust)
  write.table(PSMs_PY_10[[j]][xx,xx],file=sprintf("PSMs_PY_10_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY_20[[j]][xx,xx],file=sprintf("PSMs_PY_20_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY_30[[j]][xx,xx],file=sprintf("PSMs_PY_30_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
}


#computing PSMs for 10 subsampled chains per capture time, but different numbers of chains
ARI_PY <- matrix(0,6,5)
NMI_PY <- matrix(0,6,5)
FMI_PY <- matrix(0,6,5)
ARI_DPM <- matrix(0,6,5)
NMI_DPM <- matrix(0,6,5)
FMI_DPM <- matrix(0,6,5)

xx = c(2,4,8,12,48,96)
PSMs_PY2 <- list()
PSMs_PY4 <- list()
PSMs_PY8 <- list()
PSMs_PY12 <- list()
PSMs_PY48 <- list()
PSMs_PY96 <- list()
for (j in 1:5){
  A <- readMat(sprintf("simLargeMoreNoise%d.mat",j))
  trueClust <- A$z[1,]
  PSMS <- readMat(sprintf("PSMsLargeMoreNoise%d_%d.mat",10,j))$PSMs
  PSMs_PY <- list()
  for (kk in 1:6){
    k = xx[kk]
    PSMs <- PSMS[,,1:k]
    summaryResults <- processPSMs(PSMs)
    ARI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "adjusted_rand_index")
    ARI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "adjusted_rand_index")
    FMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "fowlkes_mallows_index")
    FMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "fowlkes_mallows_index")
    NMI_PY[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR,
                                       method = "nmi")
    NMI_DPM[kk,j] = external_validation(trueClust,summaryResults$sumClustPEAR_DP,
                                        method = "nmi")
    if (kk==1)
      PSMs_PY2[[j]] <- summaryResults$weightedPSM
    else if (kk==2)
      PSMs_PY4[[j]] <- summaryResults$weightedPSM
    else if (kk==3)
      PSMs_PY8[[j]] <- summaryResults$weightedPSM
    else if (kk==4)
      PSMs_PY12[[j]] <- summaryResults$weightedPSM
    else if (kk==5)
      PSMs_PY48[[j]] <- summaryResults$weightedPSM
    else
      PSMs_PY96[[j]] <- summaryResults$weightedPSM
  }
}
save(ARI_PY,ARI_DPM,FMI_PY,FMI_DPM,NMI_PY,NMI_DPM,PSMs_PY2,PSMs_PY4,PSMs_PY8,PSMs_PY12,
     PSMs_PY48,PSMs_PY96,file="ClustLargeMoreNoise10.rda")
load("ClustLargeMoreNoise10.rda")
xtable(ARI_PY,digits=2,align=rep("c",6))
xtable(ARI_DPM,digits=2,align=rep("c",6))


for (j in 1:5)
{
  #save PSMs ordered by true cluster labels
  A <- readMat(sprintf("simLargeMoreNoise%d.mat",j))
  trueClust <- A$z[1,]
  xx <- order(trueClust)
  write.table(PSMs_PY2[[j]][xx,xx],file=sprintf("PSMs_PY2_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY4[[j]][xx,xx],file=sprintf("PSMs_PY4_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY8[[j]][xx,xx],file=sprintf("PSMs_PY8_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY12[[j]][xx,xx],file=sprintf("PSMs_PY12_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY48[[j]][xx,xx],file=sprintf("PSMs_PY48_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
  write.table(PSMs_PY96[[j]][xx,xx],file=sprintf("PSMs_PY96_MoreNoise%d.csv",j),sep=",",col.names = F,row.names = F)
}

#convergence check

convLarge <- list()
for (j in 1:10){
  A <- readMat(sprintf("simLarge%d.mat",j))
  trueClust <- A$z[1,]
  for (kk in 1:1){
    k = kk*10
    PSMs <- readMat(sprintf("PSMsLarge%d_%d.mat",k,j))$PSMs[,,1:12]
    convLarge[[j]] <- checkConvergence(PSMs)
  }}
save(convLarge,file="convSimLarge.rda")

convLargeMoreNoise <- list()
for (j in 1:5){
  A <- readMat(sprintf("simLargeMoreNoise%d.mat",j))
  trueClust <- A$z[1,]
  for (kk in 1:3){
    k = kk*10
    PSMs <- readMat(sprintf("PSMsLargeMoreNoise%d_%d.mat",k,j))$PSMs[,,1:12]
    convLargeMoreNoise[[j]] <- checkConvergence(PSMs)
  }}

save(convLargeMoreNoise,file='convSimLargeMoreNoise.rda')




#compute Rhat statistic for alpha across subsampled chains
RhatStats12Chains <- matrix(nrow=10,ncol=10)
for (j in 1:10){
  alphas = as.matrix(read.table(sprintf("alphasSimLarge%d.csv",j),sep=","))
  RhatStats12Chains[,j] <- RhatConcentration(alphas,12)
}

RhatStats4Chains <- matrix(nrow=10,ncol=10)
for (j in 1:10){
  alphas = as.matrix(read.table(sprintf("alphasSimLarge%d.csv",j),sep=","))
  RhatStats4Chains[,j] <- RhatConcentration(alphas,4)
}

RhatStats8Chains <- matrix(nrow=10,ncol=10)
for (j in 1:10){
  alphas = as.matrix(read.table(sprintf("alphasSimLarge%d.csv",j),sep=","))
  RhatStats8Chains[,j] <- RhatConcentration(alphas,8)
}

RhatStats24Chains <- matrix(nrow=10,ncol=10)
for (j in 1:10){
  alphas = as.matrix(read.table(sprintf("alphasSimLarge%d.csv",j),sep=","))
  RhatStats24Chains[,j] <- RhatConcentration(alphas,24)
}

RhatStats48Chains <- matrix(nrow=10,ncol=10)
for (j in 1:10){
  alphas = as.matrix(read.table(sprintf("alphasSimLarge%d.csv",j),sep=","))
  RhatStats48Chains[,j] <- RhatConcentration(alphas,48)
}


RhatStats12ChainsMoreNoise <- matrix(nrow=10,ncol=5)
for (j in 1:5){
alphas = as.matrix(read.table(sprintf("alphasSimLargeMoreNoise%d.csv",j),sep=","))
RhatStats12ChainsMoreNoise[,j] <- RhatConcentration(alphas,12)
}

RhatStats4ChainsMoreNoise <- matrix(nrow=10,ncol=5)
for (j in 1:5){
  alphas = as.matrix(read.table(sprintf("alphasSimLargeMoreNoise%d.csv",j),sep=","))
  RhatStats4ChainsMoreNoise[,j] <- RhatConcentration(alphas,4)
}

RhatStats8ChainsMoreNoise <- matrix(nrow=10,ncol=5)
for (j in 1:5){
  alphas = as.matrix(read.table(sprintf("alphasSimLargeMoreNoise%d.csv",j),sep=","))
  RhatStats8ChainsMoreNoise[,j] <- RhatConcentration(alphas,8)
}

RhatStats24ChainsMoreNoise <- matrix(nrow=10,ncol=5)
for (j in 1:5){
  alphas = as.matrix(read.table(sprintf("alphasSimLargeMoreNoise%d.csv",j),sep=","))
  RhatStats24ChainsMoreNoise[,j] <- RhatConcentration(alphas,24)
}


#figures
library(ggplot2)
load("convSimLarge.rda")
dist_Frobenius <- rep(0,100)
for (j in 1:10){
dist_Frobenius[((j-1)*10+1):(j*10)] <- convLarge[[j]]$distPY/nrow(convLarge[[1]]$results[[1]]$weightedPSM)^2
}
dataSet <- c()
number_Chains <- c()
for (j in 1:10)
{
  dataSet <- c(dataSet, rep(toString(j),10))
  number_Chains <- c(number_Chains,3:12)
}
dist_FrobeniusDF <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(dist_FrobeniusDF,aes(x=number_Chains,y=dist_Frobenius)) + geom_point(alpha=0.3,color="blue",size=5)+
  scale_y_continuous(name = "Frobenius distance/# matrix elements")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains",limits=c(3, 12),breaks=c(3,seq(4,12,1)))+theme(legend.position="top")+theme_classic()
p

cophPY <- rep(0,110)
for (j in 1:10){
cophPY[((j-1)*11+1):(j*11)] <- convLarge[[j]]$coph[[length(convLarge[[j]]$coph)]]
}
dataSet <- c()
number_Chains <- c()
for (j in 1:10)
{
  dataSet <- c(dataSet, rep(toString(j),11))
  number_Chains <- c(number_Chains,2:12)
}
number_Chains <- seq(2,12,1)
cophDF <- data.frame(coph=cophPY,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=dataSet)) + geom_point(aes(color=dataSet))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 12),breaks=c(2,seq(3,12,1)))+scale_color_brewer(palette="Spectral")+
  guides(color=guide_legend(title="data set"))+theme(legend.position="top")
p

cophPYRel <- rep(0,110)
for (j in 1:10){
  a <- convLarge[[j]]$coph[[length(convLarge[[j]]$coph)]]
  cophPYRel[((j-1)*11+1):(j*11)] <- a/a[11]
}
cophDF <- data.frame(coph=cophPYRel,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=dataSet)) + geom_point(aes(color=dataSet))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 12),breaks=c(2,seq(3,12,1)))+theme(legend.position="top")+scale_color_brewer(palette="Spectral")+
  guides(color=guide_legend(title="data set"))
p


nChains <- c(rep(4,100),rep(8,100),rep(12,100),rep(24,100),rep(48,100))
dataSet <- c()
for (j in 1:10)
{
  dataSet <- c(dataSet, rep(toString(j),10))
} 
dataSetInd <- factor(rep(dataSet,5),levels=c("1","2","3","4","5","6","7","8","9","10"))
convSimAlpha <- c(as.vector(RhatStats4Chains),as.vector(RhatStats8Chains),as.vector(RhatStats12Chains),as.vector(RhatStats24Chains),as.vector(RhatStats48Chains))
group <- factor(c(rep("4",100),rep("8",100),rep("12",100),rep("24",100),rep("48",100)),
                levels = c("4","8","12","24","48"))
convSimAlphaDF <- data.frame(convSimAlpha=convSimAlpha,nChains=nChains,group=group,dataSetInd=dataSetInd)
p <- ggplot(convSimAlphaDF,aes(x = group,y=convSimAlpha, group=nChains,color=dataSetInd)) +geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_color_brewer(palette="Spectral")+guides(color=guide_legend(title="data set"))+theme(legend.position="top")
p

#for simulated data sets with higher noise levels
load("convSimLargeMoreNoise.rda")
dist_Frobenius <- rep(0,50)
for (j in 1:5){
  dist_Frobenius[((j-1)*10+1):(j*10)] <- convLargeMoreNoise[[j]]$distPY/nrow(convLargeMoreNoise[[1]]$results[[1]]$weightedPSM)^2
}
dataSet <- c()
number_Chains <- c()
for (j in 1:5)
{
  dataSet <- c(dataSet, rep(toString(j),10))
  number_Chains <- c(number_Chains,3:12)
}
dist_FrobeniusDF <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(dist_FrobeniusDF,aes(x=number_Chains,y=dist_Frobenius)) + geom_point(alpha=0.3,color="blue",size=5)+
  scale_y_continuous(name = "Frobenius distance/# matrix elements")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains",limits=c(3, 12),breaks=c(3,seq(4,12,1)))+theme(legend.position="top")+theme_classic()
p

cophPY <- rep(0,55)
for (j in 1:10){
  cophPY[((j-1)*11+1):(j*11)] <- convLargeMoreNoise[[j]]$coph[[length(convLargeMoreNoise[[j]]$coph)]]
}
dataSet <- c()
number_Chains <- c()
for (j in 1:10)
{
  dataSet <- c(dataSet, rep(toString(j),11))
  number_Chains <- c(number_Chains,2:12)
}
number_Chains <- seq(2,12,1)
cophDF <- data.frame(coph=cophPY,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=dataSet)) + geom_point(aes(color=dataSet))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 12),breaks=c(2,seq(3,12,1)))+scale_color_brewer(palette="Spectral")+
  guides(color=guide_legend(title="data set"))+theme(legend.position="top")
p

cophPYRel <- rep(0,55)
for (j in 1:5){
  a <- convLargeMoreNoise[[j]]$coph[[length(convLargeMoreNoise[[j]]$coph)]]
  cophPYRel[((j-1)*11+1):(j*11)] <- a/a[11]
}
cophDF <- data.frame(coph=cophPYRel,number_Chains=number_Chains,dataSet=dataSet)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=dataSet)) + geom_point(aes(color=dataSet))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 12),breaks=c(2,seq(3,12,1)))+theme(legend.position="top")+scale_color_brewer(palette="Spectral")+
  guides(color=guide_legend(title="data set"))
p


nChains <- c(rep(4,50),rep(8,50),rep(12,50),rep(24,50))
dataSet <- c()
for (j in 1:5)
{
  dataSet <- c(dataSet, rep(toString(j),10))
} 
dataSetInd <- factor(rep(dataSet,4),levels=c("1","2","3","4","5","6","7","8","9","10"))
convSimAlpha <- c(as.vector(RhatStats4ChainsMoreNoise),as.vector(RhatStats8ChainsMoreNoise),as.vector(RhatStats12ChainsMoreNoise),as.vector(RhatStats24ChainsMoreNoise))
group <- factor(c(rep("4",50),rep("8",50),rep("12",50),rep("24",50)),
                levels = c("4","8","12","24"))
convSimAlphaDF <- data.frame(convSimAlpha=convSimAlpha,nChains=nChains,group=group,dataSetInd=dataSetInd)
p <- ggplot(convSimAlphaDF,aes(x = group,y=convSimAlpha, group=nChains,color=dataSetInd)) +geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_color_brewer(palette="Spectral")+guides(color=guide_legend(title="data set"))+theme(legend.position="top")
p
