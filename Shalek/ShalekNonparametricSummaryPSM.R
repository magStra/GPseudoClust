library(nonparametricSummaryPSM)
library(R.matlab)
library(mcclust)

PSMs_Shalek <- readMat("PSMsShalek.mat")
PSMs_Shalek <- PSMs_Shalek$PSMsShalek
npSummaryPSMShalek <- processPSMs(PSMs_Shalek)
writeMat('npSummaryPSMShalek.mat',npSummaryPSMShalek=npSummaryPSMShalek)
#Shalek data, first 24 chains
PSMs_Shalek24 <- PSMs_Shalek[,,1:24]
npSummaryPSMShalek24 <- processPSMs(PSMs_Shalek24)
writeMat('npSummaryPSMShalek24.mat',npSummaryPSMShalek24=npSummaryPSMShalek24)
#Shalek data, first 4 chains
PSMs_Shalek4 <- PSMs_Shalek[,,1:4]
npSummaryPSMShalek4 <- processPSMs(PSMs_Shalek4)
writeMat('npSummaryPSMShalek4.mat',npSummaryPSMShalek4=npSummaryPSMShalek4)

PSMs_Shalek <- readMat("PSMsShalek.mat")
convShalek <- checkConvergence(PSMs_Shalek$PSMsShalek)
save(convShalek,file="convShalek.rda")
alphasShalek <- as.matrix(read.table("alphasShalek.csv",sep=","))
convShalekAlpha4 <- RhatConcentration(alphasShalek,4)
convShalekAlpha8 <- RhatConcentration(alphasShalek,8)
convShalekAlpha12 <- RhatConcentration(alphasShalek,12)
convShalekAlpha24 <- RhatConcentration(alphasShalek,24)
convShalekAlpha32 <- RhatConcentration(alphasShalek,32)
convShalekAlpha48 <- RhatConcentration(alphasShalek,48)


#plots
library(ggplot2)
load("convShalek.rda")
dist_Frobenius <- c(convShalek$distPY,convShalek$distDP)/nrow(convShalek$results[[1]]$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convShalek$distPY)),rep("PY+PEAR",length(convShalek$distPY)))
number_Chains <- seq(3,96,1)
dist_Frobenius <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_Frobenius,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements",breaks = c(0.0001,0.0002,0.0003,0.0004,0.0005,0.001,0.0015))+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(3, 96),breaks=c(3,seq(12,96,12)))+theme(legend.position="top")
p

cophPY <- convShalek$coph[[length(convShalek$coph)]]
cophDPM <- convShalek$coph_DP[[length(convShalek$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p

cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,96,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains", limits=c(2, 96),breaks=c(2,seq(12,96,12)))+theme(legend.position="top")
p

nChains <- c(rep(4,10),rep(8,10),rep(12,10),rep(24,10),rep(32,10),rep(48,10))
convShalekAlpha <- c(convShalekAlpha4,convShalekAlpha8,convShalekAlpha12,convShalekAlpha24,convShalekAlpha32,convShalekAlpha48)
group <- factor(c(rep("4",10),rep("8",10),rep("12",10),rep("24",10),rep("32",10),rep("48",10)),
                levels = c("4","8","12","24","32","48"))
convShalekAlphaDF <- data.frame(convShalekAlpha=convShalekAlpha,nChains=nChains,group=group)
p <- ggplot(convShalekAlphaDF,aes(x = group,y=convShalekAlpha, group=nChains)) +geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("number of chains") + ylab("extended GR statistic")+theme(text=element_text(size=12),axis.text=element_text(size=12))
p

